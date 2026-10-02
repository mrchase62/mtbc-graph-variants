#!/usr/bin/env python3
"""A projection store that persists between runs, keyed by position.

`odgi position` answers a question with no isolate and no cohort in it: where
does H37Rv position p land in reference R, and at what distance. Caching that
inside a run took gwas1000's p5states from about 1,600 CPU-hours to 240.
Keeping it BETWEEN runs removes it from the critical path for every reference
already seen.

KEYED BY POSITION, NOT BY KEY SET. The within-run cache was keyed on a
checksum of the key column, so gwas1000 got no benefit from the 150 references
scale200 had already projected -- a different cohort has a different key set
and therefore a different checksum, even though the positions overlap almost
entirely. This store is keyed on (reference, position), so a new cohort only
pays for the positions nobody has projected yet.

APPEND-ONLY PARTS, so fifty concurrent tasks cannot corrupt it. Each writer
drops a new file into <store>/<REFID>/ under a temporary name and renames it;
readers concatenate every part and take the first value for each position. No
writer ever rewrites another's file, so a killed task can leave at worst an
unnamed temporary.

BUILD-SCOPED. The store lives under the build directory, so a panel rebuild
gets a new store rather than being served coordinates from the old graph.

    missing  positions in --need that the store does not cover, for odgi
    add      fold a projection result into the store
    emit     the store's rows for the positions in --need
"""
import argparse, os, sys, time


def src_of(line):
    """The source position an odgi/frame_convert row reports, 1-based."""
    try:
        return int(line.split("\t", 1)[0].rsplit(",", 2)[-2]) + 1
    except (ValueError, IndexError):
        return None


def load(store, ref):
    d = os.path.join(store, ref)
    out = {}
    if not os.path.isdir(d):
        return out
    for f in sorted(os.listdir(d)):
        if not f.endswith(".pos"):
            continue
        for line in open(os.path.join(d, f)):
            if line.startswith("#"):
                continue
            p = src_of(line)
            if p is not None:
                out.setdefault(p, line)
    return out


def need_positions(path):
    """The anchor file odgi is given: one `path,offset,strand` per line."""
    out = []
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        try:
            out.append((int(line.rsplit(",", 2)[-2]) + 1, line))
        except (ValueError, IndexError):
            pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("missing", "add", "emit"))
    ap.add_argument("--store", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--need", default="")
    ap.add_argument("--result", default="")
    ap.add_argument("--out", default="")
    ap.add_argument("--allow-missing", action="store_true",
                    help="emit what the store has instead of refusing, for a "
                         "caller that checks its own completeness floor")
    a = ap.parse_args()

    if a.mode == "add":
        d = os.path.join(a.store, a.ref)
        os.makedirs(d, exist_ok=True)
        stem = f"{os.getpid()}.{int(time.time()*1000)}"
        tmp = os.path.join(d, f".{stem}.tmp")
        n = 0
        with open(tmp, "w") as o:
            for line in open(a.result):
                if line.startswith("#") or src_of(line) is None:
                    continue
                o.write(line)
                n += 1
        # Build the final name from the stem, not by string surgery on the
        # whole path: `.replace("/.", "/")` also rewrote any "/./" or "/../"
        # in the STORE path, so a store given as a relative path renamed
        # into a directory that does not exist.
        os.rename(tmp, os.path.join(d, f"{stem}.pos"))
        print(f"  stored {n:,} projected positions for {a.ref}")
        return 0

    have = load(a.store, a.ref)
    need = need_positions(a.need)
    if a.mode == "missing":
        miss = [ln for p, ln in need if p not in have]
        with open(a.out or "/dev/stdout", "w") as o:
            for ln in miss:
                o.write(ln + "\n")
        print(f"  {a.ref}: {len(need) - len(miss):,} of {len(need):,} positions "
              f"already in the store, {len(miss):,} to project",
              file=sys.stderr)
        return 0

    miss = [p for p, _ in need if p not in have]
    if miss and not a.allow_missing:
        sys.exit(f"FATAL: {a.ref}: {len(miss):,} needed positions are not in "
                 f"the store; run `missing` and `add` first, or pass "
                 f"--allow-missing if the caller checks its own floor")
    # SOME POSITIONS NEVER PROJECT, and demanding them is a deadlock. odgi drops
    # the roughly 0.19% of positions that have no equivalent in a given path, so
    # `missing` keeps asking for them, odgi keeps returning nothing, and there
    # is nothing to `add`. GCF_965121955 sat at 39,620 of 39,685 stored and
    # failed 12 tasks that way, once per run, forever. With --allow-missing the
    # store emits what it has and the caller decides; p5_svgt.sh already
    # requires 95% of the probe set, which is the question that matters, and
    # 39,620 of 39,685 is 99.8%.
    with open(a.out or "/dev/stdout", "w") as o:
        for p, _ in need:
            if p in have:
                o.write(have[p])
    n = len(need) - len(miss)
    print(f"  {a.ref}: emitted {n:,} of {len(need):,} positions from the store"
          + (f", {len(miss):,} with no equivalent in this path" if miss else ""),
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
