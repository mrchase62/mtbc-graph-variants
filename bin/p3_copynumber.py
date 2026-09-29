#!/usr/bin/env python3
"""Copy-number genotyping by depth, for the 674 copy-number accessory loci.

These loci are present in H37Rv (HSP union coverage > 0.90); what varies between
genomes is how many copies. Unplaced reads say nothing about them, because the
reads place perfectly well -- they simply pile onto the copies the reference has.
The signal is depth: a sample carrying three copies where its reference carries
one shows roughly three times its median depth across that locus.

IS6110 is 501 of the 674, and reliable IS6110 detection is a standing requirement
for this project, so this is the arm that matters most and the one that has never
been run.

Two things this refuses to do:

  - report a locus whose anchor did not project into the reference's coordinates
    as "zero copies". `odgi` reports dist.to.ref != 0 when a position has no
    equivalent on the target path, and that means "cannot say", not "absent".
    Conflating the two is the encoding hazard T6 measured at 508,039 cells.
  - normalise by a mean. Median depth is used throughout, because a genome
    carrying a high-copy element has a mean inflated by the very thing being
    measured.
"""
import argparse, csv, statistics, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loci", required=True)
    ap.add_argument("--positions", required=True, help="odgi position output")
    ap.add_argument("--depth", required=True, help="samtools depth -aa on the R BAM")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--flank", type=int, default=0,
                    help="extra bp either side of the locus window")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    loci = [r for r in csv.DictReader(open(a.loci), delimiter="\t")
            if r.get("novelty") == "copy_number"]
    if not loci:
        print("no copy_number loci in the catalogue", file=sys.stderr); return 1

    # Keyed by SOURCE position, not file order. `odgi position -t 4` returns
    # results in thread-completion order: feeding 1000..8000 returns
    # 4000 5000 6000 7000 8000 1000 2000 3000. Pairing by order mispairs every
    # locus with a different locus's coordinate, and because the counts still
    # match, a count check cannot detect it.
    by_src = {}
    for line in open(a.positions):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        try:
            src = int(f[0].rsplit(",", 2)[-2]) + 1
            rpos = int(f[1].rsplit(",", 2)[-2]) + 1
            dist = int(f[2])
        except (ValueError, IndexError):
            continue
        strand = f[3].strip() if len(f) > 3 and f[3].strip() in "+-" else "+"
        by_src.setdefault(src, (rpos, dist, strand))
    proj = [by_src.get(int(l["pos"])) for l in loci]
    nmiss = sum(1 for x in proj if x is None)
    if nmiss > len(loci) * 0.5:
        print(f"  WARNING: {nmiss} of {len(loci)} anchors have no projection; "
              f"odgi output does not cover the input -- refusing", file=sys.stderr)
        return 1
    proj = [x if x is not None else (0, 1, "+") for x in proj]

    depth = {}
    for line in open(a.depth):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        depth[int(f[1])] = int(f[2])
    if not depth:
        print("empty depth file", file=sys.stderr); return 1
    genome_median = statistics.median(depth.values())
    if genome_median <= 0:
        print("genome median depth is zero", file=sys.stderr); return 1

    rows, unprojected = [], 0
    for locus, (rpos, dist, strand) in zip(loci, proj):
        L = int(locus["rep_len"])
        if dist != 0:
            unprojected += 1
            rows.append(dict(sample=a.sample, reference=a.reference,
                             locus=locus["locus_id"], h37rv_anchor=locus["pos"],
                             rep_len=L, r_pos="", window_bp=0,
                             median_depth="", depth_ratio="",
                             copies="", call="no_projection"))
            continue
        # The locus follows its anchor in H37Rv's direction. Where R runs
        # reverse to H37Rv here (strand `-`), it lies on the OTHER side of the
        # projected anchor in R, and a window read forward from rpos measured
        # unrelated flank. The window is L bases plus the flank on each side.
        # (The forward window is kept exactly as before, so forward-strand
        # results do not move; the mirror image is used for `-`.)
        if strand == "-":
            lo, hi = max(1, rpos - L - a.flank), rpos + a.flank
        else:
            lo, hi = max(1, rpos - a.flank), rpos + L + a.flank
        vals = [depth.get(p, 0) for p in range(lo, hi + 1)]
        med = statistics.median(vals) if vals else 0
        ratio = med / genome_median
        # Deliberately coarse. A ratio is not a copy number until it has been
        # calibrated against loci of known copy number, which has not been done,
        # so the call is a band and the ratio is reported for calibration later.
        if ratio < 0.15:
            call = "absent"
        elif ratio < 1.5:
            call = "single_or_reference"
        elif ratio < 2.5:
            call = "gain_2x"
        else:
            call = "gain_multi"
        rows.append(dict(sample=a.sample, reference=a.reference,
                         locus=locus["locus_id"], h37rv_anchor=locus["pos"],
                         rep_len=L, r_pos=rpos, window_bp=len(vals),
                         median_depth=med, depth_ratio=round(ratio, 4),
                         copies=round(ratio, 2), call=call))

    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    called = [r for r in rows if r["call"] != "no_projection"]
    print(f"  genome median depth {genome_median:.0f}x; {len(called)} of "
          f"{len(rows)} loci projected ({unprojected} without an equivalent "
          f"position in {a.reference})")
    import collections as _c
    for k, n in _c.Counter(r["call"] for r in rows).most_common():
        print(f"    {k:<22s} {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
