#!/usr/bin/env python3
"""Stage 2: earn REF for the samples that did not report an IS6110 site.

Stage 1 (`is6110_p5_merge.py`) leaves every non-carrier NOCALL, which asserts
nothing but leaves 94.7% of the IS6110 block uninformative. This asks each
non-carrier's own element-free alignment whether it looked and found nothing.

    for each key, take a carrier's own (reference, position) as the source
    project it onto the non-carrier's reference path            -> ABSENT if
                                                                   dist != 0
    convert that reference position into the non-carrier's clean coordinate
    read depth there, and ask whether a junction candidate sits within --window

    depth >= --min-dp and no candidate   REF      it looked and found nothing
    a junction candidate is present      reported, NOT silently called ALT
    otherwise                            NOCALL

WHY THE PROJECTION GOES CARRIER -> TARGET AND NOT VIA H37Rv. 68 of the 170 keys
are off the H37Rv path, so a route through H37Rv cannot carry them -- those are
exactly the sites in sequence H37Rv does not have. Projecting straight from a
carrier's path to the target's path keeps them. It is batched by TARGET path,
one `odgi position -r <target>` per distinct reference with every key's source
in a single input file, so this is ~18 graph loads rather than one per pair.

FRAMES. The carrier position comes from a refs-frame alignment and `odgi` reads
panel coordinates; the returned target position is panel and the crossmap is
refs. Both conversions go through graphframe, and for a target the panel stores
reverse complemented the returned coordinate is converted with the strand, which
graph_frame.to_refs does. See graphframe/docs/GRAPH_FRAME_RESOLUTION.md.

DEPTH FLOOR. `--min-dp` defaults to 5, the same floor `p5_states.py` applies to
a SNP REF call, so an IS6110 REF means the same strength of evidence as every
other REF in the matrix.
"""
import argparse, bisect, collections, csv, importlib.util, os, subprocess, sys

def _load(name, path):
    s = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m

graph_frame = _load("graph_frame", os.path.join("graphframe", "bin", "graph_frame.py"))


def read_crossmap(path):
    """(orig_end, cum_deleted) rows, so an original coordinate converts to the
    element-free one by subtracting the deletions that precede it."""
    out = []
    for r in csv.DictReader(open(path, newline=""), delimiter="\t"):
        out.append((int(r["orig_end"]), int(r["cum_deleted"])))
    out.sort()
    return out


def orig_to_clean(cm, pos):
    i = bisect.bisect_right([c[0] for c in cm], pos)
    return pos - (cm[i - 1][1] if i else 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort-keys", default="is6110/results/p1i_cohort_keys.tsv")
    ap.add_argument("--stage1-states", default="is6110/results/p5_is6110_states.tsv")
    ap.add_argument("--refmap", default="refbias/p1/refmap.tsv")
    ap.add_argument("--paths", default="refbias/build/7713a8d71d8e/assets/paths.txt")
    ap.add_argument("--graph", default=None)
    ap.add_argument("--odgi", default=os.environ.get("MTB_ODGI", "odgi"))
    ap.add_argument("--samtools", default=os.environ.get("MTB_SAMTOOLS", "samtools"))
    ap.add_argument("--isclean-dir", default="refbias/p1i")
    ap.add_argument("--crossmap-dir", default="is6110/assets/isclean_matched")
    ap.add_argument("--min-dp", type=int, default=5)
    ap.add_argument("--window", type=int, default=10)
    ap.add_argument("--workdir", default="refbias/p1i/p5stage2")
    ap.add_argument("--out", default="is6110/results/p5_is6110_states.stage2.tsv")
    a = ap.parse_args()

    og = a.graph
    if og is None:
        import glob
        g = sorted(glob.glob("graphs/CX333.s10k.k23.K15/*.smooth.final.og"))
        og = g[0] if g else sys.exit("no graph; pass --graph")
    os.makedirs(a.workdir, exist_ok=True)
    fr = graph_frame.Frames()
    paths = {l.split("#")[0]: l.strip() for l in open(a.paths) if l.strip()}
    ref_of = {r["sample"]: r["reference"] for r in
              csv.DictReader(open(a.refmap, newline=""), delimiter="\t")}

    # one carrier per key: its own reference and its own coordinate
    src = {}
    for x in csv.DictReader(open(a.cohort_keys, newline=""), delimiter="\t"):
        if x["frame"] == "h37rv":
            k = f'h37rv:{x["h37rv_pos"]}:'
        else:
            k = f'node:{x["node"]}:'
        src.setdefault((x["sample"], x["reference"], int(x["r_pos"])), []).append(k)

    states = list(csv.DictReader(open(a.stage1_states, newline=""), delimiter="\t"))
    # map the stage-1 key back to the same short form the carrier table uses
    def short(k):
        f = k.split(":")
        return f'{f[0]}:{f[1]}:' if f[0] == "h37rv" else f'node:{f[1]}:'

    need = collections.defaultdict(set)     # target sample -> short keys
    for r in states:
        if r["state"] == "NOCALL":
            need[r["sample"]].add(short(r["key"]))

    # by target reference, so one odgi load serves every sample sharing it
    by_target = collections.defaultdict(set)
    for s, ks in need.items():
        by_target[ref_of[s]] |= ks
    carrier_of = {}
    for (cs, cref, cpos), ks in src.items():
        for k in ks:
            carrier_of.setdefault(k, (cref, cpos))

    proj = {}          # (target ref, short key) -> (refs pos or None, dist)
    for tref in sorted(by_target):
        tpath = paths.get(tref)
        if tpath is None:
            continue
        want = sorted(by_target[tref])
        qf = os.path.join(a.workdir, f"{tref}.pos")
        order = []
        with open(qf, "w") as fh:
            for k in want:
                cref, cpos = carrier_of.get(k, (None, None))
                if cref is None or cref not in paths:
                    continue
                fh.write(f"{paths[cref]},{fr.to_panel(cref, cpos - 1)},+\n")
                order.append((k, cref, cpos))
        if not order:
            continue
        p = subprocess.run([a.odgi, "position", "-i", og, "-F", qf, "-r", tpath,
                            "-t", "4"], capture_output=True, text=True)
        got = {}
        for line in p.stdout.split("\n"):
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 3:
                continue
            sp = f[0].rsplit(",", 2)
            tp = f[1].rsplit(",", 2)
            try:
                got[(sp[0], int(sp[1]))] = (int(tp[1]), int(f[2]))
            except (ValueError, IndexError):
                continue
        for k, cref, cpos in order:
            v = got.get((paths[cref], fr.to_panel(cref, cpos - 1)))
            if v is None:
                proj[(tref, k)] = (None, -1)
            else:
                tp, dist = v
                proj[(tref, k)] = (fr.to_refs(tref, tp) + 1, dist)
        print(f"  {tref}: {len(order)} keys projected onto its path")

    # per-sample evidence: depth from the element-free alignment, and whether a
    # junction candidate already sits there
    cand = {}

    def candidates(sample):
        if sample in cand:
            return cand[sample]
        p = os.path.join(a.isclean_dir, f"{sample}.junctions.tsv")
        v = []
        if os.path.exists(p):
            for r in csv.DictReader(open(p, newline=""), delimiter="\t"):
                try:
                    v.append(int(r["pos"]))
                except (ValueError, KeyError):
                    pass
        v.sort(); cand[sample] = v
        return v

    def depths_for(sample, contig, positions):
        """Every position this sample needs, in ONE samtools call.

        This used to be one `samtools depth -r` subprocess per position. Over
        the 23-isolate pilot that was ~3,700 calls and tolerable; over
        scale100 it was ~56,000 and took 33 minutes, nearly all of it process
        startup. A BED of the sample's positions and a single `-b` call is the
        same query with one process, and the result is identical: -a reports
        every position in the regions, including the zero-depth ones, so a
        position with no coverage still comes back rather than being missing
        and silently defaulting.
        """
        out = {}
        if not positions:
            return out
        bed = os.path.join(a.workdir, f"{sample}.depth.bed")
        with open(bed, "w") as fh:
            for p in sorted(positions):
                fh.write(f"{contig}\t{p - 1}\t{p}\n")
        bam = os.path.join(a.isclean_dir, f"{sample}.isclean.bam")
        r = subprocess.run([a.samtools, "depth", "-a", "-b", bed, bam],
                           capture_output=True, text=True)
        for line in r.stdout.split("\n"):
            f = line.split("\t")
            if len(f) == 3:
                try:
                    out[int(f[1])] = int(f[2])
                except ValueError:
                    pass
        os.remove(bed)
        return {p: out.get(p, 0) for p in positions}

    contig_of, cm_of = {}, {}
    for s, R in ref_of.items():
        cmp_ = os.path.join(a.crossmap_dir, f"{R}.crossmap.tsv")
        if os.path.exists(cmp_) and R not in cm_of:
            cm_of[R] = read_crossmap(cmp_)
        bam = os.path.join(a.isclean_dir, f"{s}.isclean.bam")
        if os.path.exists(bam) and s not in contig_of:
            r = subprocess.run([a.samtools, "idxstats", bam],
                               capture_output=True, text=True)
            for line in r.stdout.split("\n"):
                f = line.split("\t")
                if len(f) > 1 and f[0].endswith("_isclean"):
                    contig_of[s] = f[0]; break

    # pass one: resolve each row to either an immediate verdict or a position
    pending = {}                       # row index -> (sample, clean position)
    verdict = {}                       # row index -> (state, evidence)
    want = collections.defaultdict(set)
    for i, r in enumerate(states):
        if r["state"] != "NOCALL":
            continue
        s = r["sample"]; R = ref_of[s]; k = short(r["key"])
        pos, dist = proj.get((R, k), (None, -1))
        if pos is None or dist != 0:
            verdict[i] = ("ABSENT", "no projection" if pos is None
                          else f"off the path, dist {dist}")
        elif R not in cm_of or s not in contig_of:
            verdict[i] = ("NOCALL", "no crossmap or no element-free alignment")
        else:
            cpos = orig_to_clean(cm_of[R], pos)
            pending[i] = (s, cpos)
            want[s].add(cpos)

    # pass two: one samtools call per sample
    depth = {}
    for s in sorted(want):
        depth[s] = depths_for(s, contig_of[s], want[s])

    # pass three: decide
    counts = collections.Counter()
    missed = []
    out = []
    for i, r in enumerate(states):
        if r["state"] != "NOCALL":
            counts[r["state"]] += 1
            out.append(dict(r, evidence="reported by this sample"))
            continue
        if i in verdict:
            st, ev = verdict[i]
        else:
            s, cpos = pending[i]
            v = candidates(s)
            j = bisect.bisect_left(v, cpos - a.window)
            near = j < len(v) and v[j] <= cpos + a.window
            dp = depth[s].get(cpos, 0)
            if near:
                st, ev = "NOCALL", f"a junction candidate sits within {a.window} bp"
                missed.append((s, r["key"], cpos, dp))
            elif dp >= a.min_dp:
                st, ev = "REF", f"depth {dp}, no junction candidate"
            else:
                st, ev = "NOCALL", f"depth {dp} below {a.min_dp}"
        counts[st] += 1
        out.append(dict(r, state=st,
                        allele=r["allele"] if st == r["state"] else "",
                        evidence=ev))

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(out)
    n = sum(counts.values())
    print(f"\n  {n} cells after stage 2")
    for k in ("ALT", "REF", "ABSENT", "NOCALL"):
        if counts[k]:
            print(f"    {k:8} {counts[k]:6} ({100*counts[k]/n:5.1f}%)")
    print(f"  {len(missed)} cells hold a junction candidate the detector did not "
          f"call a site;\n  they are left NOCALL and listed, not promoted to ALT")
    print(f"  written: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
