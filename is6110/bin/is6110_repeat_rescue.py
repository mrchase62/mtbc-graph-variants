#!/usr/bin/env python3
"""Rescue IS6110 copies inside the DR (CRISPR) array, where short reads cannot
place a junction to one repeat.

THE PROBLEM. The DR locus is 30-40 copies of one 36 bp direct repeat separated
by short spacers. A junction read whose chromosomal half falls in it aligns
equally well to many repeats, so bwa gives the supplementary alignment MAPQ 0
and places it at one copy at random. The element-side scan then sees the
junction scattered over many stacks about 36 bp apart, none with reads of
mapping quality >= 20 -- and reads_q is what promotion and reconciliation
threshold on (--min-reads-q 10). The copy is never called.

Measured on gwas1000: twelve L1.1.1 / L1.1.1.1 isolates matched to a zero-copy
reference had no IS6110 site at all, while reads covered the whole 1355 bp
element at chromosomal depth (element/chromosome depth ratio 0.9-2.4) and 50-100
reads per isolate were clipped at the element termini with supplementary
alignments into the DR array. Across the cohort 24 isolates are in this state.

WHAT THIS DOES. For each sample, the DR array of its own reference is located
in the element-free reference sequence. If the sample has element-side stacks
inside it, none of which reaches --min-reads with mapping quality, and together
they carry at least --min-reads reads, they are replaced by ONE stack:

  clean_pos     the first base of the array. The insertion point within the
                array is not resolvable from short reads; a fixed position
                means every isolate rescued against the same reference shares
                one key, rather than one key per random placement
  reads_q       all reads in the array. The ambiguity is between identical
                repeats of ONE locus -- the DR consensus occurs nowhere else in
                the genome, which is checked per reference below -- so the
                locus-level placement is not ambiguous
  el_start/end  summed, so both_el_termini, and with it the two-sided /
                one-sided evidence class, comes from the reads as usual
  span          the array's length, so the geometry reads two_sided_wide and
                never claims a target-site duplication it did not measure
  repeat_locus  "DR" on the rescued row, "" elsewhere

A sample with a DR stack that already passes is left exactly as it is, so no
existing call moves.

IN PLACE, AND IDEMPOTENT. <sample>.elstacks.tsv is what promote, reconcile and
the writer read, so the rescued table replaces it and the scan's own output is
kept as <sample>.elstacks.raw.tsv. A table without the repeat_locus column is
a fresh scan output (P1i rewrote it) and becomes the new raw; one with the
column is regenerated from its raw.

    is6110_repeat_rescue.py --dir refbias/<cohort>/p1i --refmap <refmap>
"""
import argparse, collections, csv, os, sys

DR = "GTCGTCAGACCCAAAACCCCGAGAGGGGACGGAAAC"          # M. tuberculosis CRISPR DR
COMP = str.maketrans("ACGTN", "TGCAN")


def read_seq(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def dr_hits(seq, max_mm):
    """1-based starts of DR matches, either strand, at most max_mm mismatches."""
    out = set()
    for q in (DR, DR.translate(COMP)[::-1]):
        core = q[12:24]
        i = seq.find(core)
        while i != -1:
            st = i - 12
            if 0 <= st <= len(seq) - len(q) and \
                    sum(a != b for a, b in zip(seq[st:st + len(q)], q)) <= max_mm:
                out.add(st + 1)
            i = seq.find(core, i + 1)
    return sorted(out)


def dr_locus(seq, max_mm=4, max_gap=2000):
    """(first, last) base of the DR array, or None.

    max_gap joins repeats across the reference's own excised element: an array
    with an IS6110 copy in it is split by about 1.4 kb in the original sequence
    and not at all in the element-free one. More than one cluster is refused,
    because the rescue's premise is a single locus.
    """
    h = dr_hits(seq, max_mm)
    if not h:
        return None, "no DR array"
    clusters = [[h[0]]]
    for p in h[1:]:
        if p - clusters[-1][-1] > max_gap:
            clusters.append([p])
        else:
            clusters[-1].append(p)
    big = [c for c in clusters if len(c) >= 3]
    if len(big) != 1:
        return None, f"{len(big)} DR clusters"
    c = big[0]
    return (c[0], c[-1] + len(DR) - 1), ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="the cohort's P1i directory")
    ap.add_argument("--refmap", required=True)
    ap.add_argument("--clean-dir", default="is6110/assets/isclean_matched",
                    help="element-free reference FASTAs, <ref>.isclean.fasta")
    ap.add_argument("--min-reads", type=int, default=10,
                    help="the --min-reads-q promotion and reconciliation use")
    ap.add_argument("--slop", type=int, default=100,
                    help="bases either side of the array still counted in it")
    a = ap.parse_args()

    ref_of = {r["sample"]: r["reference"] for r in
              csv.DictReader(open(a.refmap, newline=""), delimiter="\t")}
    loci, why = {}, {}
    counts = collections.Counter()
    rescued = []
    for s in sorted(ref_of):
        cur = os.path.join(a.dir, f"{s}.elstacks.tsv")
        raw = os.path.join(a.dir, f"{s}.elstacks.raw.tsv")
        if not os.path.exists(cur) and not os.path.exists(raw):
            counts["no element-side table"] += 1
            continue
        with open(cur, newline="") as fh:
            hdr = next(csv.reader(fh, delimiter="\t"), [])
        if "repeat_locus" not in hdr:
            # a fresh scan output: it is the raw table from now on
            os.replace(cur, raw)
        rows = list(csv.DictReader(open(raw, newline=""), delimiter="\t"))
        fields = (list(rows[0]) if rows else hdr) + ["repeat_locus"]
        for r in rows:
            r["repeat_locus"] = ""

        R = ref_of[s]
        if R not in loci:
            fa = os.path.join(a.clean_dir, f"{R}.isclean.fasta")
            loci[R], why[R] = (dr_locus(read_seq(fa)) if os.path.exists(fa)
                               else (None, "no element-free reference"))
        L = loci[R]
        out = rows
        if L is None:
            counts[f"reference: {why[R]}"] += 1
        else:
            lo, hi = L[0] - a.slop, L[1] + a.slop
            inside = [r for r in rows if lo <= int(r["clean_pos"]) <= hi]
            tot = sum(int(r["reads"]) for r in inside)
            if not inside:
                counts["no element evidence in the DR array"] += 1
            elif any(int(r["reads_q"]) >= a.min_reads for r in inside):
                counts["already called in the DR array; unchanged"] += 1
            elif tot < a.min_reads:
                counts[f"fewer than {a.min_reads} reads in the DR array; unchanged"] += 1
            else:
                near = min(inside, key=lambda r: abs(int(r["clean_pos"]) - L[0]))
                off = int(near["orig_pos"]) - int(near["clean_pos"])
                agg = dict(near)
                es = sum(int(r["el_start"]) for r in inside)
                ee = sum(int(r["el_end"]) for r in inside)
                cs = sum(int(r["chr_start"]) for r in inside)
                ce = sum(int(r["chr_end"]) for r in inside)
                mq = sum(float(r["sa_mapq_mean"]) * int(r["reads"]) for r in inside)
                agg.update(
                    clean_pos=L[0], orig_pos=L[0] + off, reads=tot, reads_q=tot,
                    positions=sum(int(r["positions"]) for r in inside),
                    span=L[1] - L[0] + 1,
                    sa_mapq_max=max(int(r["sa_mapq_max"]) for r in inside),
                    sa_mapq_mean=round(mq / tot, 1),
                    el_start=es, el_end=ee,
                    el_internal=sum(int(r["el_internal"]) for r in inside),
                    chr_start=cs, chr_end=ce,
                    both_el_termini=int(es > 0 and ee > 0),
                    both_chr_sides=int(cs > 0 and ce > 0),
                    fwd=sum(int(r["fwd"]) for r in inside),
                    rev=sum(int(r["rev"]) for r in inside),
                    repeat_locus="DR")
                out = [r for r in rows if r not in inside] + [agg]
                out.sort(key=lambda r: int(r["clean_pos"]))
                kind = "two-sided" if agg["both_el_termini"] else "one-sided"
                counts[f"RESCUED, {kind}"] += 1
                rescued.append((s, R, tot, len(inside), kind))

        tmp = cur + ".tmp"
        with open(tmp, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                               lineterminator="\n")
            w.writeheader(); w.writerows(out)
        os.replace(tmp, cur)

    print(f"  DR-array rescue over {len(ref_of)} samples, {len(loci)} references")
    for k, v in counts.most_common():
        print(f"    {v:5d}  {k}")
    for s, R, tot, n, kind in rescued:
        print(f"    rescued {s} ({R}): {tot} reads from {n} stacks, {kind}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
