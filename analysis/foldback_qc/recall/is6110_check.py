#!/usr/bin/env python3
"""IS6110 junction check (user, 2026-10-12): does fold-back read removal
take away the reads that show a new IS6110 insertion? Report only. Run from
runroot with $MTB_PY after 06_is6110_sites.sbatch.

IS6110's two ends are 28 bp imperfect inverted repeats, so a read crossing
the edge of a copy could in principle look like a fold-back. In the
excluded-isolates check one IS6110 insertion (GCF_902459825, 2,248,575) was
lost after removal. Three checks, existing outputs only:

1. Labelling on a shared reference (as production would do it). For every
   true IS6110-sized insertion (truth INS/REPL, 1,340-1,370 bp; the isolate
   has the copy, its part D reference does not), the fold-back reads with a
   junction within the insertion point +-50 bp (06_is6110_sites.sbatch,
   fbrecall BAMs), against chance: the sample's fold-back reads genome-wide
   times the window's share of the genome. Other true events (DEL, other
   INS) are counted the same way for comparison.
2. Labelling on the own genome (as the removal test did). For the 12
   excluded isolates and the 9 chimeric controls, of the reads in the
   fbrecall BAM that overlap the insertion point +-50 bp, the share whose
   name is on the removed list, against the same share in 300 random 101 bp
   windows of the same BAM.
3. IS6110 copies in the own genomes (IS6110 mapped with minimap2). Fold-back
   junctions (selfcal54 sites, own genome) within +-50 bp of a copy end,
   against chance.
Plus the two insertions of GCF_902459825: reads at each point, all reads vs
fold-back pairs removed.
Writes out/is6110_check.tsv and prints a summary.
"""
import collections
import csv
import gzip
import os
import random
import subprocess

import pysam

D = "../analysis/foldback_qc"
R = f"{D}/recall"
REFS = "refbias/build/7713a8d71d8e-fix1/refs"
IS = "../analysis/detettore_eval/detettore6110/detettore6110/resources/is_targets/IS6110.fasta"
MM2 = "/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin/minimap2"
PAD = 50
CUT = 4.9


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def sites(path):
    """[(pos, reads)] of a fold-back sites file, and the total reads."""
    out = []
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out.append((int(r["pos"]), int(r["reads"])))
    return out, sum(n for _, n in out)


def near(st, lo, hi):
    return sum(n for p, n in st if lo - PAD <= p <= hi + PAD)


def genome_len(g):
    return sum(pysam.FastaFile(f"{REFS}/{g}.fasta").lengths)


def is6110_kind(t):
    return t["type"] in ("INS", "REPL") and 1340 <= int(t["alt_len"]) <= 1370


def kind(t):
    return "IS6110-sized INS" if is6110_kind(t) else f"other {t['type']}"


def overlap_reads(bam, ctg, lo, hi):
    return {r.query_name for r in bam.fetch(ctg, max(0, lo), hi)
            if not (r.is_secondary or r.is_supplementary or r.is_unmapped)}


def main():
    rate = {p["sample"]: float(p["rate"]) for p in tsv(f"{D}/selfcal/out/calib_points.tsv")}
    rm = {r["sample"]: r["reference"] for r in tsv("refbias/fbrecall/run/p1/refmap.tsv")}
    group = {s: "chimeric" if s.startswith("mar_") else "kept" if rate[s] <= CUT else "excluded" for s in rm}
    truth = tsv(f"{R}/out/truth.tsv")
    glen = {g: genome_len(g) for g in set(rm.values())}
    rows = []

    # 1. shared reference
    print("1. Fold-back reads at true events, labelled on the part D reference (+-50 bp)")
    agg = collections.defaultdict(lambda: [0, 0, 0.0, 0])  # events, at, expected, events with any
    cache = {}
    for t in truth:
        s = t["sample"]
        p = f"{R}/out/sites_fbrecall/{s}.tsv.gz"
        if s not in rm or not os.path.exists(p):
            continue
        if s not in cache:
            cache[s] = sites(p)
        st, tot = cache[s]
        lo, hi = int(t["r_start"]), int(t["r_end"])
        at = near(st, lo, hi)
        exp = tot * (hi - lo + 2 * PAD + 1) / glen[rm[s]]
        k = kind(t)
        rows.append(dict(check="shared_ref", sample=s, group=group[s], event=f"{t['type']}:{lo}:{t['alt_len'] if t['type'] != 'DEL' else t['ref_len']}",
                         kind=k, at=at, expected=round(exp, 3)))
        for key in ((group[s], k), ("all", k)):
            a = agg[key]
            a[0] += 1; a[1] += at; a[2] += exp; a[3] += at > 0
    print("group\tevent kind\tevents\tfold-back reads at events\texpected\tevents with any")
    for (g, k), a in sorted(agg.items()):
        print(f"{g}\t{k}\t{a[0]}\t{a[1]}\t{a[2]:.1f}\t{a[3]}")

    # 2. own-genome labels (the removal test)
    print("\n2. Reads at IS6110 insertions removed by own-genome labelling, against random windows")
    print("sample\tgroup\tIS6110 insertions\treads at them\tremoved\tshare\trandom windows share")
    pooled = [0, 0, 0, 0]
    for s in sorted(rm):
        if group[s] == "kept":
            continue
        nf = f"{R}/names/{s}.fb100.names" if s.startswith("ps_") else f"{D}/removal/names/{s}.fb100.names"
        if not os.path.exists(nf):
            continue
        names = set(open(nf).read().split())
        bam = pysam.AlignmentFile(f"refbias/work/fbrecall_p2/{s}.bam")
        ctg, L = bam.references[0], bam.lengths[0]
        ev = [t for t in truth if t["sample"] == s and is6110_kind(t)]
        if not ev:
            continue
        a = b = 0
        for t in ev:
            lo, hi = int(t["r_start"]), int(t["r_end"])
            rd = overlap_reads(bam, ctg, lo - PAD, hi + PAD)
            rem = len(rd & names)
            a += len(rd); b += rem
            rows.append(dict(check="own_genome_labels", sample=s, group=group[s], event=f"{t['type']}:{lo}",
                             kind="IS6110-sized INS", reads=len(rd), removed=rem))
        rng = random.Random(s)
        ra = rb = 0
        for _ in range(300):
            x = rng.randrange(PAD, L - PAD)
            rd = overlap_reads(bam, ctg, x - PAD, x + PAD)
            ra += len(rd); rb += len(rd & names)
        rows.append(dict(check="own_genome_random", sample=s, group=group[s], reads=ra, removed=rb))
        pooled[0] += a; pooled[1] += b; pooled[2] += ra; pooled[3] += rb
        print(f"{s}\t{group[s]}\t{len(ev)}\t{a}\t{b}\t{b / max(1, a):.4f}\t{rb / max(1, ra):.4f}")
    print(f"pooled\t\t\t{pooled[0]}\t{pooled[1]}\t{pooled[1] / max(1, pooled[0]):.4f}\t"
          f"{pooled[3] / max(1, pooled[2]):.4f}")

    # 3. IS6110 copies in the own genomes
    print("\n3. Fold-back reads at IS6110 copy ends in the own genome (selfcal54, +-50 bp)")
    print("sample\tcopies\tfold-back reads genome\tat copy ends\texpected")
    own = {r["sample"]: r["reference"] for r in tsv("refbias/selfcal54/run/p1/refmap.tsv")}
    ta = te = tc = 0
    for s, g in sorted(own.items()):
        p = f"{D}/selfcal/out/sites/{s}.tsv.gz"
        if not os.path.exists(p):
            continue
        paf = subprocess.run([MM2, "-c", "-N", "100", "-p", "0.5", "--secondary=yes", f"{REFS}/{g}.fasta", IS],
                             capture_output=True, text=True, check=True).stdout
        ends = set()
        for ln in paf.splitlines():
            f = ln.split("\t")
            if int(f[3]) - int(f[2]) >= 1300:  # near-full-length copy
                ends.update((int(f[7]), int(f[8])))
        st, tot = sites(p)
        at = sum(near(st, e, e) for e in ends)
        exp = tot * len(ends) * (2 * PAD + 1) / genome_len(g)
        ta += at; te += exp; tc += len(ends) // 2
        rows.append(dict(check="own_copy_ends", sample=s, kind="IS6110 copy ends", copies=len(ends) // 2,
                         at=at, expected=round(exp, 3)))
        print(f"{s}\t{len(ends) // 2}\t{tot}\t{at}\t{exp:.2f}")
    print(f"pooled\t{tc}\t\t{ta}\t{te:.1f}")

    # the lost event
    print("\nGCF_902459825: reads overlapping each IS6110 insertion point +-50 bp")
    s = "ps_GCF_902459825"
    nf = set(open(f"{R}/names/{s}.fb100.names").read().split())
    b0 = pysam.AlignmentFile(f"refbias/work/fbrecall_p2/{s}.bam")
    b1 = pysam.AlignmentFile(f"refbias/work/fbrecallx100_p2/{s}_fb100.bam")
    st, _ = cache.get(s, sites(f"{R}/out/sites_fbrecall/{s}.tsv.gz"))
    for t in truth:
        if t["sample"] == s and is6110_kind(t):
            x = int(t["r_start"])
            r0 = overlap_reads(b0, b0.references[0], x - PAD, x + PAD)
            r1 = overlap_reads(b1, b1.references[0], x - PAD, x + PAD)
            print(f"{x}\tall reads {len(r0)}\tafter removal {len(r1)}\ton removed list {len(r0 & nf)}\t"
                  f"fold-back junctions on shared ref {near(st, x, x)}")
    cols = []
    for r in rows:
        cols += [k for k in r if k not in cols]
    with open(f"{R}/out/is6110_check.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
