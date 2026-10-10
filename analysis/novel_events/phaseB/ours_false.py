#!/usr/bin/env python3
"""Our caller's (breakpoint caller + assembly) false calls in setE on clean
reads, sorted by cause; and a read-vs-assembly identity check for every
compared isolate. Report only, existing outputs (run from runroot after
rerun_compare.py; needs pysam, so $MTB_PY).

Isolates: the 46 of rerun_compare.py. Calls: typed, callable (the artifact
mask NOT applied), false by score.py's at_event.

Causes, tested in this order:
  palindrome      inside artifact_mask.py's site (+-50 bp)
  short_anchor    half or more of the reads clipped at the call (within 6 bp,
                  MAPQ >= 20, clip >= 10 bp) have fewer than 30 aligned
                  bases: reads of sequence absent from the reference pinned
                  by a chance short match, clipped on both sides, so they
                  form both clip clusters of an "insertion" themselves
  small_indel     the assembly-to-reference alignment (truth_paf, cs tag)
                  has an indel of 10-49 bp within 60 bp: real, but below the
                  truth set's 50 bp
  truth_mismatch  the isolate fails the identity check below
  other

Identity check: SNPs of the isolate's assembly against R (truth_paf cs
tag, MAPQ 60 alignments) vs SNPs from its reads (P2 haploid VCF, PASS,
alt genotype), in callable sequence. Discordant = in one set only. Isolates
above IDENT_MAX discordant are flagged: the reads and the assembly that
provides their truth are not the same strain.

Writes out/qc51/ours_false_setE.tsv and out/qc51/identity.tsv.
"""
import collections
import csv
import gzip
import re
import sys

import pysam

N = "../analysis/novel_events"
sys.path.insert(0, N)
import score  # noqa: E402

O = f"{N}/phaseB/out"
Q = f"{O}/qc51"
RM = "refbias/marinQC51/run/p1/refmap.tsv"
RM_BEFORE = "refbias/marinB63/run/p1/refmap.tsv"
P2 = "refbias/marinQC51/run/p2"
BAM = "refbias/work/marinQC51_p2"
IDENT_MAX = 20


def cs_ops(s):
    for line in open(f"{O}/truth_paf/{s}.paf"):
        f = line.rstrip("\n").split("\t")
        cs = next(x for x in f[12:] if x.startswith("cs:Z:"))[5:]
        yield int(f[7]), int(f[8]), int(f[11]), re.findall(r"([:*+\-~])([0-9]+|[a-z]+)", cs)


def asm_snps(s):
    out = set()
    for ts, _, mq, ops in cs_ops(s):
        if mq < 60:
            continue
        t = ts
        for op, b in ops:
            if op == ":":
                t += int(b)
            elif op == "*":
                out.add((t + 1, b[1].upper()))
                t += 1
            elif op == "-":
                t += len(b)
    return out


def small_indels(s, p, w=60):
    out = []
    for ts, te, _, ops in cs_ops(s):
        if not ts - w <= p <= te + w:
            continue
        t = ts
        for op, b in ops:
            if op == ":":
                t += int(b)
            elif op == "*":
                t += 1
            elif op in "-+":
                if abs(t - p) <= w and 10 <= len(b) < 50:
                    out.append(f"{'del' if op == '-' else 'ins'}{len(b)}")
                if op == "-":
                    t += len(b)
    return out


def read_snps(s):
    out = set()
    for line in gzip.open(f"{P2}/{s}.vcf.gz", "rt"):
        if line.startswith("#"):
            continue
        c = line.split("\t")
        if c[6] not in ("PASS", ".") or len(c[3]) != 1 or len(c[4]) != 1:
            continue
        if c[9].split(":")[0] in ("1", "1/1", "1|1"):
            out.add((int(c[1]), c[4]))
    return out


def short_anchor(bam, ctg, c):
    n = short = 0
    for q in {c["start"], c["end"]}:
        for a in bam.fetch(ctg, max(0, q - 6), q + 6):
            if a.is_secondary or a.is_supplementary or a.mapping_quality < 20:
                continue
            ct = a.cigartuples
            if (ct[0][0] == 4 and ct[0][1] >= 10 and abs(a.reference_start + 1 - q) <= 6) or \
                    (ct[-1][0] == 4 and ct[-1][1] >= 10 and abs(a.reference_end - q) <= 6):
                n += 1
                short += a.query_alignment_length < 30
    return short, n


def main():
    ref = {r["sample"]: r["reference"] for r in csv.DictReader(open(RM), delimiter="\t")}
    ref_b = {r["sample"]: r["reference"] for r in csv.DictReader(open(RM_BEFORE), delimiter="\t")}
    groups = {r["sample"]: r["group"] for r in csv.DictReader(open(f"{O}/samples.tsv"), delimiter="\t")}
    keep = [s for s in ref if ref[s] == ref_b.get(s)]
    masks = score.load_masks(RM, f"{N}/out/masks")
    art = {s: [tuple(map(int, x.split("\t")[1:3])) for x in open(f"{N}/out/masks/{ref[s]}.artifact.bed")]
           for s in keep}
    truth = collections.defaultdict(list)
    for g in ("novel", "setE_polished", "in_panel"):
        for r in csv.DictReader(open(f"{Q}/truth.{g}.tsv"), delimiter="\t"):
            truth[r["sample"]].append(dict(type=r["type"], start=int(r["r_start"]), end=int(r["r_end"]),
                                           size=max(int(r["ref_len"]), int(r["alt_len"]))))
    ident = []
    for s in keep:
        a = {x for x in asm_snps(s) if not score.masked(masks[s], x[0], x[0])}
        r = {x for x in read_snps(s) if not score.masked(masks[s], x[0], x[0])}
        ident.append(dict(sample=s, group=groups[s], reference=ref[s], assembly_snps=len(a), read_snps=len(r),
                          shared=len(a & r), assembly_only=len(a - r), reads_only=len(r - a),
                          discordant=len(a ^ r), flag="mismatch" if len(a ^ r) > IDENT_MAX else ""))
    ident.sort(key=lambda x: -x["discordant"])
    bad = {x["sample"] for x in ident if x["flag"]}
    rows = []
    for s in keep:
        if groups[s] != "setE_polished":
            continue
        bam = pysam.AlignmentFile(f"{BAM}/{s}.bam")
        ctg = bam.references[0]
        dys = score.read_vcf(f"{P2}/{s}.dysgu.vcf", True)
        dys20 = score.read_vcf(f"{P2}/{s}.dysgu.vcf", True, 0.2)
        for r in csv.DictReader(open(f"{Q}/asm/{s}.events.tsv"), delimiter="\t"):
            if r["type"] == "BND":
                continue
            c = dict(type=r["type"], start=int(r["start"]), end=int(r["end"]),
                     size=int(r["size"]) if r["size"].isdigit() and int(r["size"]) > 0 else None)
            if score.masked(masks[s], c["start"], c["end"]) or any(score.at_event(c, t) for t in truth[s]):
                continue
            sh, n = short_anchor(bam, ctg, c)
            si = small_indels(s, c["start"])
            if any(lo - 50 < c["end"] and c["start"] < hi + 50 for lo, hi in art[s]):
                cause = "palindrome"
            elif n and sh / n >= 0.5:
                cause = "short_anchor"
            elif si:
                cause = "small_indel"
            elif s in bad:
                cause = "truth_mismatch"
            else:
                cause = "other"
            rows.append(dict(sample=s, type=c["type"], start=c["start"], end=c["end"], size=r["size"],
                             support=r["support"], source=r["source"], cause=cause,
                             short_anchor_reads=f"{sh}/{n}", asm_small_indel=",".join(si),
                             dysgu=any(abs(x["start"] - c["start"]) <= 50 for x in dys),
                             dysgu_af20=any(abs(x["start"] - c["start"]) <= 50 for x in dys20)))
    for path, data in ((f"{Q}/ours_false_setE.tsv", rows), (f"{Q}/identity.tsv", ident)):
        with open(path, "w") as fo:
            w = csv.DictWriter(fo, fieldnames=list(data[0]), delimiter="\t")
            w.writeheader()
            w.writerows(data)
    print(collections.Counter(r["cause"] for r in rows).most_common())
    print("identity flagged:", ", ".join(f"{x['sample']} ({x['discordant']})" for x in ident if x["flag"]))


if __name__ == "__main__":
    main()
