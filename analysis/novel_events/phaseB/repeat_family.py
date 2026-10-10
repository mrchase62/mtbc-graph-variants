#!/usr/bin/env python3
"""Short-read deletion signal at every member of the tatC-helY element's
repeat family, in the 46 Phase B isolates compared on clean reads (report
only; run from runroot with $MTB_PY: needs pysam and ViennaRNA's RNA module).

Members: the 41 H37Rv copies in MtbPangenome/loci/hely_tatc/
repeat_family_loci.tsv (68 bp windows from `pos`), plus the two copies H37Rv
lacks: the tatC-helY element itself and the fadE22 copy (found by sequence).

1. Each member is located in each isolate's matched reference R by its two
   flanking 25-mers from H37Rv (exact, either strand, right flank within
   300 bp of the left). The two extra copies by their own sequence.
2. Member property: ViennaRNA MFE per nt of the window (H37Rv, or the
   sequence for the extra copies). Direct repeats are not scored here: the
   68 bp windows are BLAST-hit windows, not element boundaries (see
   loci/hely_tatc/REPEAT_FAMILY.md for the boundary-aware screen).
3. Per isolate and member:
   - depth_ratio: mean depth inside the member (10 bp in from each end) over
     the mean of 250 bp of flank each side (MAPQ >= 20, primary);
   - gap_reads: reads (MAPQ >= 20) with a CIGAR deletion of 40-120 bp
     overlapping the member (direct evidence of molecules lacking it);
   - clip_reads: reads soft-clipped >= 10 bp within 3 bp of either member
     end (weaker: also what a coverage dropout at a hairpin produces);
   - junction_reads: reads soft-clipped >= 15 bp within the member +-15 bp
     whose 15 clipped bases next to the alignment match R exactly 30-150 bp
     further on (right clip: downstream; left clip: upstream), i.e. reads of
     a molecule with that stretch deleted; junction_len is the commonest
     deletion length among them;
   - dysgu: the largest AF among dysgu DEL records (any FILTER) of 40-120 bp
     overlapping the member +-20 bp, and whether one is PASS;
   - ours: a typed call by our caller + assembly starting within 50 bp of it;
   - truth: a truth-set event within 50 bp (a real difference from R);
   - callable: outside callable_mask.py's mask.

Writes out/qc51/repeat_family.isolates.tsv and repeat_family.members.tsv.
"""
import collections
import csv
import re
import statistics
import sys

import pysam
import RNA

N = "../analysis/novel_events"
sys.path.insert(0, N)
import score  # noqa: E402

LOCI = "/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome/loci/hely_tatc/repeat_family_loci.tsv"
REFS = "refbias/build/7713a8d71d8e-fix1/refs"
Q = f"{N}/phaseB/out/qc51"
RM = "refbias/marinQC51/run/p1/refmap.tsv"
RM_BEFORE = "refbias/marinB63/run/p1/refmap.tsv"
P2 = "refbias/marinQC51/run/p2"
BAM = "refbias/work/marinQC51_p2"
FL = 25
EXTRA = {  # copies H37Rv lacks, located by sequence
    "tatC-helY": "CTCGCCTGGGCTGGCGAGCAGACGCAAAATCCCCCGCACGCCCGGCGTGTCGGGGGATTTTGCGTCTG",
    "fadE22": None,  # filled from a reference (see find_fade22)
}
COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    return "".join(x.strip() for x in open(path) if not x.startswith(">")).upper()


def locate(seq, left, right):
    """Member interval [a, b) in seq between the two flanks, + strand."""
    hits = []
    for strand, lq, rq in (("+", left, right), ("-", rc(right), rc(left))):
        for m in re.finditer(lq, seq):
            a = m.end()
            j = seq.find(rq, a, a + 300)
            if j >= 0:
                hits.append((a, j))
    return hits[0] if len(hits) == 1 else None


def locate_seq(seq, s):
    hits = [(m.start(), m.start() + len(s)) for q in (s, rc(s)) for m in re.finditer(q, seq)]
    return hits[0] if len(hits) == 1 else None


def depth(b, ctg, s, e):
    if e <= s:
        return 0.0
    cov = b.count_coverage(ctg, s, e, quality_threshold=0,
                           read_callback=lambda r: not (r.is_secondary or r.is_supplementary
                                                        or r.mapping_quality < 20))
    return sum(sum(x[i] for x in cov) for i in range(e - s)) / (e - s)


def del_reads(b, ctg, s, e):
    gap = clip = 0
    for r in b.fetch(ctg, max(0, s - 150), e + 150):
        if r.is_secondary or r.is_supplementary or r.mapping_quality < 20:
            continue
        p, hit, cl = r.reference_start, False, False
        for op, ln in r.cigartuples:
            if op == 2 and 40 <= ln <= 120 and p < e and p + ln > s:
                hit = True
            if op in (0, 2, 3, 7, 8):
                p += ln
        ct = r.cigartuples
        if (ct[-1][0] == 4 and ct[-1][1] >= 10 and min(abs(r.reference_end - s), abs(r.reference_end - e)) <= 3) or \
                (ct[0][0] == 4 and ct[0][1] >= 10 and min(abs(r.reference_start - s), abs(r.reference_start - e)) <= 3):
            cl = True
        gap += hit
        clip += cl and not hit
    return gap, clip


def junctions(b, ctg, sq, s, e):
    lens = collections.Counter()
    for r in b.fetch(ctg, max(0, s - 160), e + 160):
        if r.is_secondary or r.is_supplementary or r.mapping_quality < 20:
            continue
        ct, q = r.cigartuples, r.query_sequence
        if ct[-1][0] == 4 and ct[-1][1] >= 15 and s - 15 <= r.reference_end <= e + 15:
            k = q[len(q) - ct[-1][1]:][:15]
            j = sq.find(k, r.reference_end + 30, r.reference_end + 165)
            if j >= 0:
                lens[j - r.reference_end] += 1
                continue
        if ct[0][0] == 4 and ct[0][1] >= 15 and s - 15 <= r.reference_start <= e + 15:
            k = q[:ct[0][1]][-15:]
            j = sq.find(k, max(0, r.reference_start - 165), r.reference_start - 30)
            if j >= 0:
                lens[r.reference_start - (j + 15)] += 1
    return sum(lens.values()), (lens.most_common(1)[0][0] if lens else "")


def dysgu_dels(path):
    out = []
    for line in open(path):
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        info = dict(kv.split("=", 1) if "=" in kv else (kv, "") for kv in c[7].split(";"))
        if info.get("SVTYPE") != "DEL":
            continue
        ln = abs(int(info.get("SVLEN", 0) or 0))
        if not 40 <= ln <= 120:
            continue
        f = dict(zip(c[8].split(":"), c[9].split(":")))
        out.append((int(c[1]), int(info.get("END", c[1])), float(f.get("AF", 0)), c[6] == "PASS"))
    return out


def main():
    h37 = read_fasta(__import__("os").environ["MTB_H37RV"])
    members = []
    for r in csv.DictReader(open(LOCI), delimiter="\t"):
        a = int(r["pos"]) - 1
        b = a + 68
        members.append(dict(name=f"{r['pos']}_{r['gene']}", h37_pos=r["pos"], gene=r["gene"],
                            identity=r["identity_to_query"], left=h37[a - FL:a], right=h37[b:b + FL],
                            mfe=RNA.fold(h37[a:b])[1] / 68))
    ref = {r["sample"]: r["reference"] for r in csv.DictReader(open(RM), delimiter="\t")}
    ref_b = {r["sample"]: r["reference"] for r in csv.DictReader(open(RM_BEFORE), delimiter="\t")}
    groups = {r["sample"]: r["group"] for r in csv.DictReader(open(f"{N}/phaseB/out/samples.tsv"), delimiter="\t")}
    keep = [s for s in ref if ref[s] == ref_b.get(s)]
    masks = score.load_masks(RM, f"{N}/out/masks")
    truth = collections.defaultdict(list)
    for g in ("novel", "setE_polished", "in_panel"):
        for r in csv.DictReader(open(f"{Q}/truth.{g}.tsv"), delimiter="\t"):
            truth[r["sample"]].append((int(r["r_start"]), int(r["r_end"])))
    refseq = {g: read_fasta(f"{REFS}/{g}.fasta") for g in set(ref[s] for s in keep)}
    # the fadE22 copy: 68 bp in fadE22 (REPEAT_FAMILY.md section 3), absent from H37Rv;
    # anchored on its distinctive 20-mer and taken as the 68 bp around it
    core = "GATTTTGAGGCGATTCTGCG"
    for g, sq in refseq.items():
        hits = [m.start() for q in (core, rc(core)) for m in re.finditer(q, sq)]
        if len(hits) == 1:
            p = hits[0]
            w = sq[p - 30:p + 38] if core in sq else rc(sq[p - 18:p + 50])
            EXTRA["fadE22"] = w
            break
    for name, s in EXTRA.items():
        if s:
            members.append(dict(name=name, h37_pos="absent", gene=name, identity="", seq=s,
                                mfe=RNA.fold(s)[1] / len(s)))
    rows = []
    for s in keep:
        g = ref[s]
        sq = refseq[g]
        bam = pysam.AlignmentFile(f"{BAM}/{s}.bam")
        ctg = bam.references[0]
        dys = dysgu_dels(f"{P2}/{s}.dysgu.vcf")
        ours = [c for c in score.read_proto(f"{Q}/asm/{s}.events.tsv") if c["type"] != "BND"]
        for m in members:
            iv = locate_seq(sq, m["seq"]) if "seq" in m else locate(sq, m["left"], m["right"])
            if iv is None:
                rows.append(dict(sample=s, group=groups[s], member=m["name"], found=0))
                continue
            a, b = iv
            fl = (depth(bam, ctg, a - 260, a - 10) + depth(bam, ctg, b + 10, b + 260)) / 2
            inside = depth(bam, ctg, a + 10, b - 10) if b - a > 30 else None
            d = [x for x in dys if x[0] <= b + 20 and x[1] >= a - 20]
            rows.append(dict(
                sample=s, group=groups[s], member=m["name"], found=1, r_start=a + 1, r_end=b,
                length=b - a, callable=int(not score.masked(masks[s], a + 1, b)),
                depth_ratio=round(inside / fl, 3) if inside is not None and fl else "",
                flank_depth=round(fl, 1), **dict(zip(("gap_reads", "clip_reads"), del_reads(bam, ctg, a, b))),
                **dict(zip(("junction_reads", "junction_len"), junctions(bam, ctg, sq, a, b))),
                dysgu_af=max((x[2] for x in d), default=0.0), dysgu_pass=int(any(x[3] for x in d)),
                ours=int(any(a - 50 <= c["start"] <= b + 50 for c in ours)),
                truth=int(any(x <= b + 50 and y >= a - 50 for x, y in truth[s]))))
    with open(f"{Q}/repeat_family.isolates.tsv", "w") as fo:
        cols = ["sample", "group", "member", "found", "r_start", "r_end", "length", "callable", "depth_ratio",
                "flank_depth", "gap_reads", "clip_reads", "junction_reads", "junction_len", "dysgu_af", "dysgu_pass", "ours", "truth"]
        w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
    summ = []
    for m in members:
        x = [r for r in rows if r["member"] == m["name"] and r["found"] and r["length"] > 30]
        dr = [r["depth_ratio"] for r in x if r["depth_ratio"] != ""]
        summ.append(dict(
            member=m["name"], gene=m["gene"], h37rv_pos=m["h37_pos"], identity_to_element=m["identity"],
            mfe_per_nt=round(m["mfe"], 2), isolates_carrying=len(x),
            callable=sum(r["callable"] for r in x),
            depth_ratio_median=round(statistics.median(dr), 2) if dr else "",
            depth_ratio_min=min(dr) if dr else "",
            isolates_gap_reads=sum(r["gap_reads"] > 0 for r in x),
            gap_reads_total=sum(r["gap_reads"] for r in x),
            isolates_junction_reads=sum(r["junction_reads"] > 0 for r in x),
            junction_reads_total=sum(r["junction_reads"] for r in x),
            junction_per_100x_median=round(statistics.median(100 * r["junction_reads"] / r["flank_depth"]
                                                             for r in x if r["flank_depth"]), 2) if x else "",
            junction_len=collections.Counter(r["junction_len"] for r in x if r["junction_len"] != "").most_common(1)[0][0]
            if any(r["junction_len"] != "" for r in x) else "",
            isolates_clip_reads=sum(r["clip_reads"] > 0 for r in x),
            clip_reads_total=sum(r["clip_reads"] for r in x),
            isolates_dysgu_del=sum(r["dysgu_af"] > 0 for r in x),
            isolates_dysgu_pass=sum(r["dysgu_pass"] for r in x),
            dysgu_af_max=max((r["dysgu_af"] for r in x), default=0),
            isolates_ours=sum(r["ours"] for r in x), isolates_truth=sum(r["truth"] for r in x)))
    summ.sort(key=lambda r: (-r["junction_reads_total"], -r["clip_reads_total"]))
    with open(f"{Q}/repeat_family.members.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(summ[0]), delimiter="\t")
        w.writeheader()
        w.writerows(summ)


if __name__ == "__main__":
    main()
