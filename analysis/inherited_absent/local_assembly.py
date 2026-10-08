#!/usr/bin/env python3
"""Local reassembly of ABSENT regions for one sample: is the region at its
H37Rv locus or not?

For each event (make_events.py) in this sample:
1. reads from the sample's H37Rv BAM in [start - 1,500, end + 1,500], any
   MAPQ, with their mates wherever those map, and unmapped mates (bwa places
   an unmapped mate at its partner's position, so the window fetch includes
   them);
2. SPAdes --isolate, k 21,33,55,77;
3. contigs >= 300 bp aligned to H37Rv with minimap2 -x asm5 -c;
4. each contig classified against the region [start, end], with F = 300 bp
   of flank required on both sides:
     present_here  one alignment runs from start - F to end + F with at most
                   10% of the region's bases deleted (counted over +-10 bp
                   for regions under 50 bp)
     deleted_here  that alignment deletes >= 80% of the region, or the same
                   contig anchors both flanks (F bp each, same strand) and its
                   alignments cover < 20% of the region
     complex       anchors both flanks, but neither of the above
     (no verdict)  a contig that does not anchor both flanks
   The event's verdict joins its contigs' verdicts: present_here,
   deleted_here, mixed (both), complex, or unresolved (no contig anchors
   both flanks).

A contig must reach unique flank on both sides, so a copy of the region
elsewhere in the genome cannot pass as present here. Report only.
"""
import argparse
import csv
import os
import re
import shutil
import subprocess
import sys

import pysam

PAD = 1500
F = 300
MIN_CONTIG = 300
MAX_REGION = 15000


def collect_reads(bam, ctg, lo, hi):
    names, mate_pos = set(), set()
    for r in bam.fetch(ctg, max(0, lo - 1), hi):
        if r.is_secondary:
            continue
        names.add(r.query_name)
        if (r.is_paired and not r.mate_is_unmapped
                and (r.next_reference_name != ctg
                     or not lo - 1 <= r.next_reference_start <= hi)):
            mate_pos.add((r.next_reference_name, r.next_reference_start))
    reads = {}

    def keep(r):
        if r.is_secondary or r.is_supplementary or r.query_name not in names:
            return
        seq = r.get_forward_sequence()
        qual = r.get_forward_qualities()
        if not seq:
            return
        q = "".join(chr(min(x, 41) + 33) for x in qual) if qual is not None \
            else "I" * len(seq)
        mate = 1 if (r.is_read1 or not r.is_paired) else 2
        reads.setdefault(r.query_name, {})[mate] = (seq, q)

    for r in bam.fetch(ctg, max(0, lo - 1), hi):
        keep(r)
    for c, p in list(mate_pos)[:5000]:
        for r in bam.fetch(c, p, p + 1):
            keep(r)
    return reads


def write_fastq(reads, d):
    p1, p2, ps = (os.path.join(d, x) for x in ("r1.fq", "r2.fq", "s.fq"))
    n_pair = n_single = 0
    with open(p1, "w") as f1, open(p2, "w") as f2, open(ps, "w") as fs:
        for name, m in reads.items():
            if 1 in m and 2 in m:
                f1.write(f"@{name}/1\n{m[1][0]}\n+\n{m[1][1]}\n")
                f2.write(f"@{name}/2\n{m[2][0]}\n+\n{m[2][1]}\n")
                n_pair += 1
            else:
                for k, (s, q) in m.items():
                    fs.write(f"@{name}/{k}\n{s}\n+\n{q}\n")
                    n_single += 1
    return p1, p2, ps, n_pair, n_single


def ref_cover(ts, cigar, s, e):
    """Aligned ref bases of [s, e] (1-based) and deleted bases inside it."""
    p = ts  # 0-based
    cov = dele = 0
    big_ins = 0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=X":
            a, b = max(p + 1, s), min(p + n, e)
            cov += max(0, b - a + 1)
            p += n
        elif op in "DN":
            a, b = max(p + 1, s), min(p + n, e)
            dele += max(0, b - a + 1)
            p += n
        elif op == "I":
            if s <= p <= e and n >= 50:
                big_ins += n
    return cov, dele, big_ins


def classify(alns, s, e):
    L = e - s + 1
    # A short deletion can be placed a few bases off in a repeat or a
    # homopolymer, so for regions under 50 bp the deleted bases are counted
    # over the region +-10 bp. Present means at most 10% of the region
    # deleted: none at all below 10 bp. (The first version allowed up to 20
    # deleted bases, which called 1-22 bp deletions present.)
    pad = 10 if L < 50 else 0
    for a in alns:
        if a["ts"] + 1 <= s - F and a["te"] >= e + F:
            cov, dele, ins = ref_cover(a["ts"], a["cg"], s - pad, e + pad)
            if dele <= 0.1 * L and ins == 0:
                return "present_here"
            if dele >= 0.8 * L:
                return "deleted_here"
    for strand in "+-":
        sa = [a for a in alns if a["strand"] == strand]
        left = any(a["ts"] + 1 <= s - F and a["te"] >= s - 1 - 10 for a in sa)
        right = any(a["ts"] + 1 <= e + 1 + 10 and a["te"] >= e + F for a in sa)
        if left and right:
            cov = sum(ref_cover(a["ts"], a["cg"], s, e)[0] for a in sa)
            if cov < 0.2 * L:
                return "deleted_here"
            return "complex"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--h37rv-mmi", required=True)
    ap.add_argument("--spades", required=True, help="python and spades.py, space-separated")
    ap.add_argument("--minimap2", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--reclassify", action="store_true",
                    help="reuse each event's saved contigs.paf; no reassembly")
    a = ap.parse_args()
    evs = [r for r in csv.DictReader(open(a.events), delimiter="\t")
           if r["sample"] == a.sample]
    bam = pysam.AlignmentFile(a.bam)
    ctg = bam.references[0]
    rows = []
    for ev in evs:
        s, e = int(ev["start"]), int(ev["end"])
        d = os.path.join(a.workdir, ev["event"])
        os.makedirs(d, exist_ok=True)
        rec = dict(ev, n_pairs=0, n_single=0, n_contigs=0, anchoring=0,
                   contig_verdicts="", verdict="unresolved", note="")
        if e - s + 1 > MAX_REGION:
            rec["note"] = "region too large"
            rows.append(rec)
            continue
        paf = os.path.join(d, "contigs.paf")
        if a.reclassify and os.path.exists(paf):
            rec["n_contigs"] = sum(1 for _ in pysam.FastxFile(
                os.path.join(d, "contigs.fasta")))
            rec["note"] = "reclassified"
            classify_event(rec, paf, ctg, s, e)
            rows.append(rec)
            continue
        reads = collect_reads(bam, ctg, s - PAD, e + PAD)
        p1, p2, ps, npair, nsing = write_fastq(reads, d)
        rec.update(n_pairs=npair, n_single=nsing)
        sp = os.path.join(d, "spades")
        cmd = a.spades.split() + ["--isolate", "-k", "21,33,55,77", "-t", "1",
                                  "-m", "4", "-o", sp]
        if npair:
            cmd += ["-1", p1, "-2", p2]
        if nsing:
            cmd += ["-s", ps]
        r = subprocess.run(cmd, capture_output=True, text=True)
        contigs = os.path.join(sp, "contigs.fasta")
        if r.returncode != 0 or not os.path.exists(contigs):
            rec["note"] = f"spades exit {r.returncode}"
            rows.append(rec)
            continue
        keep = os.path.join(d, "contigs.fasta")
        n = 0
        with open(keep, "w") as fo:
            for c in pysam.FastxFile(contigs):
                if len(c.sequence) >= MIN_CONTIG:
                    fo.write(f">{c.name}\n{c.sequence}\n")
                    n += 1
        rec["n_contigs"] = n
        # SPAdes intermediates are large; keep only the filtered contigs
        if os.path.realpath(sp).startswith(os.path.realpath(a.workdir) + os.sep):
            shutil.rmtree(sp)
        for f in (p1, p2, ps):
            os.remove(f)
        with open(paf, "w") as fo:
            subprocess.run([a.minimap2, "-x", "asm5", "-c", "--secondary=no",
                            a.h37rv_mmi, keep], stdout=fo,
                           stderr=subprocess.DEVNULL, check=True)
        classify_event(rec, paf, ctg, s, e)
        rows.append(rec)
        print(a.sample, ev["event"], s, e, rec["verdict"], f"{npair}p/{nsing}s",
              f"{n} contigs", file=sys.stderr)
    cols = list(evs[0].keys()) + ["n_pairs", "n_single", "n_contigs",
                                  "anchoring", "contig_verdicts", "verdict", "note"]
    with open(a.out, "w") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def classify_event(rec, paf, ctg, s, e):
        alns = {}
        for line in open(paf):
            c = line.rstrip("\n").split("\t")
            if c[5] != ctg:
                continue
            cg = next((x[5:] for x in c[12:] if x.startswith("cg:Z:")), "")
            alns.setdefault(c[0], []).append(dict(
                strand=c[4], ts=int(c[7]), te=int(c[8]), cg=cg))
        verdicts = {}
        for name, al in alns.items():
            v = classify(al, s, e)
            if v:
                verdicts[name] = v
        rec["anchoring"] = len(verdicts)
        rec["contig_verdicts"] = ",".join(f"{k.split('_length')[0]}:{v}"
                                          for k, v in verdicts.items())
        vs = set(verdicts.values())
        if {"present_here", "deleted_here"} <= vs:
            rec["verdict"] = "mixed"
        elif "present_here" in vs:
            rec["verdict"] = "present_here"
        elif "deleted_here" in vs:
            rec["verdict"] = "deleted_here"
        elif "complex" in vs:
            rec["verdict"] = "complex"


if __name__ == "__main__":
    main()
