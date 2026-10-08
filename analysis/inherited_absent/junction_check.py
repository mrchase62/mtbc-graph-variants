#!/usr/bin/env python3
"""For one sample, decide whether each ABSENT region is really missing at its
locus, from soft clips and split reads at the matched reference's deletion
boundaries in the sample's H37Rv alignment.

count_false_absent.py found ABSENT cells where the sample has normal depth.
Depth alone cannot tell "the sequence is here" from "the sequence is
elsewhere in this genome and its reads map here". The junctions can:

  present_here  reads run straight across both boundaries of R's deletion,
                and fewer than 3 reads join the boundaries: the sample has the
                sequence at this locus, so ABSENT is wrong
  deleted_here  >= 3 reads join the two boundaries (a split read whose SA lands
                at the other boundary, or a CIGAR deletion from one to the
                other) and few reads bridge them: the locus is deleted, and the
                depth comes from a copy elsewhere
  mixed         both: junction reads and bridging reads (a duplication, a
                mixture, or a repeat)
  unresolved    neither, e.g. a boundary in a repeat with no MAPQ>=20 reads

R's deletion is read off the build's collapsed graph VCF, which genotypes
every panel genome: the largest R-carried deletion overlapping the cell.
Cells with no R deletion overlapping are reported as R_no_genotype when R
has no genotype at an overlapping graph site (its path does not cross the
site in place: a rearrangement, an inversion or an assembly break, which no
single junction can test), else no_R_deletion.

Bridging at a boundary B: a MAPQ>=20 primary read aligned without a gap from
B-10 to B+10. Clipped at a boundary: a soft clip of >= 5 bp ending or starting
within 10 bp of it. Report only; nothing in the pipeline reads this.
"""
import argparse
import collections
import csv
import gzip
import sys

import pysam

TOL = 10
MIN_CLIP = 5
MIN_MAPQ = 20


def r_deletions(vcf, acc):
    """Deleted H37Rv spans (first, last base, 1-based) carried by acc, and
    the spans of sites where acc has no genotype (its path does not cross
    the site in place)."""
    out, missing = [], []
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                cols = line.rstrip("\n").split("\t")[9:]
                if acc not in cols:
                    # R is the VCF's own reference path (H37Rv): it carries
                    # no deletion relative to itself
                    return out, missing
                col = 9 + cols.index(acc)
                continue
            c = line.rstrip("\n").split("\t")
            gt = c[col].split(":")[0]
            if gt == ".":
                missing.append((int(c[1]), int(c[1]) + len(c[3]) - 1))
                continue
            if gt == "0":
                continue
            alt = c[4].split(",")[int(gt) - 1]
            ref = c[3]
            if len(ref) <= len(alt):
                continue
            pos = int(c[1])
            # trim the shared prefix; what remains of REF is deleted
            k = 0
            while k < len(alt) and ref[k] == alt[k]:
                k += 1
            s, e = pos + k, pos + len(ref) - 1 - (len(alt) - k)
            if e >= s:
                out.append((s, e))
    out.sort()
    missing.sort()
    return out, missing


def aligned_across(read, a, b):
    """True if the read is aligned without a gap over [a, b] (1-based)."""
    for bs, be in read.get_blocks():  # 0-based half-open
        if bs + 1 <= a and be >= b:
            return True
    return False


def boundary_evidence(bam, ctg, s, e):
    """Reads bridging s, bridging e, and reads joining s-1 to e+1."""
    left, right = s - 1, e + 1          # last kept base, first kept base
    bridge_s = bridge_e = 0
    clip_l = clip_r = 0
    joined = set()
    lo, hi = max(0, s - 300), e + 300
    windows = [(lo, s + 300)] if e - s < 600 else [(lo, s + 300), (e - 300, hi)]
    seen = set()
    for wlo, whi in windows:
        for r in bam.fetch(ctg, max(0, wlo), whi):
            if (r.is_unmapped or r.is_secondary or r.is_supplementary
                    or r.mapping_quality < MIN_MAPQ):
                continue
            key = (r.query_name, r.is_read1)
            if key in seen:
                continue
            seen.add(key)
            if aligned_across(r, s - TOL, s + TOL):
                bridge_s += 1
            if aligned_across(r, e - TOL, e + TOL):
                bridge_e += 1
            ct = r.cigartuples
            # CIGAR deletion from the left boundary to the right one
            p = r.reference_start
            for op, n in ct:
                if op == 2 and n >= 1 and abs(p - left) <= TOL \
                        and abs(p + n + 1 - right) <= TOL:
                    joined.add(r.query_name)
                if op in (0, 2, 3, 7, 8):
                    p += n
            sa = r.get_tag("SA").rstrip(";").split(";") if r.has_tag("SA") else []
            sa_pos = [int(x.split(",")[1]) for x in sa if x.split(",")[0] == ctg]
            if ct[-1][0] == 4 and ct[-1][1] >= MIN_CLIP \
                    and abs(r.reference_end - left) <= TOL:
                clip_r += 1
                if any(abs(q - right) <= 2 * TOL + ct[-1][1] for q in sa_pos):
                    joined.add(r.query_name)
            if ct[0][0] == 4 and ct[0][1] >= MIN_CLIP \
                    and abs(r.reference_start + 1 - right) <= TOL:
                clip_l += 1
                if any(abs(q - left) <= 2 * TOL + ct[0][1] for q in sa_pos):
                    joined.add(r.query_name)
    return bridge_s, bridge_e, clip_r, clip_l, len(joined)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--records", required=True, help="count_false_absent records.tsv")
    ap.add_argument("--summary", required=True, help="count_false_absent summary.tsv")
    ap.add_argument("--refmap", required=True)
    ap.add_argument("--graph-vcf", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--verdicts", default="contradicted,partial,supported")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    want = set(a.verdicts.split(","))
    recs = [r for r in csv.DictReader(open(a.records), delimiter="\t")
            if r["sample"] == a.sample and r["verdict"] in want]
    base = next(float(r["baseline_depth"]) for r in
                csv.DictReader(open(a.summary), delimiter="\t")
                if r["sample"] == a.sample)
    acc = next(r["reference"] for r in
               csv.DictReader(open(a.refmap), delimiter="\t")
               if r["sample"] == a.sample)
    dels, missing = r_deletions(a.graph_vcf, acc)
    bam = pysam.AlignmentFile(a.bam)
    ctg = bam.references[0]

    # each cell -> the largest R deletion overlapping it
    events, cell_event = {}, []
    for r in recs:
        s, e = int(r["start"]), int(r["end"])
        hits = [d for d in dels if d[0] <= e and d[1] >= s]
        ev = max(hits, key=lambda d: d[1] - d[0]) if hits else None
        cell_event.append((r, ev))
        if ev is not None:
            events.setdefault(ev, None)
        elif any(m[0] <= e and m[1] >= s for m in missing):
            cell_event[-1] = (r, "R_no_genotype")
    min_bridge = max(5, 0.25 * base)
    for ev in events:
        bs, be, cr, cl, jn = boundary_evidence(bam, ctg, *ev)
        bridge = min(bs, be)
        if jn >= 3 and bridge >= min_bridge:
            v = "mixed"
        elif jn >= 3:
            v = "deleted_here"
        elif bridge >= min_bridge:
            v = "present_here"
        else:
            v = "unresolved"
        events[ev] = (bs, be, cr, cl, jn, v)

    with open(a.out, "w") as fo:
        fo.write("sample\treference\tid\tregion\tdepth_verdict\tmean_depth"
                 "\tr_del_start\tr_del_end\tr_del_len\tbridge_start\tbridge_end"
                 "\tclip_at_start\tclip_at_end\tjoined_reads\tjunction_verdict\n")
        for r, ev in cell_event:
            if ev is None or isinstance(ev, str):
                tail = ["", "", "", "", "", "", "", "",
                        ev if isinstance(ev, str) else "no_R_deletion"]
            else:
                bs, be, cr, cl, jn, v = events[ev]
                tail = [ev[0], ev[1], ev[1] - ev[0] + 1, bs, be, cr, cl, jn, v]
            fo.write("\t".join(map(str, [a.sample, acc, r["id"], r["region"],
                                         r["verdict"], r["mean_depth"]] + tail))
                     + "\n")
    c = collections.Counter((r["verdict"], (events[ev][5] if isinstance(ev, tuple)
                                            else ev or "no_R_deletion"))
                            for r, ev in cell_event)
    print(a.sample, acc, f"baseline {base:.0f}", f"{len(events)} R deletions",
          dict(c), file=sys.stderr)


if __name__ == "__main__":
    main()
