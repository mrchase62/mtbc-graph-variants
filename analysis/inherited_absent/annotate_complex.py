#!/usr/bin/env python3
"""What is the extra sequence in a complex local assembly?

local_assembly.py calls a region complex when a contig anchors both flanks
on H37Rv but is neither continuous H37Rv nor a clean deletion, and
unresolved when no contig anchors both. This takes every complex contig,
and for unresolved events every contig aligned within 1,500 bp of the
region (typically one flank anchored, then sequence H37Rv does not
continue), and finds the stretches of it, 50 bp or more, that its
minimap2 asm5 alignments to H37Rv leave uncovered. Each stretch is blasted
(blastn -task blastn) against three subjects:
  - H37Rv: short rearranged H37Rv pieces asm5 does not align;
  - the sample's matched reference R;
  - the build's accessory catalogue (one representative per locus).
The best hit per subject is reported: identity, length, query coverage,
position. Report only.
"""
import argparse
import csv
import os
import subprocess
import tempfile

import pysam

MIN_SEG = 50
PAD = 1500


def uncovered(qlen, spans):
    spans = sorted(spans)
    out, p = [], 0
    for a, b in spans:
        if a - p >= MIN_SEG:
            out.append((p, a))
        p = max(p, b)
    if qlen - p >= MIN_SEG:
        out.append((p, qlen))
    return out


def best_hits(blastn, query, subject):
    r = subprocess.run([blastn, "-task", "blastn", "-query", query,
                        "-subject", subject, "-evalue", "1e-10", "-outfmt",
                        "6 qseqid sseqid pident length qstart qend sstart send bitscore qlen"],
                       capture_output=True, text=True, check=True)
    best = {}
    for line in r.stdout.splitlines():
        c = line.split("\t")
        q, bits = c[0], float(c[8])
        if q not in best or bits > best[q][0]:
            cov = int(c[3]) / int(c[9])
            best[q] = (bits, f"{c[1]}:{c[6]}-{c[7]} id={c[2]} len={c[3]} qcov={cov:.2f}")
    return {q: v[1] for q, v in best.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True, help="local_assembly results.tsv")
    ap.add_argument("--workdir", required=True, help="local_assembly --workdir")
    ap.add_argument("--refmap", required=True)
    ap.add_argument("--refs", required=True, help="directory of <accession>.fasta")
    ap.add_argument("--h37rv", required=True)
    ap.add_argument("--accessory", required=True, help="accessory_catalogue.fasta")
    ap.add_argument("--blastn", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    ref_of = {r["sample"]: r["reference"] for r in
              csv.DictReader(open(a.refmap), delimiter="\t")}
    rows = []
    for ev in csv.DictReader(open(a.results), delimiter="\t"):
        if ev["verdict"] not in ("complex", "unresolved"):
            continue
        d = os.path.join(a.workdir, ev["event"])
        if not os.path.exists(os.path.join(d, "contigs.paf")):
            continue
        s0, e0 = int(ev["start"]), int(ev["end"])
        spans, qlen, near = {}, {}, set()
        for line in open(os.path.join(d, "contigs.paf")):
            c = line.split("\t")
            name = c[0].split("_length")[0]
            qlen[name] = int(c[1])
            spans.setdefault(name, []).append((int(c[2]), int(c[3])))
            if int(c[7]) < e0 + PAD and int(c[8]) > s0 - PAD:
                near.add(name)
        # complex contigs, or for an unresolved event any contig placed near
        # the region (one flank anchored, the rest not H37Rv-continuous)
        cx = [x.split(":")[0] for x in str(ev["contig_verdicts"]).split(",")
              if x.endswith(":complex")]
        if ev["verdict"] == "unresolved":
            cx = sorted(near)
        if not cx:
            continue
        seqs = {c.name.split("_length")[0]: c.sequence
                for c in pysam.FastxFile(os.path.join(d, "contigs.fasta"))}
        segs = []
        for name in cx:
            for s, e in uncovered(qlen.get(name, len(seqs[name])), spans.get(name, [])):
                segs.append((f"{name}:{s}-{e}", seqs[name][s:e]))
        if not segs:
            continue
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fq:
            for n, s in segs:
                fq.write(f">{n}\n{s}\n")
        R = ref_of.get(ev["sample"], "")
        hits = {"h37rv": best_hits(a.blastn, fq.name, a.h37rv),
                "accessory": best_hits(a.blastn, fq.name, a.accessory)}
        rf = os.path.join(a.refs, f"{R}.fasta")
        hits["matched_ref"] = best_hits(a.blastn, fq.name, rf) if os.path.exists(rf) else {}
        os.remove(fq.name)
        for n, s in segs:
            rows.append(dict(event=ev["event"], sample=ev["sample"], start=ev["start"],
                             end=ev["end"], reference=R, segment=n, seg_len=len(s),
                             h37rv=hits["h37rv"].get(n, ""),
                             matched_ref=hits["matched_ref"].get(n, ""),
                             accessory=hits["accessory"].get(n, "")))
    cols = ["event", "sample", "start", "end", "reference", "segment", "seg_len",
            "h37rv", "matched_ref", "accessory"]
    with open(a.out, "w") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} uncovered segments in "
          f"{len({r['event'] for r in rows})} complex events -> {a.out}")


if __name__ == "__main__":
    main()
