#!/usr/bin/env python3
"""Build the accessory panel from the T2-validated candidates, classified by
how much of each locus is genuinely absent from H37Rv.

A two-class split (reference gap vs polymorphic) was the original design and the
data rejected it. Blasting each locus representative against H37Rv shows that
**98% of the non-distinguishing loci have >= 90% of their sequence already in
H37Rv** (median coverage 0.999). They are not sequence the reference lacks --
they are extra copies of sequence it already has. That is also why they failed to
separate carriers from non-carriers in T2: if the reference carries the sequence,
every genome blast-hits it.

So the classification that matters is novelty relative to H37Rv, and it implies
three different methods rather than one:

  NOVEL (H37Rv coverage < 0.05) -- genuinely absent from the reference. These
      need ALT contigs; nothing else can make their interior callable.

  MOSAIC (0.05 - 0.90) -- partly novel. Rearrangement junctions and partial
      insertions. ALT contigs help but the boundaries need care.

  COPY_NUMBER (> 0.90) -- the sequence is already in H37Rv; the variation is in
      how many copies exist. Adding a contig is the WRONG method: these want
      depth-based genotyping against the existing reference. Many sit in the
      PE_PGRS cluster at 3.93-3.95 Mb.

Loci are also flagged for whether the representative separates carriers from
non-carriers (T2 "confirmed"), which is orthogonal to novelty and decides whether
there is a genotype worth calling.
"""
import argparse, collections, csv, os, subprocess, sys


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].strip().split()[0]; buf = []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--validation", default="refbias/T2.accessory_validation.tsv")
    ap.add_argument("--candidates", default="refbias/t2/candidates.tsv")
    ap.add_argument("--fasta", default="refbias/t2/candidates.fasta")
    ap.add_argument("--h37rv", required=True)
    ap.add_argument("--blastn", required=True)
    ap.add_argument("--min-sep", type=float, default=0.50)
    ap.add_argument("--min-carrier-full", type=float, default=0.90)
    ap.add_argument("--novel-max", type=float, default=0.05,
                    help="H37Rv coverage below this is genuinely novel sequence")
    ap.add_argument("--cn-min", type=float, default=0.90,
                    help="H37Rv coverage above this is copy-number, not new sequence")
    ap.add_argument("--min-carriers", type=int, default=17,
                    help="carriers needed to call a locus common (~5%% of 332)")
    ap.add_argument("--outdir", default="refbias/panel")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    seqs = read_fasta(a.fasta)
    carriers = {r["allele_id"]: r["carriers"].split(",")
                for r in csv.DictReader(open(a.candidates), delimiter="\t")}

    by_locus = collections.defaultdict(list)
    for r in csv.DictReader(open(a.validation), delimiter="\t"):
        by_locus[int(r["pos"])].append(r)

    rows = []
    for pos, alleles in sorted(by_locus.items()):
        confirmed = [r for r in alleles
                     if float(r["carrier_frac_fulllength"]) >= a.min_carrier_full
                     and float(r["separation"]) >= a.min_sep]
        cls = "polymorphic" if confirmed else "reference_gap"
        # representative: the longest allele, preferring a confirmed one so the
        # emitted sequence is the haplotype that actually separates carriers.
        pool = confirmed or alleles
        rep = max(pool, key=lambda r: (int(r["len"]), int(r["AC"])))
        allc = set()
        for r in alleles:
            allc |= set(carriers[r["allele_id"]])
        rows.append(dict(locus_id=f"ACC_{pos:07d}", pos=pos, klass=cls,
                         rep_allele=rep["allele_id"], rep_len=int(rep["len"]),
                         n_alleles=len(alleles), n_confirmed=len(confirmed),
                         carriers_any=len(allc),
                         carrier_frac=round(len(allc) / 332, 4),
                         rep_separation=float(rep["separation"])))

    # Homology to H37Rv decides whether an entry can be a clean ALT contig at all.
    # An entry resembling reference sequence steals reads from it and raises the
    # core-genome false-positive rate -- the opposite of the intended effect.
    tmp = os.path.join(a.outdir, "_reps.fasta")
    with open(tmp, "w") as fh:
        for r in rows:
            fh.write(f">{r['locus_id']}\n{seqs[r['rep_allele']]}\n")
    # Coverage is the UNION of all HSPs, not the single best one. With
    # -max_hsps 1 a contig whose homology is spread over many short alignments
    # reads as novel because no individual HSP is long: measured that way, 2 of
    # 51 "novel" contigs actually reach 5% once HSPs are unioned. Union also
    # confirms how hard the mosaic class is -- median union coverage 0.981,
    # nearly the whole contig sitting somewhere in H37Rv.
    res = subprocess.run(
        [a.blastn, "-query", tmp, "-subject", a.h37rv, "-outfmt",
         "6 qseqid qstart qend qlen", "-evalue", "1e-10", "-dust", "no"],
        capture_output=True, text=True, check=True).stdout
    spans = collections.defaultdict(list)
    qlen = {}
    for line in res.strip().split("\n"):
        if not line:
            continue
        q, qs, qe, ql = line.split("\t")
        lo, hi = sorted((int(qs), int(qe)))
        spans[q].append((lo, hi))
        qlen[q] = int(ql)

    def union_len(iv):
        tot = cs = ce = 0
        started = False
        for lo, hi in sorted(iv):
            if not started or lo > ce:
                if started:
                    tot += ce - cs + 1
                cs, ce, started = lo, hi, True
            else:
                ce = max(ce, hi)
        return tot + (ce - cs + 1 if started else 0)

    hom = {q: union_len(v) / qlen[q] for q, v in spans.items()}
    for r in rows:
        f = min(hom.get(r["locus_id"], 0.0), 1.0)
        r["h37rv_cov"] = round(f, 4)
        r["novelty"] = ("novel" if f < a.novel_max
                        else "copy_number" if f > a.cn_min else "mosaic")
    os.remove(tmp)

    with open(os.path.join(a.outdir, "panel_manifest.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    for cls, fn in (("novel", "accessory_novel.fasta"),
                    ("mosaic", "accessory_mosaic.fasta")):
        with open(os.path.join(a.outdir, fn), "w") as fh:
            for r in rows:
                if r["novelty"] == cls:
                    fh.write(f">{r['locus_id']} pos={r['pos']} len={r['rep_len']} "
                             f"carriers={r['carriers_any']} class={r['klass']}\n"
                             f"{seqs[r['rep_allele']]}\n")

    print(f"  {len(rows)} loci from {sum(r['n_alleles'] for r in rows)} alleles")
    print(f"\n  {'novelty':<12}{'loci':>6}{'bp':>12}{'common':>8}{'common bp':>12}"
          f"{'method':>26}")
    METHOD = {"novel": "ALT contig",
              "mosaic": "ALT contig, check bounds",
              "copy_number": "depth genotyping"}
    for cls in ("novel", "mosaic", "copy_number"):
        sel = [r for r in rows if r["novelty"] == cls]
        com = [r for r in sel if r["carriers_any"] >= a.min_carriers]
        print(f"  {cls:<12}{len(sel):>6}{sum(r['rep_len'] for r in sel):>12,}"
              f"{len(com):>8}{sum(r['rep_len'] for r in com):>12,}"
              f"{METHOD[cls]:>26}")
    print(f"\n  polymorphic (separates carriers) x novelty:")
    for cls in ("novel", "mosaic", "copy_number"):
        sel = [r for r in rows if r["novelty"] == cls and r["klass"] == "polymorphic"]
        print(f"    {cls:<12}{len(sel):>6} loci")
    print(f"\n  written: {a.outdir}/panel_manifest.tsv, "
          f"accessory_novel.fasta, accessory_mosaic.fasta")


if __name__ == "__main__":
    sys.exit(main())
