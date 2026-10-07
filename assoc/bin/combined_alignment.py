#!/usr/bin/env python3
"""The cohort alignment plus every panel genome, for the cohort's tree.

WHY. A cohort tree built from the cohort's isolates alone places its deep
branches -- the ones the association tests and the ancestral reconstruction
lean on hardest -- from whatever sample of lineages the cohort happens to hold.
The panel's genomes span the complex, and their SNPs fix those branches. So the
tree is built from the cohort and the panel together, and only then cut down to
the tips the event writer reads (assoc/bin/prune_for_cohort.py). This used to
be analysis/combined_tree/build_alignment.py, run by hand; it carried its own
copy of the panel-genotype rule with faults the outgroup step had already had
fixed (audit TP-5: the last duplicate record won, missing calls were read as
REF), applied to all 332 panel genomes. Review 2 (R2-TREES-1) found the chain
still built cohort-only trees; this is the chain's step 2 now.

ONE RULE FOR EVERY PANEL GENOME. Each panel genome's base at a column comes from
assoc/bin/add_outgroup.py's outgroup_allele() -- the rule written and tested
for the outgroup (any exact record carrying the ALT is ALT; a missing call is
N; a record whose REF span covers the site is consulted; a base inside an indel
or a third allele is N). The outgroup is simply one of the panel rows.

COLUMNS.
  cohort  every column of the cohort alignment (bin/vcf_to_alignment.py's
          output, which carries the reference row --ref-sample).
  panel   biallelic SNPs in the panel VCF that no cohort column holds. The
          isolates are written REF there ONLY where no cohort record's REF
          span covers the position: the merged VCF holds every site any
          isolate varies at, so no record means every isolate matches the
          reference. Where a cohort record covers the position but is not in
          the cohort alignment (filtered there, or another allele), the
          isolates' states are not known, and the column is left out.
Off the reference contig (node-frame columns) the panel rows are N.

The panel's own copy of the reference (--h37rv-name) is left out: the cohort
alignment's reference row already stands for it, under the name the event
writer expects.

FILTERS over all taxa, as for any alignment the tree step reads: a column more
than --max-missing N, or constant among the called taxa, is dropped.

    $MTB_PY assoc/bin/combined_alignment.py \\
        --cohort-alignment data/trees/<c>.snps.fasta \\
        --cohort-sites data/trees/<c>.sites.tsv \\
        --vcf refbias/<c>/p5/merged.vcf.gz \\
        --panel-vcf <build>/assets/graph_collapsed.vcf.gz \\
        --outgroup <outgroup> \\
        --out data/trees/<c>.combined.fasta \\
        --sites-out data/trees/<c>.combined.sites.tsv
"""
import argparse, bisect, collections, csv, os, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from add_outgroup import outgroup_allele, read_fasta, PANEL_CONTIG

ACGT = set("ACGT")
MISSING = {".", "./.", ".|."}


def gt_index(gt):
    if gt in MISSING:
        return None
    return int(gt.replace("|", "/").split("/")[0])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort-alignment", required=True)
    ap.add_argument("--cohort-sites", required=True)
    ap.add_argument("--vcf", required=True, help="the cohort's merged VCF")
    ap.add_argument("--panel-vcf", required=True,
                    help="the build's collapsed graph VCF")
    ap.add_argument("--panel-contig", default=PANEL_CONTIG)
    ap.add_argument("--cohort-contig", default="NC_000962.3")
    ap.add_argument("--outgroup", required=True,
                    help="must be a panel genome; checked, and emitted like "
                         "every other panel row")
    ap.add_argument("--ref-sample", default="H37Rv",
                    help="the cohort alignment's reference row")
    ap.add_argument("--h37rv-name", default="GCF_000195955",
                    help="the panel's copy of the reference, left out")
    ap.add_argument("--max-missing", type=float, default=0.10)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sites-out", required=True)
    a = ap.parse_args()
    bcftools = os.environ.get("MTB_BCFTOOLS", "bcftools")

    # ---- the cohort alignment
    cseqs = read_fasta(a.cohort_alignment)
    csites = list(csv.DictReader(open(a.cohort_sites), delimiter="\t"))
    width = len(next(iter(cseqs.values())))
    if width != len(csites):
        sys.exit(f"FATAL: {a.cohort_alignment} has {width} columns but "
                 f"{a.cohort_sites} has {len(csites)} rows")
    if a.ref_sample not in cseqs:
        sys.exit(f"FATAL: the cohort alignment has no reference row "
                 f"{a.ref_sample} (vcf_to_alignment.py --ref-sample)")
    ckeys = {(s["chrom"], int(s["pos"]), s["ref"].upper(), s["alt"].upper())
             for s in csites}

    # ---- positions any cohort record covers (its REF span)
    covered = set()
    q = subprocess.run([bcftools, "query", "-r", a.cohort_contig,
                        "-f", "%POS\t%REF\n", a.vcf],
                       capture_output=True, text=True)
    if q.returncode != 0:
        sys.exit(f"FATAL: bcftools query failed on {a.vcf}: {q.stderr.strip()}")
    for line in q.stdout.splitlines():
        p, r = line.split("\t")
        p = int(p)
        covered.update(range(p, p + max(1, len(r))))

    # ---- the panel: every sample, one sweep of the reference contig
    panel = subprocess.run([bcftools, "query", "-l", a.panel_vcf],
                           capture_output=True, text=True, check=True).stdout.split()
    if a.outgroup not in panel:
        sys.exit(f"FATAL: outgroup {a.outgroup} is not a sample of {a.panel_vcf}")
    if a.outgroup == a.h37rv_name:
        sys.exit("FATAL: the outgroup cannot be the reference genome")
    emit = [i for i, s in enumerate(panel) if s != a.h37rv_name]
    clash = [panel[i] for i in emit if panel[i] in cseqs]
    if clash:
        sys.exit(f"FATAL: panel genomes already in the cohort alignment: {clash[:5]}")
    recs = []                    # (pos, ref, [alts], [gt index per sample])
    p_snps = []                  # panel-only candidate columns
    q = subprocess.run([bcftools, "query", "-r", a.panel_contig,
                        "-f", "%POS\t%REF\t%ALT[\t%GT]\n", a.panel_vcf],
                       capture_output=True, text=True)
    if q.returncode != 0:
        sys.exit(f"FATAL: bcftools query failed on {a.panel_vcf}: {q.stderr.strip()}")
    for line in q.stdout.splitlines():
        f = line.split("\t")
        pos, ref, alts = int(f[0]), f[1].upper(), f[2].upper().split(",")
        gts = [gt_index(g) for g in f[3:]]
        recs.append((pos, ref, alts, gts))
        if (len(ref) == 1 and len(alts) == 1 and len(alts[0]) == 1
                and ref in ACGT and alts[0] in ACGT
                and any(g == 1 for g in gts)):
            p_snps.append((pos, ref, alts[0]))

    n_skip = n_dup = 0
    extra = []
    for pos, ref, alt in p_snps:
        k = (a.cohort_contig, pos, ref, alt)
        if k in ckeys:
            n_dup += 1
            continue
        if pos in covered:
            n_skip += 1              # isolates' state there is not known
            continue
        extra.append(k)
    extra = sorted(set(extra))

    # ---- columns, and the panel records covering each
    columns = [(s["chrom"], int(s["pos"]), s["ref"].upper(), s["alt"].upper(),
                "cohort", int(s["column"])) for s in csites] \
        + [(c, p, r, x, "panel", None) for c, p, r, x in extra]
    site_pos = sorted({c[1] for c in columns if c[0] == a.cohort_contig})
    cover = collections.defaultdict(list)
    for rec in recs:
        pos, ref = rec[0], rec[1]
        i = bisect.bisect_left(site_pos, pos)
        while i < len(site_pos) and site_pos[i] <= pos + len(ref) - 1:
            cover[site_pos[i]].append(rec)
            i += 1

    names = list(cseqs) + [panel[i] for i in emit]
    out = {n: [] for n in names}
    kept = []
    drop = collections.Counter()
    ntax = len(names)
    for chrom, pos, ref, alt, src, ccol in columns:
        if src == "cohort":
            col = [cseqs[n][ccol] for n in cseqs]
        else:
            col = [ref] * len(cseqs)
        if chrom != a.cohort_contig:
            col += ["N"] * len(emit)
        else:
            rs = cover.get(pos, ())
            if (len(rs) == 1 and rs[0][0] == pos and rs[0][1] == ref
                    and rs[0][2] == [alt]):
                g = rs[0][3]         # the common case: one exact record
                col += [ref if g[i] == 0 else alt if g[i] == 1 else "N"
                        for i in emit]
            else:
                for i in emit:
                    call, _ = outgroup_allele(
                        pos, ref, alt, [(r[0], r[1], r[2], r[3][i]) for r in rs])
                    col.append({"REF": ref, "ALT": alt}.get(call, "N"))
        if col.count("N") / ntax > a.max_missing:
            drop[f"more than {a.max_missing:.0%} missing"] += 1
            continue
        if len({c for c in col if c != "N"}) < 2:
            drop["constant among called taxa"] += 1
            continue
        for n, c in zip(names, col):
            out[n].append(c)
        kept.append((chrom, pos, ref, alt, src))

    with open(a.out, "w") as fh:
        for n in names:
            fh.write(f">{n}\n{''.join(out[n])}\n")
    with open(a.sites_out, "w") as fh:
        fh.write("column\tchrom\tpos\tref\talt\tsource\n")
        for i, k in enumerate(kept):
            fh.write(f"{i}\t" + "\t".join(map(str, k)) + "\n")
    src = collections.Counter(k[4] for k in kept)
    print(f"  {len(cseqs)} cohort rows + {len(emit)} panel genomes "
          f"({a.h37rv_name} left out; outgroup {a.outgroup} among them)")
    print(f"  panel-only SNP columns: {len(extra):,} added, {n_dup:,} already "
          f"cohort columns, {n_skip:,} left out (a cohort record covers them)")
    print(f"  {ntax} taxa x {len(kept):,} columns ("
          + ", ".join(f"{k} {v:,}" for k, v in sorted(src.items())) + ")")
    for k, v in drop.most_common():
        print(f"    dropped, {k}: {v:,}")
    print(f"  -> {a.out}\n  -> {a.sites_out}")


if __name__ == "__main__":
    main()
