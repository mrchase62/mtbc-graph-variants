#!/usr/bin/env python3
"""Part D of ../FOLDBACK_SV_EFFECT_PLAN.md (fold-back and recall): runner
tables for cohort fbrecall. Run from runroot with $MTB_PY.

Isolates: selfcal54 QC PASS (reads refbias/selfcal54/reads) and the 9
chimeric controls (all reads, refbias/fbremove/reads). Each isolate's own
genome is a panel genome. The matched reference R is the nearest OTHER panel
genome by SNP-profile distance over the panel SNP matrix (P0 asset
panel_snps.vcf.gz, sites where both genotypes are called), excluding the own
genome and genomes within 50 SNPs of it, as stages 1-2 did when the sample was
itself a panel genome. Every assembly difference between the own genome and R
is then a real event the callers should find (truth.py).

Writes refbias/fbrecall/{cohort.tsv,crams.tsv,run/p1/refmap.tsv,run/.mtb_build}
and ../analysis/foldback_qc/recall/genomes.tsv (truth.py input).
"""
import csv
import datetime
import gzip
import os

import numpy as np

BUILD = "7713a8d71d8e-fix1"
B = f"refbias/build/{BUILD}"
C = "refbias/fbrecall"
O = "../analysis/foldback_qc/recall"
EXCL = 50
COLS = ("sample sample_lineage sample_sublineage meandepth reference ref_lineage ref_sublineage "
        "snp_distance inside_label rank2_reference rank2_distance margin tie_break ref_intervals "
        "rank2_intervals interval_delta").split()


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def panel_matrix():
    names, rows = None, []
    with gzip.open(f"{B}/assets/panel_snps.vcf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            c = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                names = c[9:]
                continue
            rows.append([-1 if g.startswith(".") else int(g.split("/")[0].split("|")[0]) for g in c[9:]])
    return names, np.array(rows, dtype=np.int8).T          # genomes x sites


def main():
    names, M = panel_matrix()
    idx = {n: i for i, n in enumerate(names)}

    def nearest(g):
        a = M[idx[g]]
        both = (a >= 0) & (M >= 0)
        d = ((a != M) & both).sum(axis=1)
        order = [j for j in np.argsort(d, kind="stable") if names[j] != g and d[j] > EXCL]
        j, k = order[0], order[1]
        return names[j], int(d[j]), names[k], int(d[k])

    pts = []
    sc_g = {r["cohort_sample"]: r["genome"] for r in tsv("../analysis/foldback_qc/selfcal/genomes.tsv")}
    for r in tsv("refbias/selfcal54/cohort.tsv"):
        pts.append((r["sample"], sc_g[r["sample"]], r, f"refbias/selfcal54/reads/{r['sample']}.bam"))
    own = {r["sample"]: r["reference"] for r in tsv("refbias/fbremove/run/p1/refmap.tsv")}
    fbc = {r["sample"]: r for r in tsv("refbias/fbremove/cohort.tsv")}
    for s in (x.split()[0] for x in open("../analysis/foldback_qc/removal/isolates.txt") if x.strip()):
        pts.append((s, own[f"{s}_all"], fbc[f"{s}_all"], f"refbias/fbremove/reads/{s}.bam"))

    os.makedirs(f"{C}/run/p1", exist_ok=True)
    root = os.environ["MTB_CRAM_ROOT"]
    with open(f"{C}/run/.mtb_build", "w") as fo:
        fo.write(f"build_id\t{BUILD}\nrecorded\t{datetime.datetime.now().astimezone().isoformat(timespec='seconds')}\n")
    with open(f"{C}/cohort.tsv", "w") as co, open(f"{C}/crams.tsv", "w") as cr, \
            open(f"{C}/run/p1/refmap.tsv", "w") as rm, open(f"{O}/genomes.tsv", "w") as ge:
        co.write("sample\tlineage\tdeepest\tmeandepth\tcoverage\tmapping_rate\terror_rate\tduplicate_rate\n")
        cr.write("sample\tcram_relpath\n")
        rm.write("\t".join(COLS) + "\n")
        ge.write("sample\town_genome\tfile\treference\tsnp_distance\n")
        for s, g, r, bam in pts:
            R, d, R2, d2 = nearest(g)
            assert os.path.exists(f"{B}/refs/{R}.fasta.bwt"), R
            assert os.path.exists(bam), bam
            co.write(f"{s}\t{r['lineage']}\t{r['deepest']}\t{r['meandepth']}\t100\t1\t0.002\t0\n")
            cr.write(f"{s}\t{os.path.relpath(os.path.realpath(bam), root)}\n")
            v = dict.fromkeys(COLS, "")
            v.update(sample=s, sample_lineage=r["lineage"], sample_sublineage=r["deepest"],
                     meandepth=r["meandepth"], reference=R, snp_distance=str(d),
                     inside_label="nearest_other_panel_genome", rank2_reference=R2,
                     rank2_distance=str(d2), margin=str(d2 - d))
            rm.write("\t".join(v[c] for c in COLS) + "\n")
            ge.write(f"{s}\t{g}\t{g}.fasta.gz\t{R}\t{d}\n")
            print(f"{s}\town {g}\tR {R}\td {d}")
    print(f"{len(pts)} isolates")


if __name__ == "__main__":
    main()
