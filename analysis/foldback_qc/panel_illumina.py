#!/usr/bin/env python3
"""Step 1 of the fold-back cutoff calibration (user, 2026-10-11): which CX333
panel genomes (finished assemblies) have Illumina reads of the same
BioSample? Mapped to their own genome, such reads give every call as noise
(phaseB/self_noise.py), so isolates whose reads fall at 1-6 fold-backs per
1,000 can measure the false-call rate where no control exists yet.

Metadata only, no reads downloaded:
  - panel: <build>/assets/accessions.txt (333 GCF accessions);
  - BioSample: data/ncbi/assembly_summary_refseq.txt (column biosample);
  - Illumina runs: ENA portal API read_run search by sample_accession
    (as foldback_sim/ena_library.py), instrument_platform ILLUMINA,
    library_strategy WGS;
  - whether the BioSample is one of the 54,461 production CRAMs
    ($MTB_CRAM_ROOT/metadata/samples.json.gz; reads already aligned to
    H37Rv, no download needed), and its fold-back rate if it was in the 500
    measured (novel_events/controls/production_chimera.tsv, old count).
Run from runroot with $MTB_PY. Writes foldback_qc/panel_illumina.tsv (one
row per panel genome) and prints a summary.
"""
import csv
import gzip
import json
import os
import sys

sys.path.insert(0, "../analysis/foldback_sim")
from ena_library import fetch  # noqa: E402

B = "refbias/build/7713a8d71d8e-fix1"
Q = "../analysis/foldback_qc"


def main():
    panel = [x.strip() for x in open(f"{B}/assets/accessions.txt") if x.strip()]
    bios = {}
    for line in open("data/ncbi/assembly_summary_refseq.txt"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        bios[f[0].split(".")[0]] = (f[2], f[7], f[8], f[9])
    prod = json.load(gzip.open(os.path.join(os.environ["MTB_CRAM_ROOT"], "metadata/samples.json.gz"), "rt"))
    measured = {r["sample"]: r["inv_chimera_per_1k"] for r in
                csv.DictReader(open("../analysis/novel_events/controls/production_chimera.tsv"), delimiter="\t")}
    want = sorted({bios[g][0] for g in panel if g in bios and bios[g][0].startswith("SAM")})
    runs = fetch("sample_accession", want)
    by = {}
    for r in runs:
        if r["instrument_platform"] == "ILLUMINA" and r["library_strategy"] == "WGS":
            by.setdefault(r["sample_accession"], []).append(r)
    rows = []
    for g in panel:
        b, org, strain, iso = bios.get(g, ("", "", "", ""))
        rs = by.get(b, [])
        rows.append(dict(
            genome=g, biosample=b, strain=(strain or iso).replace("strain=", ""),
            illumina_runs=len(rs),
            runs=",".join(r["run_accession"] for r in rs),
            instruments=",".join(sorted({r["instrument_model"] for r in rs})),
            layouts=",".join(sorted({r["library_layout"] for r in rs})),
            gbases=round(sum(int(r["base_count"] or 0) for r in rs) / 1e9, 2),
            studies=",".join(sorted({r["study_accession"] for r in rs})),
            in_production=int(b in prod), production_cram=prod.get(b, {}).get("path", ""),
            foldback_old_per_1k=measured.get(b, "")))
    with open(f"{Q}/panel_illumina.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    n_b = sum(1 for r in rows if r["biosample"])
    n_ill = sum(1 for r in rows if r["illumina_runs"])
    n_prod = sum(1 for r in rows if r["in_production"])
    n_pe = sum(1 for r in rows if "PAIRED" in r["layouts"])
    print(f"panel genomes {len(rows)}; with a BioSample {n_b}; with Illumina WGS runs {n_ill} "
          f"(paired-end {n_pe}); BioSample among the production CRAMs {n_prod}; "
          f"of those, fold-back already measured {sum(1 for r in rows if r['foldback_old_per_1k'])}")
    ins = {}
    for r in rows:
        for m in r["instruments"].split(","):
            if m:
                ins[m] = ins.get(m, 0) + 1
    print("instruments (genomes):", dict(sorted(ins.items(), key=lambda x: -x[1])))


if __name__ == "__main__":
    main()
