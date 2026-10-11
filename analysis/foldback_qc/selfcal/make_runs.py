#!/usr/bin/env python3
"""Fold-back cut-off calibration, step 2 (user approved all 54, 2026-10-11;
plan in ../CALIBRATION_PLAN.md): one Illumina run per panel genome.

From ../panel_illumina.tsv, every panel genome with paired-end Illumina WGS
runs of its own BioSample, except the Phase B isolates (already controls,
novel_events/phaseB/genomes.tsv same_isolate_in_cx333). One run per genome:
the paired-end ILLUMINA WGS run with the most bases (runs are not pooled, so
each sample is one library). ENA portal API for FASTQ links, md5s and the
library fields (kept for the read-quality table: library_selection,
library_construction_protocol, nominal_length, instrument). Metadata only.
Run from the repository root; writes selfcal/genomes.tsv (sample = ps_<GCF>).
"""
import csv
import io
import urllib.parse
import urllib.request

Q = "analysis/foldback_qc"
API = "https://www.ebi.ac.uk/ena/portal/api/search"
F = ("run_accession,sample_accession,study_accession,instrument_platform,instrument_model,"
     "library_layout,library_strategy,library_source,library_selection,library_name,"
     "library_construction_protocol,nominal_length,read_count,base_count,first_public,"
     "fastq_ftp,fastq_md5,fastq_bytes")


def ena(runs):
    out = []
    for i in range(0, len(runs), 40):
        q = " OR ".join(f'run_accession="{r}"' for r in runs[i:i + 40])
        data = urllib.parse.urlencode(dict(result="read_run", query=q, fields=F, format="tsv",
                                           limit=0)).encode()
        with urllib.request.urlopen(urllib.request.Request(API, data=data), timeout=120) as r:
            out += list(csv.DictReader(io.StringIO(r.read().decode()), delimiter="\t"))
    return {r["run_accession"]: r for r in out}


def main():
    own = {r["same_isolate_in_cx333"] for r in
           csv.DictReader(open("analysis/novel_events/phaseB/genomes.tsv"), delimiter="\t")}
    panel = [r for r in csv.DictReader(open(f"{Q}/panel_illumina.tsv"), delimiter="\t")
             if int(r["illumina_runs"]) and r["genome"] not in own]
    meta = ena(sorted({x for r in panel for x in r["runs"].split(",")}))
    rows = []
    for p in panel:
        cand = [meta[x] for x in p["runs"].split(",") if x in meta]
        cand = [m for m in cand if m["library_layout"] == "PAIRED" and m["library_strategy"] == "WGS"
                and m["instrument_platform"] == "ILLUMINA" and len(m["fastq_ftp"].split(";")) >= 2]
        if not cand:
            continue
        m = max(cand, key=lambda m: int(m["base_count"] or 0))
        # paired runs sometimes carry a third (unpaired) file first: keep _1/_2
        fq = [x for x in m["fastq_ftp"].split(";") if x.endswith(("_1.fastq.gz", "_2.fastq.gz"))]
        k = [m["fastq_ftp"].split(";").index(x) for x in fq]
        if len(fq) != 2:
            print("skip, no _1/_2 pair:", p["genome"], m["run_accession"], m["fastq_ftp"])
            continue
        md5 = [m["fastq_md5"].split(";")[i] for i in k]
        size = [m["fastq_bytes"].split(";")[i] for i in k]
        rows.append(dict(cohort_sample=f"ps_{p['genome']}", genome=p["genome"], biosample=p["biosample"],
                         strain=p["strain"], run=m["run_accession"], runs_available=len(cand),
                         r1=fq[0].rsplit("/", 1)[1], r2=fq[1].rsplit("/", 1)[1],
                         fastq_ftp=";".join(fq), fastq_md5=";".join(md5), fastq_bytes=";".join(size),
                         depth_raw=round(int(m["base_count"]) / 4411532, 1),
                         **{c: m[c].replace("\t", " ").replace("\n", " ") for c in (
                             "study_accession", "instrument_model", "library_selection",
                             "library_source", "library_name", "library_construction_protocol",
                             "nominal_length", "read_count", "base_count", "first_public")}))
    with open(f"{Q}/selfcal/genomes.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    gb = sum(int(x) for r in rows for x in r["fastq_bytes"].split(";")) / 1e9
    print(f"{len(rows)} genomes, one run each; {gb:.1f} GB of FASTQ to download")


if __name__ == "__main__":
    main()
