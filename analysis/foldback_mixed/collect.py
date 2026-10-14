#!/usr/bin/env python3
"""Step 1 of the fold-back vs mixed-call proposal (Claude Doc "Fold-back
chimeras and mixed calls"; user, 2026-10-13: "run 1-7"). Analysis only.

One pass over the collection's Mutect2 calls (tumour-only,
--min-base-quality 20, per the user) for every isolate that passes the
standard read QC and has a fold-back rate. Only calls whose whole REF span is
outside the build's repeat mask are counted.

Classes by AF: fixed >= 0.9; mixed_high 0.5-0.9; mixed_mid 0.25-0.5;
mixed_low 0.1-0.25; trace < 0.1 (counted and reported only: 100x does not
support it, user). Types: snv (1 bp REF and ALT), indel (anything else
under 50 bp, MNVs included). PASS and non-PASS apart.

Outputs (--outdir):
  isolates.tsv        one row per isolate: QC covariates, fold-back rates,
                      TB-Profiler mixed fields, counts per type x filter x class
  mixed_calls.tsv.gz  every unmasked call with 0.1 <= AF < 0.9: position,
                      alleles, PASS, AF, DP, alt reads, TLOD, ECNT, alt reads
                      per strand, and whether the position is a TB-Profiler
                      barcode (lineage-defining) site

    collect.py --outdir out --threads 16
"""
import argparse, csv, gzip, os, pickle
from multiprocessing import Pool
import numpy as np

REPO = "/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants"
MASK = f"{REPO}/runroot/refbias/build/7713a8d71d8e-fix1/assets/repeat_mask.bed"
FB = f"{REPO}/analysis/foldback_qc"
BARCODE = "/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/tb-profiler/share/tbprofiler/tbdb.barcode.bed"
PO = "/n/netscratch/sfortune_lab/Lab/pculviner/notebooks/260728_sv_exploration/pipeline_outputs/{0}/results/{0}"
L = 4411532
CLASSES = ("fixed", "mixed_high", "mixed_mid", "mixed_low", "trace")
RES_COLS = ("n_mixed", "mix_freq", "conflict_lineages", "m2_pass", "m2_fail", "m2_mix_call", "m2_low_alt")

mask = np.zeros(L + 200, bool)
for l in open(MASK):
    c, s, e = l.split("\t")[:3]; mask[int(s):int(e)] = True
cmask = np.concatenate([[0], np.cumsum(mask)])
barcode = set()
for l in open(BARCODE):
    f = l.split("\t")
    if len(f) > 2 and f[1].isdigit():
        barcode.add(int(f[2]))  # bed end = 1-based position


def cls(af):
    if af >= 0.9: return "fixed"
    if af >= 0.5: return "mixed_high"
    if af >= 0.25: return "mixed_mid"
    if af >= 0.1: return "mixed_low"
    return "trace"


def one(s):
    try:
        o = pickle.load(gzip.open(PO.format(s) + ".m2variants.pkl.gz"))
    except (OSError, EOFError, pickle.UnpicklingError):
        return s, None, None
    cnt = {f"{t}_{p}_{c}": 0 for t in ("snv", "indel") for p in ("pass", "fail") for c in CLASSES}
    cnt["masked"] = cnt["sv50"] = 0
    rows = []
    for i in range(len(o["POS"])):
        pos, ref, alt = int(o["POS"][i]), o["REF"][i], o["ALT"][i]
        if max(len(ref), len(alt)) > 50:
            cnt["sv50"] += 1; continue
        if cmask[pos - 1 + len(ref)] - cmask[pos - 1] > 0:
            cnt["masked"] += 1; continue
        af = float(o["AF"][i]); ok = bool(o["PASSF"][i])
        t = "snv" if len(ref) == 1 and len(alt) == 1 else "indel"
        c = cls(af)
        cnt[f"{t}_{'pass' if ok else 'fail'}_{c}"] += 1
        if c.startswith("mixed"):
            rows.append((s, pos, ref, alt, t, int(ok), round(af, 4), int(o["DP"][i]), int(o["AD"][i]),
                         float(o["TLOD"][i]), int(o["ECNT"][i]), int(o["ASSB_F"][i]), int(o["ASSB_R"][i]),
                         int(pos in barcode)))
    res = {}
    try:
        r = next(csv.DictReader(open(PO.format(s) + ".results.csv")))
        res = {k: r.get(k, "") for k in RES_COLS}
    except (OSError, StopIteration):
        pass
    return s, (cnt, res), rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    fb = {}
    for f in (f"{FB}/screen50k/rates_all.tsv", f"{FB}/cohorts_in_use_rates.tsv"):
        for r in csv.DictReader(open(f), delimiter="\t"):
            fb[r["sample"]] = r
    qc = [r for r in csv.DictReader(open(f"{FB}/screen50k/pass_standard.tsv"), delimiter="\t") if r["sample"] in fb]
    print(f"isolates: {len(qc)}", flush=True)
    qcols = [k for k in qc[0] if k != "sample"]
    fcols = ["reads", "dup_like", "inv_chimera_per_1k", "other_split_per_1k"]
    with Pool(a.threads) as P, open(f"{a.outdir}/isolates.tsv", "w") as fi, \
            gzip.open(f"{a.outdir}/mixed_calls.tsv.gz", "wt") as fm:
        fm.write("sample\tpos\tref\talt\ttype\tpass\taf\tdp\talt_reads\ttlod\tecnt\talt_fwd\talt_rev\tbarcode\n")
        head = None; missing = 0
        for s, cr, rows in P.imap(one, [r["sample"] for r in qc], chunksize=50):
            if cr is None:
                missing += 1; continue
            cnt, res = cr
            if head is None:
                head = ["sample"] + qcols + fcols + list(RES_COLS) + list(cnt)
                fi.write("\t".join(head) + "\n")
            fi.write("\t".join([s] + [qcd[s][k] for k in qcols] + [fb[s].get(k, "") for k in fcols]
                               + [str(res.get(k, "")) for k in RES_COLS] + [str(cnt[k]) for k in cnt]) + "\n")
            for r in rows:
                fm.write("\t".join(map(str, r)) + "\n")
    print(f"done; {missing} without Mutect2 calls", flush=True)


qcd = {}
if __name__ == "__main__":
    for r in csv.DictReader(open(f"{FB}/screen50k/pass_standard.tsv"), delimiter="\t"):
        qcd[r["sample"]] = r
    main()
