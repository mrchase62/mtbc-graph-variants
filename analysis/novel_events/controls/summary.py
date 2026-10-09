#!/usr/bin/env python3
"""One row per in-panel control: false calls (false_calls.tsv), fold-back
chimeras and the duplicate proxy (01_artifacts), oxidative damage (GATK
pre-adapter G>T and deamination C>T, Phred-scaled, lower is worse), Picard
PCT_CHIMERAS and PCT_ADAPTER (02_picard), and kraken2 contamination (phaseB
09_kraken, parsed as qc_summary.py does). Missing inputs are left blank."""
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phaseB"))
from qc_summary import kraken  # noqa: E402

D = os.path.dirname(os.path.abspath(__file__))
O = os.path.join(D, "out")
KR = os.path.join(D, "..", "..", "..", "runroot", "refbias", "marinB63", "qc")


def table(path, start):
    if not os.path.exists(path):
        return []
    lines = [x for x in open(path) if x.strip() and not x.startswith("#")]
    i = next((j for j, x in enumerate(lines) if x.startswith(start)), None)
    return [] if i is None else list(csv.DictReader(lines[i:], delimiter="\t"))


rows = []
for fc in csv.DictReader(open(os.path.join(D, "false_calls.tsv")), delimiter="\t"):
    s = fc["sample"]
    r = dict(sample=s, ours_false=fc["ours_asm"], dysgu_false=fc["dysgu_pass"])
    sp = table(os.path.join(O, s, "split_reads.tsv"), "sample")
    r.update(inv_chimera_per_1k=sp[0]["inv_chimera_per_1k"] if sp else "",
             other_split_per_1k=sp[0]["other_split_per_1k"] if sp else "",
             dup_like=sp[0]["dup_like"] if sp else "")
    pa = {x["REF_BASE"] + x["ALT_BASE"]: x["TOTAL_QSCORE"]
          for x in table(os.path.join(O, s, "artifacts.pre_adapter_summary_metrics.txt"), "SAMPLE_ALIAS")}
    r.update(oxoG_GT_q=pa.get("GT", ""), deam_CT_q=pa.get("CT", ""))
    al = {x["CATEGORY"]: x for x in table(os.path.join(O, s, "alignment_summary.txt"), "CATEGORY")}
    p = al.get("PAIR", {})
    r.update(picard_pct_chimeras=p.get("PCT_CHIMERAS", ""), picard_pct_adapter=p.get("PCT_ADAPTER", ""))
    kp = os.path.join(KR, s, "kraken2.tsv")
    if os.path.exists(kp) and os.path.getsize(kp):
        tf, on, of = kraken(kp)
        r.update(mtbc_frac=f"{tf:.4f}", top_other=on, top_other_frac=f"{of:.4f}")
    else:
        r.update(mtbc_frac="", top_other="not run", top_other_frac="")
    rows.append(r)
with open(os.path.join(D, "summary.tsv"), "w") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
