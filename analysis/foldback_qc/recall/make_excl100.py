#!/usr/bin/env python3
"""Excluded-isolates check: runner tables for cohort fbrecallx100, the 12
excluded selfcal54 isolates (v2 > 4.9) with every fold-back read pair removed
(04_excl_sets.sbatch), each on the same reference as its all-reads sample in
fbrecall. Same library, reference and truth; the only difference is the
fold-back reads (as make_fb100.py for the chimeric controls). Run from
runroot with $MTB_PY. --list prints the 12 isolates.
Writes refbias/fbrecallx100/{cohort.tsv,crams.tsv,run/p1/refmap.tsv,run/.mtb_build}.
"""
import csv
import datetime
import os
import sys

BUILD = "7713a8d71d8e-fix1"
C = "refbias/fbrecallx100"
CUT = 4.9


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def excluded():
    rate = {p["sample"]: float(p["rate"]) for p in tsv("../analysis/foldback_qc/selfcal/out/calib_points.tsv")}
    return [r["sample"] for r in tsv("refbias/fbrecall/run/p1/refmap.tsv")
            if r["sample"].startswith("ps_") and rate[r["sample"]] > CUT]


def main():
    ex = excluded()
    if "--list" in sys.argv:
        print("\n".join(ex))
        return
    rm = {r["sample"]: r for r in tsv("refbias/fbrecall/run/p1/refmap.tsv")}
    co = {r["sample"]: r for r in tsv("refbias/fbrecall/cohort.tsv")}
    root = os.environ["MTB_CRAM_ROOT"]
    cols = list(next(iter(rm.values())))
    os.makedirs(f"{C}/run/p1", exist_ok=True)
    with open(f"{C}/run/.mtb_build", "w") as fo:
        fo.write(f"build_id\t{BUILD}\nrecorded\t{datetime.datetime.now().astimezone().isoformat(timespec='seconds')}\n")
    with open(f"{C}/cohort.tsv", "w") as fc, open(f"{C}/crams.tsv", "w") as fr, \
            open(f"{C}/run/p1/refmap.tsv", "w") as fm:
        fc.write("\t".join(co[next(iter(co))].keys()) + "\n")
        fr.write("sample\tcram_relpath\n")
        fm.write("\t".join(cols) + "\n")
        for s in ex:
            bam = f"{C}/reads/{s}_fb100.bam"
            assert os.path.exists(bam), bam
            t = f"{s}_fb100"
            fc.write("\t".join([t] + list(co[s].values())[1:]) + "\n")
            fr.write(f"{t}\t{os.path.relpath(os.path.realpath(bam), root)}\n")
            fm.write("\t".join(dict(rm[s], sample=t)[c] for c in cols) + "\n")
            print(f"{t}\tR {rm[s]['reference']}")
    print(f"{len(ex)} fb100 samples")


if __name__ == "__main__":
    main()
