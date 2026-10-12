#!/usr/bin/env python3
"""Part D, causal arm: runner tables for cohort fbrecall100, the 9 chimeric
controls with every fold-back read pair removed (../removal, set fb100),
each mapped to the same reference R as its all-reads sample in fbrecall
(make_cohort.py). Same library, reference and truth; the only difference is
the fold-back reads. Run from runroot with $MTB_PY after make_cohort.py.
Writes refbias/fbrecall100/{cohort.tsv,crams.tsv,run/p1/refmap.tsv,run/.mtb_build}.
"""
import csv
import datetime
import os

BUILD = "7713a8d71d8e-fix1"
C = "refbias/fbrecall100"


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def main():
    rm = {r["sample"]: r for r in tsv("refbias/fbrecall/run/p1/refmap.tsv")}
    co = {r["sample"]: r for r in tsv("refbias/fbrecall/cohort.tsv")}
    root = os.environ["MTB_CRAM_ROOT"]
    cols = list(next(iter(rm.values())))
    os.makedirs(f"{C}/run/p1", exist_ok=True)
    with open(f"{C}/run/.mtb_build", "w") as fo:
        fo.write(f"build_id\t{BUILD}\nrecorded\t{datetime.datetime.now().astimezone().isoformat(timespec='seconds')}\n")
    n = 0
    with open(f"{C}/cohort.tsv", "w") as fc, open(f"{C}/crams.tsv", "w") as fr, \
            open(f"{C}/run/p1/refmap.tsv", "w") as fm:
        fc.write("\t".join(co[next(iter(co))].keys()) + "\n")
        fr.write("sample\tcram_relpath\n")
        fm.write("\t".join(cols) + "\n")
        for s in (x.split()[0] for x in open("../analysis/foldback_qc/removal/isolates.txt") if x.strip()):
            bam = f"refbias/fbremove/reads/{s}_fb100.bam"
            assert os.path.exists(bam), bam
            t = f"{s}_fb100"
            fc.write("\t".join([t] + list(co[s].values())[1:]) + "\n")
            fr.write(f"{t}\t{os.path.relpath(os.path.realpath(bam), root)}\n")
            m = dict(rm[s], sample=t)
            fm.write("\t".join(m[c] for c in cols) + "\n")
            n += 1
            print(f"{t}\tR {m['reference']}")
    print(f"{n} fb100 samples")


if __name__ == "__main__":
    main()
