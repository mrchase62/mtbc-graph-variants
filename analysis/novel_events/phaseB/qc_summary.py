#!/usr/bin/env python3
"""One QC row per isolate for the Phase B read QC (08_qc.sbatch), with the
same pass/fail rules as the user's hybrid-assembly pipeline
(mtbc-hybrid-assembly scripts/lib/parse_qc_outputs.py, config/pipeline.yaml).
The parsing is ported, not imported, because that tree is on netscratch.

  contamination  kraken2 on the raw pairs: fail if reads assigned to the
                 M. tuberculosis complex (taxid 77643, clade count) are under
                 85% of classified reads, or if the largest other species or
                 genus (outside the complex and its ancestors) is over 5%
  mixed          TB-Profiler 6.7.0 lineage calls: mixed if two different
                 top-level lineages (text before the first '.') are each at
                 10% or more; fail if no lineage
  depth          after fastp and the 100x downsample: fail below 60x
"""
import argparse
import json

GENOME = 4411532
TARGET = "77643"


def kraken(path):
    rows, stack, anc_t = [], [], set()
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        pct, n_clade, n_direct, rank, taxid, name = f[:6]
        taxid = taxid.strip()
        depth = (len(name) - len(name.lstrip(" "))) // 2
        while stack and stack[-1][0] >= depth:
            stack.pop()
        anc = {t for _, t in stack}
        if taxid == TARGET:
            anc_t = anc
        rows.append(dict(n_clade=int(n_clade), n_direct=int(n_direct), rank=rank.strip(),
                         taxid=taxid, name=name.strip(), in_t=taxid == TARGET or TARGET in anc))
        stack.append((depth, taxid))
    tot = sum(r["n_direct"] for r in rows if r["rank"] != "U")
    t = next((r["n_clade"] for r in rows if r["taxid"] == TARGET), 0)
    oth = sorted((r for r in rows if not r["in_t"] and r["taxid"] not in anc_t
                  and r["rank"] in ("S", "G") and r["n_clade"] > 0), key=lambda r: -r["n_clade"])
    return (t / tot if tot else 0.0, oth[0]["name"] if oth else "none",
            oth[0]["n_clade"] / tot if oth and tot else 0.0)


def tbprofiler(path):
    d = json.load(open(path))
    calls = d.get("lineage", [])
    tops = {c["lineage"].split(".")[0] for c in calls if c.get("fraction", 0) >= 10.0}
    main = d.get("main_lineage", "unknown")
    status = "fail" if not main or main == "unknown" or not calls else \
        "mixed" if len(tops) > 1 else "ok"
    return main, d.get("sub_lineage", "unknown"), status


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--kraken", required=True)
    ap.add_argument("--fastp", required=True, help="fastp JSON")
    ap.add_argument("--tbprofiler", required=True, help="TB-Profiler results JSON")
    ap.add_argument("--final-bases", type=int, required=True, help="bases after downsampling")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    tf, oname, ofrac = kraken(a.kraken)
    fp = json.load(open(a.fastp))
    before, after = fp["summary"]["before_filtering"], fp["summary"]["after_filtering"]
    adapt = fp.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0)
    lin, sub, tbs = tbprofiler(a.tbprofiler)
    depth = a.final_bases / GENOME
    why = []
    if tf < 0.85:
        why.append(f"MTBC {tf:.1%} < 85%")
    if ofrac > 0.05:
        why.append(f"{oname} {ofrac:.1%} > 5%")
    if tbs != "ok":
        why.append(f"TB-Profiler {tbs}")
    if depth < 60:
        why.append(f"depth {depth:.0f}x < 60x")
    cols = dict(sample=a.sample, status="FAIL" if why else "PASS", reason="; ".join(why),
                mtbc_frac=f"{tf:.4f}", top_other=oname, top_other_frac=f"{ofrac:.4f}",
                raw_depth=f"{before['total_bases'] / GENOME:.1f}",
                fastp_depth=f"{after['total_bases'] / GENOME:.1f}",
                reads_raw=before["total_reads"], reads_fastp=after["total_reads"],
                adapter_trimmed_reads=adapt, final_depth=f"{depth:.1f}",
                lineage=lin, sub_lineage=sub, tbprofiler=tbs)
    with open(a.out, "w") as fo:
        fo.write("\t".join(cols) + "\n" + "\t".join(str(v) for v in cols.values()) + "\n")
    print("\t".join(f"{k}={v}" for k, v in cols.items()))


if __name__ == "__main__":
    main()
