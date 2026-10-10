#!/usr/bin/env python3
"""What are dysgu's remaining false calls on clean reads, and can a filter on
dysgu's own fields remove them without losing true events? (run from runroot
after rerun_compare.py; report only, no compute)

Isolates: the 46 compared in rerun_compare.py (marinQC51, P1 reference
unchanged). Calls: dysgu PASS, typed (BND excluded), in callable sequence,
as in score.py.

Each call is labelled
  true    within 50 bp of a true event of its isolate (score.py's precision rule)
  false   otherwise
and carries dysgu's evidence fields: AF (share of reads at the site that
support the event), SU (supporting reads), PE / SR / SC (pairs, split reads,
soft clips), PROB (dysgu's model probability), MAPQP, COV, plus
`recurrent`: how many other isolates on the same matched reference have a
dysgu PASS call within 50 bp.

Then single-field cutoffs are swept. For each: false calls removed, true
calls lost, and typed recall / precision of dysgu alone and of our caller +
assembly + filtered dysgu, per group, with score.py's matching rules.

Writes out/qc51/dysgu_calls.tsv and out/qc51/dysgu_filters.tsv.
"""
import collections
import csv
import os
import sys

N = "../analysis/novel_events"
sys.path.insert(0, N)
from score import COMPAT, TOL, load_masks, masked, near, read_proto  # noqa: E402

OUT = f"{N}/phaseB/out/qc51"
P2 = "refbias/marinQC51/run/p2"
RM = "refbias/marinQC51/run/p1/refmap.tsv"
RM_BEFORE = "refbias/marinB63/run/p1/refmap.tsv"
FIELDS = ("AF", "SU", "PE", "SR", "SC", "PROB", "MAPQP", "COV")
CUTS = dict(AF=(0.1, 0.2, 0.3, 0.5, 0.7), SU=(4, 6, 8, 10, 15), PROB=(0.5, 0.7, 0.8, 0.9),
            MAPQP=(30, 50, 59), recurrent_max=(0, 1, 2))


def read_dysgu(path):
    out = []
    for line in open(path):
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if c[6] not in ("PASS", "."):
            continue
        info = dict(kv.split("=", 1) if "=" in kv else (kv, "") for kv in c[7].split(";"))
        fmt = dict(zip(c[8].split(":"), c[9].split(":"))) if len(c) > 9 else {}
        t = {"DUP": "INS", "BND": "REARR", "TRA": "REARR"}.get(info.get("SVTYPE", ""), info.get("SVTYPE", ""))
        pos = int(c[1])
        end = int(info["END"]) if info.get("END", "").lstrip("-").isdigit() else pos
        d = dict(type=t, start=pos, end=max(pos, end), svlen=abs(int(info.get("SVLEN", 0) or 0)),
                 kind=info.get("KIND", ""), ct=info.get("CT", ""))
        for f in FIELDS:
            try:
                d[f] = float(fmt.get(f, info.get(f, "nan")))
            except ValueError:
                d[f] = float("nan")
        out.append(d)
    return out


def main():
    groups = {r["sample"]: r["group"] for r in csv.DictReader(open(f"{N}/phaseB/out/samples.tsv"), delimiter="\t")}
    ref = {r["sample"]: r["reference"] for r in csv.DictReader(open(RM), delimiter="\t")}
    ref_b = {r["sample"]: r["reference"] for r in csv.DictReader(open(RM_BEFORE), delimiter="\t")}
    keep = [s for s in ref if ref[s] == ref_b.get(s)]
    masks = load_masks(RM, f"{N}/out/masks")
    truth = collections.defaultdict(list)
    for g in ("novel", "setE_polished", "in_panel"):
        for r in csv.DictReader(open(f"{OUT}/truth.{g}.tsv"), delimiter="\t"):
            if r["sample"] in keep:
                truth[r["sample"]].append(dict(type=r["type"], start=int(r["r_start"]), end=int(r["r_end"])))
    # as in score.py: a call is true if near any true event; recall counts callable true events only
    truth_all = {s: list(truth[s]) for s in keep}
    for s in keep:
        truth[s] = [t for t in truth[s] if not masked(masks[s], t["start"], t["end"])]

    def callable_typed(s, cs):
        return [c for c in cs if c["type"] != "BND" and not masked(masks[s], c["start"], c["end"])]

    dys = {s: callable_typed(s, read_dysgu(f"{P2}/{s}.dysgu.vcf")) for s in keep}
    ours = {s: callable_typed(s, read_proto(f"{N}/phaseB/out/qc51/asm/{s}.events.tsv")) for s in keep}
    rows = []
    for s in keep:
        for c in dys[s]:
            c["true"] = any(near(c, t) or abs(c["start"] - t["start"]) <= TOL for t in truth_all[s])
            c["recurrent"] = sum(1 for o in keep if o != s and ref[o] == ref[s]
                                 and any(abs(x["start"] - c["start"]) <= TOL for x in dys[o]))
            c["ours_also"] = any(abs(x["start"] - c["start"]) <= TOL for x in ours[s])
            rows.append(dict(sample=s, group=groups[s], reference=ref[s], label="true" if c["true"] else "false",
                             **{k: c[k] for k in ("type", "start", "end", "svlen", "kind", "ct") + FIELDS
                                + ("recurrent", "ours_also")}))
    with open(f"{OUT}/dysgu_calls.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    def score(keep_call):
        res = []
        for g in ("novel", "setE_polished", "in_panel"):
            iso = [s for s in keep if groups[s] == g]
            fd = {s: [c for c in dys[s] if keep_call(c)] for s in iso}
            n_true = sum(len(truth[s]) for s in iso)
            hit_d = sum(any(c["type"] in COMPAT[t["type"]] and near(c, t) for c in fd[s]) for s in iso for t in truth[s])
            hit_c = sum(any(c["type"] in COMPAT[t["type"]] and near(c, t) for c in fd[s] + ours[s])
                        for s in iso for t in truth[s])
            tp = sum(c["true"] for s in iso for c in fd[s])
            fp = sum(not c["true"] for s in iso for c in fd[s])
            ctp = sum(any(near(c, t) or abs(c["start"] - t["start"]) <= TOL for t in truth_all[s])
                      for s in iso for c in ours[s] + fd[s])
            cn = sum(len(ours[s]) + len(fd[s]) for s in iso)
            res.append(dict(group=g, isolates=len(iso), n_true=n_true, dysgu_true_calls=tp, dysgu_false_calls=fp,
                            dysgu_recall=round(hit_d / n_true, 3) if n_true else "",
                            dysgu_precision=round(tp / (tp + fp), 3) if tp + fp else "",
                            combo_recall=round(hit_c / n_true, 3) if n_true else "",
                            combo_precision=round(ctp / cn, 3) if cn else "",
                            false_per_isolate=",".join(str(sum(not c["true"] for c in fd[s])) for s in iso)))
        return res

    filt = [("none", "", lambda c: True)]
    for f, cuts in CUTS.items():
        for x in cuts:
            if f == "recurrent_max":
                filt.append((f, x, lambda c, x=x: c["recurrent"] <= x))
            else:
                filt.append((f, x, lambda c, f=f, x=x: c[f] >= x))
    out = []
    for f, x, fn in filt:
        for r in score(fn):
            out.append(dict(filter=f, cutoff=x, **r))
    with open(f"{OUT}/dysgu_filters.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)


if __name__ == "__main__":
    main()
