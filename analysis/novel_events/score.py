#!/usr/bin/env python3
"""Score structural-event calls against the Phase A truth set.

Callers: the prototype breakpoint caller (breakpoint_caller.py tables), and
dysgu and delly from P2 (VCFs in the matched reference's coordinates), each
scored with all calls and with PASS calls only.

A call matches a true event when:
  - the start is within 50 bp of the true start, or, for deletions, the end
    is within 50 bp of the true end, or the two overlap reciprocally by 50%
    or more;
  - and the types are compatible:
      DEL   DEL
      INS   INS
      REPL  DEL or INS
      INV   INV
      REARR REARR, INV or DEL
Typed recall counts only compatible matches. Breakpoint recall also accepts
any call within 50 bp, including the prototype's untyped BND calls.
Precision: typed calls (BND excluded) that match any true event within
50 bp, over all typed calls.

Truth events are split by size (50-150, 150-500, 500-2,000, >2,000 bp) and
by whether an insertion is IS6110-sized (1,340-1,370 bp).
"""
import argparse
import collections
import csv
import os

TOL = 50
COMPAT = {"DEL": {"DEL"}, "INS": {"INS"}, "REPL": {"DEL", "INS"},
          "INV": {"INV"}, "REARR": {"REARR", "INV", "DEL"}}


def size_class(n):
    return "50-150" if n <= 150 else "150-500" if n <= 500 else \
        "500-2k" if n <= 2000 else ">2k"


def read_vcf(path, pass_only):
    out = []
    if not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if pass_only and c[6] not in ("PASS", "."):
            continue
        info = dict(kv.split("=", 1) if "=" in kv else (kv, "") for kv in c[7].split(";"))
        t = info.get("SVTYPE", "")
        if t == "DUP":
            t = "INS"
        elif t in ("BND", "TRA"):
            t = "REARR"
        pos = int(c[1])
        end = int(info.get("END", pos)) if info.get("END", "").lstrip("-").isdigit() else pos
        out.append(dict(type=t, start=pos, end=max(pos, end)))
    return out


def read_proto(path):
    if not os.path.exists(path):
        return []
    return [dict(type=r["type"], start=int(r["start"]), end=int(r["end"]))
            for r in csv.DictReader(open(path), delimiter="\t")]


def near(c, t):
    if abs(c["start"] - t["start"]) <= TOL:
        return True
    if t["type"] in ("DEL", "REPL") and c["type"] == "DEL":
        if abs(c["end"] - t["end"]) <= TOL:
            return True
        ov = min(c["end"], t["end"]) - max(c["start"], t["start"]) + 1
        if ov > 0 and ov >= 0.5 * (t["end"] - t["start"] + 1) \
                and ov >= 0.5 * (c["end"] - c["start"] + 1):
            return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--proto-dir", required=True)
    ap.add_argument("--p2-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    truth = collections.defaultdict(list)
    for r in csv.DictReader(open(a.truth), delimiter="\t"):
        size = max(int(r["ref_len"]), int(r["alt_len"]))
        truth[r["sample"]].append(dict(type=r["type"], start=int(r["r_start"]),
                                       end=int(r["r_end"]), size=size,
                                       is6110=r["type"] in ("INS", "REPL") and 1340 <= int(r["alt_len"]) <= 1370))
    callers = {
        "prototype": lambda s: read_proto(os.path.join(a.proto_dir, f"{s}.events.tsv")),
        "dysgu_all": lambda s: read_vcf(os.path.join(a.p2_dir, f"{s}.dysgu.vcf"), False),
        "dysgu_pass": lambda s: read_vcf(os.path.join(a.p2_dir, f"{s}.dysgu.vcf"), True),
        "delly_all": lambda s: read_vcf(os.path.join(a.p2_dir, f"{s}.delly.vcf"), False),
        "delly_pass": lambda s: read_vcf(os.path.join(a.p2_dir, f"{s}.delly.vcf"), True),
    }
    calls = {k: {s: f(s) for s in truth} for k, f in callers.items()}
    calls["dysgu+delly_pass"] = {s: calls["dysgu_pass"][s] + calls["delly_pass"][s] for s in truth}
    calls["prototype+dysgu_pass"] = {s: calls["prototype"][s] + calls["dysgu_pass"][s] for s in truth}

    rows = []
    for name, by_s in calls.items():
        strata = collections.defaultdict(lambda: [0, 0, 0])  # n, typed hit, bp hit
        n_calls = n_true_calls = 0
        for s, ts in truth.items():
            cs = by_s.get(s, [])
            for t in ts:
                keys = ["all", f"type={t['type']}", f"size={size_class(t['size'])}"]
                if t["is6110"]:
                    keys.append("IS6110-sized INS")
                typed = any(c["type"] in COMPAT[t["type"]] and near(c, t) for c in cs)
                bp = any(near(c, t) or abs(c["start"] - t["start"]) <= TOL for c in cs)
                for k in keys:
                    strata[k][0] += 1
                    strata[k][1] += typed
                    strata[k][2] += bp
            for c in cs:
                if c["type"] == "BND":
                    continue
                n_calls += 1
                n_true_calls += any(near(c, t) or abs(c["start"] - t["start"]) <= TOL for t in ts)
        prec = n_true_calls / n_calls if n_calls else float("nan")
        for k, (n, h, b) in sorted(strata.items()):
            rows.append(dict(caller=name, stratum=k, n_true=n, typed_recall=round(h / n, 3),
                             breakpoint_recall=round(b / n, 3), typed_calls=n_calls,
                             precision=round(prec, 3)))
    with open(a.out, "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        if r["stratum"] == "all":
            print(r)


if __name__ == "__main__":
    main()
