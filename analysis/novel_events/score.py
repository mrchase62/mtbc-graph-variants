#!/usr/bin/env python3
"""Score structural-event calls against the Phase A truth set.

Callers: the prototype breakpoint caller (breakpoint_caller.py tables), and
dysgu from P2 (VCF in the matched reference's coordinates), scored with all
calls and with PASS calls only. delly was dropped 2026-10-10: typed recall
0.10 on the Phase B novel events, adding nothing to our caller + dysgu.

A call matches a true event when:
  - the start is within 50 bp of the true start, or, for deletions, the end
    is within 50 bp of the true end, or the two overlap reciprocally by 50%
    or more;
  - or (added 2026-10-10, truvari bench's defaults: refdist 500, pctsize
    0.7) the starts are within 500 bp and the smaller size is at least 70%
    of the larger, when both sizes are known. The same insertion or deletion
    in a repeat can be placed anywhere along it, and was scored false 70-120
    bp from its true copy. truvari's sequence-similarity test (pctseq) is
    not applied: the truth set and the callers' tables carry no allele
    sequences;
  - and the types are compatible:
      DEL   DEL
      INS   INS
      REPL  DEL or INS
      INV   INV
      REARR REARR, INV or DEL
Typed recall counts only compatible matches. Breakpoint recall also accepts
any call within 50 bp of the true start or (added 2026-10-10) the true end,
including the prototype's untyped BND calls.
Precision: typed calls (BND excluded) that match a true event (above) or lie
within 50 bp of its start or end, over all typed calls. A call at the far end
of a true deletion was counted false before.

dysgu_af20 (the dysgu calls used in the combination from 2026-10-10): PASS
and FORMAT AF >= 0.2. AF is the share of reads at the site with direct
evidence (split, clipped or discordant); true events in these clonal
isolates sit near 0.35, most false calls below 0.1 (phaseB/dysgu_false.py).

Depth-scan calls (step 2) match when the true event lies within the call's
run +-200 bp, since read-depth boundaries are only window-precise.

Truth events are split by size (50-150, 150-500, 500-2,000, >2,000 bp) and
by whether an insertion is IS6110-sized (1,340-1,370 bp).

GRIDSS (--gridss-dir, <dir>/<sample>/<sample>.gridss.vcf): breakend pairs
are typed by the usual reading of VCF breakend notation, once per pair:
  t[p[ with p > pos, or ]p]t with p < pos   deletion of the bases between,
      or an insertion when the inserted sequence is at least 70% as long
      as the deleted span (a replacement is scored as INS or DEL either way)
  t[p[ with p < pos, or ]p]t with p > pos   tandem duplication of the segment,
      scored as an insertion of its length (as the prototype does)
  t]p] or [p[t                            inversion
  other chromosome                       REARR
Single breakends (t. or .t; the partner could not be placed) are untyped
BND, counted for breakpoint recall only. gridss_pass uses FILTER=PASS;
gridss_all uses every record.

Callable sequence (--mask-dir with --refmap; Phase B): a true event or a call
whose span +-50 bp touches the matched reference's uncallable mask
(callable_mask.py, <reference>.bed; plus artifact_mask.py's
<reference>.artifact.bed with --artifact-mask, off by default: the mask is
under evaluation, user 2026-10-10) is uncallable. Recall is then measured over callable true
events only, and uncallable calls are flagged, not counted in precision.
n_true_uncallable and calls_uncallable report how many were set aside.
"""
import argparse
import collections
import csv
import bisect
import os

TOL = 50
COMPAT = {"DEL": {"DEL"}, "INS": {"INS"}, "REPL": {"DEL", "INS"},
          "INV": {"INV"}, "REARR": {"REARR", "INV", "DEL"}}


def size_class(n):
    return "50-150" if n <= 150 else "150-500" if n <= 500 else \
        "500-2k" if n <= 2000 else ">2k"


def read_vcf(path, pass_only, min_af=None):
    out = []
    if not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if pass_only and c[6] not in ("PASS", "."):
            continue
        if min_af is not None:
            fmt = dict(zip(c[8].split(":"), c[9].split(":"))) if len(c) > 9 else {}
            try:
                if float(fmt.get("AF", "nan")) < min_af:
                    continue
            except ValueError:
                pass
        info = dict(kv.split("=", 1) if "=" in kv else (kv, "") for kv in c[7].split(";"))
        t = info.get("SVTYPE", "")
        if t == "DUP":
            t = "INS"
        elif t in ("BND", "TRA"):
            t = "REARR"
        pos = int(c[1])
        end = int(info.get("END", pos)) if info.get("END", "").lstrip("-").isdigit() else pos
        svlen = info.get("SVLEN", "").lstrip("-")
        size = int(svlen) if svlen.isdigit() and int(svlen) > 0 else \
            (max(pos, end) - pos + 1 if t in ("DEL", "INV") and end > pos else None)
        out.append(dict(type=t, start=pos, end=max(pos, end), size=size))
    return out


def read_proto(path, depth=False):
    if not os.path.exists(path):
        return []
    return [dict(type=r["type"], start=int(r["start"]), end=int(r["end"]), depth=depth,
                 size=int(r["size"]) if r.get("size", "").isdigit() and int(r["size"]) > 0 else None)
            for r in csv.DictReader(open(path), delimiter="\t")]


def near(c, t):
    if c.get("depth"):  # depth runs have fuzzy ends: the event within the run +-200 bp
        return c["start"] - 200 <= t["start"] <= c["end"] + 200
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


SIZE_DIST, PCTSIZE = 500, 0.7  # truvari bench defaults (refdist, pctsize)


def matches(c, t):
    """near(), or the same event placed elsewhere in a repeat: starts within
    500 bp and sizes within 70% of each other (both sizes known)."""
    if near(c, t):
        return True
    if c.get("depth") or not c.get("size") or not t.get("size") or t["type"] == "REARR":
        return False
    return abs(c["start"] - t["start"]) <= SIZE_DIST and \
        min(c["size"], t["size"]) >= PCTSIZE * max(c["size"], t["size"])


def at_event(c, t):
    """For precision and breakpoint recall: a match, or within 50 bp of the
    true start or end."""
    return matches(c, t) or abs(c["start"] - t["start"]) <= TOL or abs(c["start"] - t["end"]) <= TOL


def load_masks(refmap, mask_dir, artifacts=False):
    ref = {}
    for r in csv.DictReader(open(refmap), delimiter="\t"):
        ref[r["sample"]] = r["reference"]
    masks = {}
    for s, g in ref.items():
        beds = [os.path.join(mask_dir, f"{g}.bed")]
        art = os.path.join(mask_dir, f"{g}.artifact.bed")  # artifact_mask.py; opt-in
        if artifacts and os.path.exists(art):
            beds.append(art)
        iv = sorted((int(f[1]), int(f[2])) for p in beds for f in (x.split("\t") for x in open(p)))
        masks[s] = ([x[0] for x in iv], iv)
    return masks


def masked(m, s, e):
    if m is None:
        return False
    starts, iv = m
    lo, hi = s - 1 - TOL, e + TOL  # 1-based event span +-50 bp, as 0-based half-open
    i = bisect.bisect_right(starts, hi)
    return any(iv[j][1] > lo for j in range(max(0, i - 8), i))


def read_gridss(path, pass_only):
    import re
    if not os.path.exists(path):
        return []
    recs = {}
    out = []
    for line in open(path):
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if pass_only and c[6] != "PASS":
            continue
        pos, alt = int(c[1]), c[4]
        info = dict(x.split("=", 1) if "=" in x else (x, "") for x in c[7].split(";"))
        if alt.startswith(".") or alt.endswith("."):
            out.append(dict(type="BND", start=pos, end=pos))
            continue
        m = re.match(r"^([A-Za-z]*)([\[\]])([^:\[\]]+):(\d+)([\[\]])([A-Za-z]*)$", alt)
        if not m:
            continue
        mate = info.get("MATEID", "")
        if mate in recs:  # the pair is typed from its first record
            continue
        recs[c[2]] = 1
        before, br, chrom, p, _, after = m.groups()
        p = int(p)
        ins = max(len(before), len(after)) - 1
        if chrom != c[0]:
            out.append(dict(type="REARR", start=pos, end=pos))
        elif (before and br == "[") and p > pos or (after and br == "]") and p < pos:
            a, b = min(pos, p), max(pos, p)
            span = b - a - 1
            if ins >= 0.7 * span:
                out.append(dict(type="INS", start=a, end=a, size=ins))
            else:
                out.append(dict(type="DEL", start=a + 1, end=b - 1, size=span))
        elif (before and br == "[") or (after and br == "]"):
            out.append(dict(type="INS", start=max(pos, p), end=max(pos, p), size=abs(pos - p) + 1))
        else:
            out.append(dict(type="INV", start=min(pos, p), end=max(pos, p)))
    # events under 50 bp are out of scope, as in the truth set
    return [x for x in out if x.get("size", 50) >= 50]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--proto-dir", required=True)
    ap.add_argument("--p2-dir", required=True)
    ap.add_argument("--asm-dir", help="assemble_candidates.py tables (step 1)")
    ap.add_argument("--depth-dir", help="depth_scan.py tables (step 2)")
    ap.add_argument("--gridss-dir", help="GRIDSS output folder, <sample>/<sample>.gridss.vcf")
    ap.add_argument("--refmap", help="P1 refmap.tsv: sample -> matched reference")
    ap.add_argument("--mask-dir", help="callable_mask.py BEDs, <reference>.bed")
    ap.add_argument("--artifact-mask", action="store_true",
                    help="also mask artifact_mask.py's sites (under evaluation)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    masks = load_masks(a.refmap, a.mask_dir, a.artifact_mask) if a.mask_dir else {}
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
        "dysgu_af20": lambda s: read_vcf(os.path.join(a.p2_dir, f"{s}.dysgu.vcf"), True, 0.2),
    }
    calls = {k: {s: f(s) for s in truth} for k, f in callers.items()}
    if a.gridss_dir:
        for nm, ps in (("gridss_pass", True), ("gridss_all", False)):
            calls[nm] = {s: read_gridss(os.path.join(a.gridss_dir, s, f"{s}.gridss.vcf"), ps)
                         for s in truth}
    calls["prototype+dysgu_pass"] = {s: calls["prototype"][s] + calls["dysgu_pass"][s] for s in truth}
    if a.asm_dir:
        calls["proto_asm"] = {s: read_proto(os.path.join(a.asm_dir, f"{s}.events.tsv")) for s in truth}
        calls["proto_asm+dysgu_pass"] = {s: calls["proto_asm"][s] + calls["dysgu_pass"][s] for s in truth}
        calls["proto_asm+dysgu_af20"] = {s: calls["proto_asm"][s] + calls["dysgu_af20"][s] for s in truth}
    if a.depth_dir:
        calls["depth"] = {s: read_proto(os.path.join(a.depth_dir, f"{s}.events.tsv"), True) for s in truth}
    if a.asm_dir and a.depth_dir:
        calls["proto_asm+depth"] = {s: calls["proto_asm"][s] + calls["depth"][s] for s in truth}
        calls["proto_asm+depth+dysgu_pass"] = {s: calls["proto_asm+depth"][s] + calls["dysgu_pass"][s]
                                               for s in truth}

    rows = []
    n_unc_truth = collections.Counter()
    for s, ts in truth.items():
        m = masks.get(s)
        for t in ts:
            t["unc"] = masked(m, t["start"], t["end"])
            if t["unc"]:
                n_unc_truth["all"] += 1
    for name, by_s in calls.items():
        strata = collections.defaultdict(lambda: [0, 0, 0])  # n, typed hit, bp hit
        n_calls = n_true_calls = n_unc_calls = 0
        for s, ts in truth.items():
            m = masks.get(s)
            cs_all = by_s.get(s, [])
            cs = [c for c in cs_all if not masked(m, c["start"], c["end"])]
            n_unc_calls += sum(1 for c in cs_all if c["type"] != "BND") - \
                sum(1 for c in cs if c["type"] != "BND")
            for t in ts:
                if t["unc"]:
                    continue
                keys = ["all", f"type={t['type']}", f"size={size_class(t['size'])}"]
                if t["is6110"]:
                    keys.append("IS6110-sized INS")
                typed = any(c["type"] in COMPAT[t["type"]] and matches(c, t) for c in cs)
                bp = any(at_event(c, t) for c in cs)
                for k in keys:
                    strata[k][0] += 1
                    strata[k][1] += typed
                    strata[k][2] += bp
            for c in cs:
                if c["type"] == "BND":
                    continue
                n_calls += 1
                n_true_calls += any(at_event(c, t) for t in ts)
        prec = n_true_calls / n_calls if n_calls else float("nan")
        for k, (n, h, b) in sorted(strata.items()):
            rows.append(dict(caller=name, stratum=k, n_true=n, typed_recall=round(h / n, 3),
                             breakpoint_recall=round(b / n, 3), typed_calls=n_calls,
                             precision=round(prec, 3),
                             n_true_uncallable=n_unc_truth[k] if k == "all" else "",
                             calls_uncallable=n_unc_calls))
    with open(a.out, "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        if r["stratum"] == "all":
            print(r)


if __name__ == "__main__":
    main()
