#!/usr/bin/env python3
"""Two variant-level scans side by side: which associations one finds and the
other does not.

WHY THIS EXISTS. The merged VCF writes an IS6110 insertion twice where the
accessory catalogue holds a locus that IS the copy: once as CLASS=is6110 and
once as a CLASS=accessory_presence record flagged ACCIS6110. How the two
should be counted against each other -- and against an insertion that sits on
top of an accessory locus, whose carrier set differs -- is not decided (HANDOFF
SV-EVAL, decision a). So cohort_assoc_tail.sh runs the scan on everything and
again with every ACCIS6110 record removed (assoc_scan.py --exclude-acc-is6110),
and this table is the difference: the associations whose call depends on that
undecided choice. Everything significant in both is not listed.

    $MTB_PY assoc/bin/scan_compare.py --full assoc/<c>/scan.tsv \\
        --sens assoc/<c>/scan.no_accis6110.tsv \\
        --out assoc/<c>/scan.is6110_sensitivity.tsv

SIGNIFICANT means assoc_scan.survives(): q < 0.05 under the branch null, the
region null (the level-2 null for a conditional record) and the lineage null.
A row is listed when it survives in one scan and not the other, or when it is
in only one of them (an excluded ACCIS6110 record, which by construction is in
the full scan alone, is listed only if it survived there). Rows are records,
so a test unit of several records appears once per member; the summary counts
both rows and units.
"""
import argparse, csv, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from assoc_scan import survives

FIELDS = ["row_key", "id", "cls", "region", "change", "in_full", "in_sens",
          "survives_full", "survives_sens", "unit_id_full", "unit_id_sens",
          "unit_members_full", "unit_members_sens",
          "gains", "carriers",
          "q_branch_full", "q_branch_sens", "q_region_full", "q_region_sens",
          "q_cond_full", "q_cond_sens", "q_lineage_full", "q_lineage_sens"]


def read_scan(path):
    if not os.path.exists(path):
        sys.exit(f"FATAL: no scan at {path}")
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        need = {"row_key", "q_branch", "q_region", "q_lineage", "q_cond"}
        if not need <= set(rd.fieldnames or []):
            sys.exit(f"FATAL: {path} lacks {sorted(need - set(rd.fieldnames or []))}; "
                     f"not an assoc_scan.py output")
        out = {}
        for r in rd:
            if r["row_key"] in out:
                sys.exit(f"FATAL: {path} repeats the row {r['row_key']}")
            out[r["row_key"]] = r
    return out


def compare(full, sens):
    """(listed rows, counts) for two {row_key: scan row} tables."""
    listed = []
    n = dict(rows_full=len(full), rows_sens=len(sens),
             rows_only_full=0, rows_only_sens=0,
             surv_rows_full=0, surv_rows_sens=0,
             surv_units_full=0, surv_units_sens=0)
    for tag, t in (("full", full), ("sens", sens)):
        for r in t.values():
            if survives(r):
                n[f"surv_rows_{tag}"] += 1
                if r.get("unit_rep", "1") in ("1", ""):
                    n[f"surv_units_{tag}"] += 1
    for k in list(full) + [k for k in sens if k not in full]:
        f, s = full.get(k), sens.get(k)
        sf = survives(f) if f else False
        ss = survives(s) if s else False
        if f is None:
            n["rows_only_sens"] += 1
        if s is None:
            n["rows_only_full"] += 1
        if sf == ss:
            continue
        if s is None:
            change = "significant_in_full_only_excluded_from_sens"
        elif f is None:
            change = "significant_in_sens_only_not_in_full"
        else:
            change = ("significant_in_full_only" if sf
                      else "significant_in_sens_only")
        b = f or s
        o = dict(row_key=k, id=b["id"], cls=b["cls"], region=b["region"],
                 change=change, in_full=int(f is not None),
                 in_sens=int(s is not None), survives_full=int(sf),
                 survives_sens=int(ss), gains=b.get("gains", ""),
                 carriers=b.get("carriers", ""))
        for tag, r in (("full", f), ("sens", s)):
            o[f"unit_id_{tag}"] = (r or {}).get("unit_id", "")
            o[f"unit_members_{tag}"] = (r or {}).get("unit_members", "")
            for q in ("q_branch", "q_region", "q_cond", "q_lineage"):
                o[f"{q}_{tag}"] = (r or {}).get(q, "")
        listed.append(o)
    listed.sort(key=lambda o: (o["change"], o["cls"], o["row_key"]))
    return listed, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", required=True, help="the scan on every record")
    ap.add_argument("--sens", required=True,
                    help="the scan with the ACCIS6110 records excluded")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    full, sens = read_scan(a.full), read_scan(a.sens)
    # The sensitivity scan is the full one MINUS records; a record in it that
    # the full scan does not hold means the two were run on different inputs.
    extra = [k for k in sens if k not in full]
    if extra:
        sys.exit(f"FATAL: {len(extra):,} rows of {a.sens} are not in {a.full}, "
                 f"e.g. {extra[:3]}; the two scans are not of one event matrix")
    listed, n = compare(full, sens)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=FIELDS)
        w.writeheader()
        w.writerows(listed)
    by = {}
    for o in listed:
        by[o["change"]] = by.get(o["change"], 0) + 1
    lines = [f"full scan         {a.full}",
             f"sensitivity scan  {a.sens}",
             f"rows              {n['rows_full']:,} full, {n['rows_sens']:,} "
             f"sensitivity ({n['rows_only_full']:,} only in full: the "
             f"excluded ACCIS6110 records that were tested)",
             f"surviving         {n['surv_units_full']:,} units "
             f"({n['surv_rows_full']:,} rows) full; {n['surv_units_sens']:,} "
             f"units ({n['surv_rows_sens']:,} rows) sensitivity",
             f"differ            {len(listed):,} rows"]
    lines += [f"  {k:<46}{v:>6,}" for k, v in sorted(by.items())]
    txt = "\n".join(lines) + "\n"
    stem = a.out[:-4] if a.out.endswith(".tsv") else a.out
    with open(f"{stem}.summary.txt", "w") as fh:
        fh.write(txt)
    print(txt, end="")
    print(f"  -> {a.out}, {stem}.summary.txt")


if __name__ == "__main__":
    main()
