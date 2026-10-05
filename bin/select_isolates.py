#!/usr/bin/env python3
"""Standing inclusion screen for real isolates, from the collection's results table.

Public SRA data varies enormously in quality, and a marginal isolate does not
give a marginally worse result -- it gives a confidently wrong one. A low-depth
sample produces missing calls that score as false negatives attributable to the
method rather than the data. A mixed sample carries real minority alleles at
positions a single-strain analysis scores as errors. Neither is visible in the
output, so the screen has to happen at selection time.

Criteria (project standing rule):

  paired-end       `is_paired`. Single-end libraries are excluded outright.
  depth >= 60x     `meandepth`. NOT `coverage`, which is percent of genome
                   covered (typically ~99) and is a different quantity.
  not mixed        `conflict_lineages` / `n_mixed` / `mix_freq`. Excluded at any
                   minor-component frequency. These flags come from the table and
                   cannot be re-derived: a mixture at a few percent is invisible
                   to a lineage call, which collapses to the majority.
  complete run     `exit_type`.
  error rate       `error_rate` <= 0.01. Added 2026-09-20. There is no natural
                   gap to cut at -- the collection runs from 0.0002 to 0.0139
                   with p50 0.0039 and p95 0.0099 -- so this is a judgement,
                   set to trim the top ~5% tail. It excludes 2,279 of the
                   47,815 isolates that pass the other four criteria (4.8%).

OPTIONAL, OFF BY DEFAULT: `--max-errors-per-read`, which is `avglen_1` x
`error_rate`. That product predicts ASSEMBLY cost far better than either factor
alone -- a 250 bp read at 0.009 carries 2.25 expected errors and every one of
them spawns spurious k-mers, which is why one truth genome took five times
longer to assemble than its peers and SPAdes reported 56.8 million distinct
k-mers for a 4.4 Mb genome. It is deliberately not a default: a noisy sample is
still a valid isolate and is only expensive, so this belongs to assembly-bound
work rather than to the inclusion rule.

NO UPPER DEPTH BOUND, AND WHY THAT IS ONLY SAFE HERE. This collection was
downsampled on ingest: over 54,461 isolates `meandepth` has a median of 88.7
and a hard maximum of 99.7, with nothing above 100. The ceiling is the ingest
pipeline, not the sequencing, so a screen with no upper bound happens to be
safe against THIS table and is not safe against raw SRA data, where depth has a
long right tail. A cohort assembled from outside this collection needs an
explicit upper cap, and any runtime or depth expectation calibrated here is a
lower bound for external data.

Every test over real isolates should run this and report both the criteria and
the number failing each screen, rather than taking whichever isolates happen to
be staged already.

  bin/select_isolates.py --out cohort.tsv                    # screen everything
  bin/select_isolates.py --samples refbias/t8/samples.tsv    # screen a cohort
"""
import argparse, csv, sys

TABLE = ("/n/netscratch/sfortune_lab/Lab/pculviner/notebooks/"
         "260728_sv_exploration/completed_results_all.csv")


def truthy(v):
    return str(v).strip().lower() in ("true", "1", "yes", "t")


def num(v):
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        return 0.0


def screen(row, min_depth, require_complete, max_error=None, max_epr=None):
    """Return the list of reasons this isolate is excluded; empty means keep."""
    bad = []
    if not truthy(row.get("is_paired")):
        bad.append("single_end")
    if num(row.get("meandepth")) < min_depth:
        bad.append("low_depth")
    mixed = bool((row.get("conflict_lineages") or "").strip()) \
        or num(row.get("n_mixed")) > 0 or num(row.get("mix_freq")) > 0
    if mixed:
        bad.append("mixed")
    if require_complete and (row.get("exit_type") or "").strip() != "complete":
        bad.append("incomplete")
    if max_error is not None and num(row.get("error_rate")) > max_error:
        bad.append("high_error_rate")
    if max_epr is not None:
        epr = num(row.get("avglen_1")) * num(row.get("error_rate"))
        if epr > max_epr:
            bad.append("high_errors_per_read")
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default=TABLE)
    ap.add_argument("--samples", help="restrict to these isolates "
                                      "(first column of a TSV)")
    ap.add_argument("--min-depth", type=float, default=60.0)
    ap.add_argument("--no-require-complete", action="store_true")
    ap.add_argument("--max-error-rate", type=float, default=0.01,
                    help="ceiling on error_rate; 0 disables")
    ap.add_argument("--max-errors-per-read", type=float, default=None,
                    help="ceiling on avglen_1 x error_rate. OFF by default. "
                         "Predicts assembly cost, not call quality; use it for "
                         "assembly-bound cohorts.")
    ap.add_argument("--out", help="TSV of passing isolates with their QC")
    a = ap.parse_args()

    keep = None
    if a.samples:
        keep = {l.split("\t")[0].strip() for l in open(a.samples) if l.strip()}
        # a header row is not an isolate; it was being reported as "not kept"
        keep.discard("sample")

    rows, reasons, n = [], {}, 0
    with open(a.table) as fh:
        rd = csv.DictReader(fh)
        idcol = rd.fieldnames[0]
        for r in rd:
            s = (r.get(idcol) or "").strip()
            if keep is not None and s not in keep:
                continue
            n += 1
            bad = screen(r, a.min_depth, not a.no_require_complete,
                         a.max_error_rate if a.max_error_rate > 0 else None,
                         a.max_errors_per_read)
            if bad:
                for b in bad:
                    reasons[b] = reasons.get(b, 0) + 1
                continue
            rows.append(dict(sample=s, lineage=r.get("call") or "",
                             deepest=r.get("deepest_call") or "",
                             meandepth=round(num(r.get("meandepth")), 1),
                             coverage=round(num(r.get("coverage")), 2),
                             mapping_rate=round(num(r.get("mapping_rate")), 4),
                             error_rate=round(num(r.get("error_rate")), 5),
                             read_len=round(num(r.get("avglen_1"))),
                             errors_per_read=round(
                                 num(r.get("avglen_1")) * num(r.get("error_rate")), 2),
                             duplicate_rate=round(num(r.get("duplicate_rate")), 4)))

    print(f"  screened {n} isolates; {len(rows)} pass "
          f"(paired-end, >= {a.min_depth:g}x, not mixed"
          f"{'' if a.no_require_complete else ', complete'}"
          f"{f', error_rate <= {a.max_error_rate:g}' if a.max_error_rate > 0 else ''}"
          f"{f', errors/read <= {a.max_errors_per_read:g}' if a.max_errors_per_read else ''})")
    if n:
        for k, v in sorted(reasons.items(), key=lambda x: -x[1]):
            print(f"    excluded {v:>6d} ({100*v/n:>5.1f}%)  {k}")
    if keep is not None:
        missing = keep - {r["sample"] for r in rows} - set()
        absent = keep - {r["sample"] for r in rows}
        if absent:
            print(f"    of the requested cohort, {len(absent)} not kept: "
                  f"{', '.join(sorted(absent))}")
    if a.out and rows:
        with open(a.out, "w") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader(); w.writerows(rows)
        print(f"  written: {a.out}")
    return 0 if rows else 1


if __name__ == "__main__":
    sys.exit(main())
