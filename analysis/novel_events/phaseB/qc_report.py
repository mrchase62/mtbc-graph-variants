#!/usr/bin/env python3
"""QC report for a QC'd Phase B cohort (run from runroot after 08_qc, 09_kraken
and make_qc_cohort.py). One row per isolate in the Phase B set, with:

  run          study, instrument, read length (ENA, foldback_sim/ena_library.tsv)
  fold-back    chimeras per 1,000 reads (foldback_qc/phaseB_rates.tsv); isolates
               above the cutoff are listed as excluded before QC
  fastp        reads and bases in and out, reads removed by reason, Q30 before and
               after, adapter-trimmed reads, overlap-corrected reads, duplication,
               insert-size peak
  TB-Profiler  lineage, sub-lineage, mixed status, drug-resistance type,
               % reads mapped, median depth
  kraken2      MTBC share of classified reads, largest other taxon
  depth        raw and final (after rasusa)
  status       qc.tsv from qc_summary.py (the user's rules)

Writes <cohort dir>/qc_report.tsv and qc_report.html (self-contained).
"""
import argparse
import csv
import html
import json
import os
import statistics as st

GENOME = 4411532


def rd(path):
    return list(csv.DictReader(open(path), delimiter="\t")) if os.path.exists(path) else []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--cutoff", type=float, default=2.0)
    a = ap.parse_args()
    N = "../analysis/novel_events/phaseB"
    base = f"refbias/{a.cohort}"
    samples = {r["sample"]: r for r in rd(f"{N}/out/samples.tsv")}
    fb = {r["sample"]: float(r["inv_chimera_per_1k"]) for r in rd("../analysis/foldback_qc/phaseB_rates.tsv")}
    ena = {}
    for r in rd("../analysis/foldback_sim/ena_library.tsv"):
        if r["set"] == "phaseB":
            ena[r["isolate"]] = r
    runs = {"mar_" + r["sample"].replace("-", "_"): r for r in rd(f"{N}/runs.tsv")}
    rows = []
    for s in sorted(samples, key=lambda x: (samples[x]["group"], x)):
        sm, e = samples[s], ena.get(s, {})
        rl = round(int(e.get("base_count") or 0) / max(1, int(e.get("read_count") or 1)) / 2) if e else ""
        row = dict(sample=s, group=sm["group"], study=e.get("study_accession", ""),
                   instrument=e.get("instrument_model", ""), read_length=rl,
                   foldback_per_1k=f"{fb.get(s, float('nan')):.2f}")
        q = f"{base}/qc/{s}"
        if fb.get(s, 0) > a.cutoff:
            row.update(status="EXCLUDED", reason=f"fold-back {fb[s]:.1f} per 1,000 > {a.cutoff:g}")
            rows.append(row)
            continue
        if os.path.exists(f"{q}/fastp.json"):
            d = json.load(open(f"{q}/fastp.json"))
            b, f, fr = d["summary"]["before_filtering"], d["summary"]["after_filtering"], d["filtering_result"]
            row.update(raw_depth=f"{b['total_bases'] / GENOME:.1f}",
                       reads_in=b["total_reads"], reads_out=f["total_reads"],
                       reads_removed_pct=f"{100 * (1 - f['total_reads'] / b['total_reads']):.2f}",
                       bases_removed_pct=f"{100 * (1 - f['total_bases'] / b['total_bases']):.2f}",
                       low_quality_reads=fr["low_quality_reads"], too_many_N_reads=fr["too_many_N_reads"],
                       too_short_reads=fr["too_short_reads"],
                       q30_before=f"{100 * b['q30_rate']:.1f}", q30_after=f"{100 * f['q30_rate']:.1f}",
                       gc=f"{100 * f['gc_content']:.1f}",
                       adapter_trimmed_pct=f"{100 * d.get('adapter_cutting', {}).get('adapter_trimmed_reads', 0) / b['total_reads']:.2f}",
                       corrected_reads_pct=f"{100 * fr.get('corrected_reads', 0) / b['total_reads']:.2f}",
                       duplication_pct=f"{100 * d.get('duplication', {}).get('rate', 0):.2f}",
                       insert_peak=d.get("insert_size", {}).get("peak", ""))
        if os.path.exists(f"{q}/tbprofiler.json"):
            t = json.load(open(f"{q}/tbprofiler.json"))
            qc = t.get("qc", {})
            row.update(lineage=t.get("main_lineage", ""), sub_lineage=t.get("sub_lineage", ""),
                       dr_type=t.get("drtype", ""),
                       tbp_pct_mapped=f"{qc.get('percent_reads_mapped', '')}",
                       tbp_median_depth=f"{qc.get('genome_median_depth', '')}")
        qt = rd(f"{q}/qc.tsv")
        if qt:
            x = qt[0]
            row.update(tbprofiler=x["tbprofiler"],
                       mtbc_pct=f"{100 * float(x['mtbc_frac']):.2f}" if x["mtbc_frac"] else "",
                       top_other=x["top_other"],
                       top_other_pct=f"{100 * float(x['top_other_frac']):.2f}" if x["top_other_frac"] else "",
                       final_depth=x["final_depth"], status=x["status"], reason=x["reason"])
        else:
            row.update(status="MISSING", reason="no QC output")
        rows.append(row)
    cols = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(f"{base}/qc_report.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
    write_html(f"{base}/qc_report.html", a, rows, cols)
    print(f"{len(rows)} isolates: " + ", ".join(
        f"{k} {sum(1 for r in rows if r.get('status') == k)}" for k in ("PASS", "FAIL", "EXCLUDED", "MISSING")))


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def write_html(path, a, rows, cols):
    qcd = [r for r in rows if r.get("status") in ("PASS", "FAIL")]
    metrics = [("raw_depth", "Raw depth (x)"), ("reads_removed_pct", "Reads removed by fastp (%)"),
               ("bases_removed_pct", "Bases removed by fastp (%)"), ("q30_before", "Q30 before (%)"),
               ("q30_after", "Q30 after (%)"), ("adapter_trimmed_pct", "Adapter-trimmed reads (%)"),
               ("duplication_pct", "Duplication (%)"), ("insert_peak", "Insert-size peak (bp)"),
               ("mtbc_pct", "MTBC, % of classified (kraken2)"), ("top_other_pct", "Largest other taxon (%)"),
               ("final_depth", "Final depth after downsampling (x)"), ("foldback_per_1k", "Fold-back per 1,000 reads")]
    srow = []
    for k, lab in metrics:
        v = [num(r.get(k)) for r in qcd]
        v = [x for x in v if x is not None]
        if v:
            srow.append(f"<tr><td>{lab}</td><td class=n>{min(v):.2f}</td><td class=n>{st.median(v):.2f}</td>"
                        f"<td class=n>{max(v):.2f}</td></tr>")
    groups = {}
    for r in rows:
        g = groups.setdefault(r["group"], {"PASS": 0, "FAIL": 0, "EXCLUDED": 0, "MISSING": 0})
        g[r.get("status", "MISSING")] += 1
    grow = "".join(f"<tr><td>{g}</td><td class=n>{c['PASS']}</td><td class=n>{c['FAIL']}</td>"
                   f"<td class=n>{c['EXCLUDED']}</td></tr>" for g, c in groups.items())
    show = ["sample", "group", "status", "reason", "study", "instrument", "read_length", "foldback_per_1k",
            "raw_depth", "reads_removed_pct", "bases_removed_pct", "q30_before", "q30_after",
            "adapter_trimmed_pct", "duplication_pct", "insert_peak", "lineage", "sub_lineage", "tbprofiler",
            "dr_type", "mtbc_pct", "top_other", "top_other_pct", "final_depth"]
    head = "".join(f"<th data-k='{i}'>{html.escape(c.replace('_', ' '))}</th>" for i, c in enumerate(show))
    body = ""
    for r in rows:
        cls = {"PASS": "ok", "FAIL": "bad", "EXCLUDED": "ex"}.get(r.get("status"), "bad")
        body += f"<tr class='{cls}'>" + "".join(
            f"<td{' class=n' if num(r.get(c)) is not None else ''}>{html.escape(str(r.get(c, '')))}</td>" for c in show) + "</tr>"
    fails = [r for r in rows if r.get("status") == "FAIL"]
    fail_text = ("<p>No isolate failed the QC rules.</p>" if not fails else
                 "<ul>" + "".join(f"<li><b>{html.escape(r['sample'])}</b>: {html.escape(r['reason'])}</li>" for r in fails) + "</ul>")
    page = TEMPLATE.format(cohort=html.escape(a.cohort), cutoff=a.cutoff, grow=grow, srow="".join(srow),
                           head=head, body=body, fail_text=fail_text, n=len(rows),
                           npass=sum(1 for r in rows if r.get("status") == "PASS"))
    open(path, "w").write(page)


TEMPLATE = """<title>Phase B Read QC</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;600&family=IBM+Plex+Mono&family=Source+Serif+4:opsz,wght@8..60,600&display=swap">
<style>
/* Layout: summary first (counts, ranges, failures), then the full per-isolate table in its own scroll box. */
:root {{ --bg:#f7f8f6; --panel:#fff; --fg:#1d2421; --muted:#5b6662; --line:#d3dad6; --ok:#127a6b; --bad:#b3261e; --ex:#8a6d1f; --okbg:#e3f2ee; --badbg:#fae3e1; --exbg:#f5ecd6;
  --sans:"IBM Plex Sans",Arial,sans-serif; --mono:"IBM Plex Mono",Menlo,monospace; --serif:"Source Serif 4",Georgia,serif; }}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{ --bg:#141917; --panel:#1b2220; --fg:#e3e9e6; --muted:#9aa7a2; --line:#34403c; --ok:#4cc3ae; --bad:#f08a80; --ex:#e0c26a; --okbg:#183a34; --badbg:#3d1d1a; --exbg:#3a3118; color-scheme:dark; }} }}
:root[data-theme="dark"] {{ --bg:#141917; --panel:#1b2220; --fg:#e3e9e6; --muted:#9aa7a2; --line:#34403c; --ok:#4cc3ae; --bad:#f08a80; --ex:#e0c26a; --okbg:#183a34; --badbg:#3d1d1a; --exbg:#3a3118; color-scheme:dark; }}
body {{ background:var(--bg); color:var(--fg); font-family:var(--sans); font-size:15px; line-height:1.55; }}
.wrap {{ max-width:1200px; margin:0 auto; padding-inline:18px; padding-block:32px 56px; }}
h1 {{ font-family:var(--serif); font-size:2rem; margin:0 0 .3em; }}
h2 {{ font-family:var(--serif); font-size:1.3rem; margin:1.8em 0 .5em; }}
p, li {{ max-width:75ch; }}
.muted {{ color:var(--muted); }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(300px,1fr)); gap:20px; align-items:start; }}
.grid > div {{ min-width:0; }}
table {{ border-collapse:collapse; font-variant-numeric:tabular-nums; font-size:.88rem; width:100%; }}
th, td {{ padding:5px 8px; border-bottom:1px solid var(--line); text-align:left; white-space:nowrap; }}
th {{ color:var(--muted); font-weight:600; font-size:.78rem; position:sticky; top:0; background:var(--panel); cursor:pointer; }}
td.n {{ text-align:right; font-family:var(--mono); font-size:.82rem; }}
.box {{ background:var(--panel); border:1px solid var(--line); border-radius:6px; padding:12px; overflow:auto; }}
.big {{ max-height:70vh; }}
tr.ok td:nth-child(3) {{ color:var(--ok); font-weight:600; }} tr.bad td:nth-child(3) {{ color:var(--bad); font-weight:600; }} tr.ex td:nth-child(3) {{ color:var(--ex); font-weight:600; }}
tr.bad {{ background:var(--badbg); }} tr.ex {{ background:var(--exbg); }}
code {{ font-family:var(--mono); font-size:.85em; }}
</style>
<div class="wrap">
<h1>Phase B read QC</h1>
<p class="muted">Cohort {cohort}: {npass} of {n} isolates pass. Rules (the lab standard): fastp; TB-Profiler with mixed samples removed (two top-level lineages each at 10% or more); kraken2 on raw reads, MTBC at least 85% of classified reads and no other species or genus above 5%; rasusa to 100x, failing below 60x. Isolates with more than {cutoff:g} fold-back chimeras per 1,000 reads were excluded before QC.</p>
<div class="grid">
<div><h2>Outcome by group</h2><div class="box"><table><tr><th>Group</th><th>Pass</th><th>Fail</th><th>Excluded (fold-back)</th></tr>{grow}</table></div></div>
<div><h2>Failures</h2>{fail_text}</div>
</div>
<h2>Ranges over the QC'd isolates</h2>
<div class="box"><table><tr><th>Metric</th><th>Min</th><th>Median</th><th>Max</th></tr>{srow}</table></div>
<p class="muted">fastp settings: <code>--detect_adapter_for_pe --correction --qualified_quality_phred 20 --length_required 36</code>. Raw depth is the input to fastp; final depth is after downsampling. Fold-back rates were measured on the earlier alignments of the raw reads.</p>
<h2>All isolates</h2>
<p class="muted">Click a column header to sort.</p>
<div class="box big"><table id="t"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>
</div>
<script>
(function(){{var t=document.getElementById('t'),dir={{}};
t.querySelectorAll('th').forEach(function(th){{th.addEventListener('click',function(){{var k=+th.dataset.k;dir[k]=!dir[k];
var rows=Array.from(t.tBodies[0].rows);rows.sort(function(a,b){{var x=a.cells[k].textContent,y=b.cells[k].textContent,nx=parseFloat(x),ny=parseFloat(y);
var c=(!isNaN(nx)&&!isNaN(ny))?nx-ny:x.localeCompare(y);return dir[k]?c:-c;}});rows.forEach(function(r){{t.tBodies[0].appendChild(r);}});}});}});}})();
</script>
"""

if __name__ == "__main__":
    main()
