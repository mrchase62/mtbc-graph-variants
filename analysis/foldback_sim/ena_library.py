#!/usr/bin/env python3
"""ENA run metadata (library preparation, instrument) for the isolates whose
fold-back chimera rate was measured: the 500 production samples
(production_chimera.tsv, by BioSample) and the 62 Phase B isolates
(phaseB/runs.tsv, by run accession). Metadata only, from the ENA portal API
search endpoint (result=read_run); no reads are downloaded. One row per run;
a BioSample with several runs gives several rows."""
import argparse
import csv
import io
import time
import urllib.parse
import urllib.request

API = "https://www.ebi.ac.uk/ena/portal/api/search"
FIELDS = ("run_accession,sample_accession,secondary_sample_accession,study_accession,"
          "instrument_platform,instrument_model,library_layout,library_strategy,"
          "library_source,library_selection,library_name,library_construction_protocol,"
          "nominal_length,read_count,base_count,center_name,first_public")


def fetch(key, ids):
    rows = []
    for i in range(0, len(ids), 40):
        q = " OR ".join(f'{key}="{x}"' for x in ids[i:i + 40])
        data = urllib.parse.urlencode(dict(result="read_run", query=q, fields=FIELDS,
                                           format="tsv", limit=0)).encode()
        for attempt in range(3):
            try:
                with urllib.request.urlopen(urllib.request.Request(API, data=data), timeout=120) as r:
                    rows += list(csv.DictReader(io.StringIO(r.read().decode()), delimiter="\t"))
                break
            except OSError:
                if attempt == 2:
                    raise
                time.sleep(5)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--production", required=True, help="production_chimera.tsv")
    ap.add_argument("--phaseb", required=True, help="phaseB/runs.tsv")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    prod = [r["sample"] for r in csv.DictReader(open(a.production), delimiter="\t")]
    pb = {r["run"]: "mar_" + r["sample"].replace("-", "_")
          for r in csv.DictReader(open(a.phaseb), delimiter="\t")}
    rows = [dict(r, set="production", isolate=r["sample_accession"]) for r in fetch("sample_accession", prod)]
    rows += [dict(r, set="phaseB", isolate=pb.get(r["run_accession"], "")) for r in fetch("run_accession", list(pb))]
    with open(a.out, "w") as fo:
        w = csv.DictWriter(fo, fieldnames=["set", "isolate"] + FIELDS.split(","), delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    got = {r["isolate"] for r in rows}
    print(f"{len(rows)} runs; production {len(set(prod) & got)}/{len(prod)}, "
          f"phaseB {len(set(pb.values()) & got)}/{len(pb)}")


if __name__ == "__main__":
    main()
