#!/usr/bin/env python3
"""Apply the agreed exclusions to the 514 selected assemblies and emit the panel.

Decisions recorded here so the panel is reproducible from the manifest plus QC:

  EXCLUDE scrambled            collinearity: reference coverage < 0.97 together
                               with large-block order violations. Seven genomes,
                               all one 2023-08-01 submission batch.
  EXCLUDE mixed lineage        tb-profiler reporting two incompatible lineages.
  EXCLUDE engineered           deliberate genetic modification (mc2 6030,
                               BCG Danish delta-sapM).
  EXCLUDE BCG and H37Ra        attenuated by serial laboratory passage. BCG is
                               matched on the tb-profiler BARCODE CALL
                               (La1.2.BCG), not on the submitter-supplied
                               organism/strain name. Name matching missed BCG
                               Danish 1331 (GCF_005156105), whose RefSeq record
                               says only "Mycobacterium tuberculosis variant
                               bovis | Danish 1331" -- nothing says BCG. The same
                               class of error dropped the M. canettii outgroup
                               earlier, because RefSeq spells it "canetti".
                               Submitter metadata is not a reliable filter;
                               the barcode is computed from the sequence.
  EXCLUDE CDC1551              per instruction.
  EXCLUDE BioSample duplicate  same isolate assembled twice; keep the first.

  KEEP all H37Rv assemblies    they differ structurally and the differences are
                               the biology: three of five lack the esxN.2 /
                               esxJ.3 paralogs that Chitale 2022 added, two
                               carry them.
  KEEP all Erdman              used in the lab.
  KEEP the lineage-2.2.1       eight genomes sharing a real ~20% inverted
  inversion clade              segment at 0.9997 reference coverage. The 2025
                               block-count QC would have failed these; they are
                               structural variation, not misassembly.

Reference stays NC_000962.3 (GCF_000195955) so the snpEff database, RD BED and
IS6110 GFF continue to apply, even though that assembly is the one MISSING the
esx paralogs.
"""
import argparse, csv, re, sys

SCRAMBLED_COV = 0.97
INVERSION_CLADE = set("""GCF_000193185 GCF_003265005 GCF_014899425 GCF_014900195
GCF_016917775 GCF_026167585 GCF_027912555 GCF_014899745""".split())
ENGINEERED = {"GCF_045345655": "mc2 6030", "GCF_005155785": "BCG Danish dsapM"}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--collinearity", required=True)
    ap.add_argument("--summary", required=True, help="RefSeq assembly_summary")
    ap.add_argument("--lineages", help="legacy tb-profiler table (semicolon calls)")
    ap.add_argument("--manual-exclusions", help="TSV of accession<TAB>reason for "
                    "genomes removed on QC grounds that no automated rule catches "
                    "-- currently SNP excess within sublineage, which needs a graph "
                    "to measure and so cannot be a rule applied before one exists")
    ap.add_argument("--barcode", help="lineages.all.tsv from tbprofiler_collect.py: "
                    "authoritative source for mixed-lineage and BCG calls, because "
                    "it is derived from the sequence rather than from metadata")
    ap.add_argument("--out", required=True)
    ap.add_argument("--excluded", required=True)
    a = ap.parse_args()

    man = {}
    for line in open(a.manifest):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3:
            man[f[0]] = f[2]

    col = {r["accession"]: r for r in csv.DictReader(open(a.collinearity), delimiter="\t")}

    manual = {}
    if a.manual_exclusions:
        for r in csv.DictReader(open(a.manual_exclusions), delimiter="\t"):
            manual[r["accession"]] = r["reason"]

    biosample, strain = {}, {}
    for line in open(a.summary):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        acc = f[0].split(".")[0]
        if acc in man:
            biosample[acc] = f[2]
            strain[acc] = f[8].replace("strain=", "") if len(f) > 8 else ""

    mixed, bcg_barcode, call = set(), set(), {}
    if a.barcode:
        for r in csv.DictReader(open(a.barcode), delimiter="\t"):
            acc = r["accession"]
            call[acc] = r.get("strain", "") or ""
            if r.get("mixed", "") == "True":
                mixed.add(acc)
            if "BCG" in call[acc]:
                bcg_barcode.add(acc)
    if a.lineages:
        with open(a.lineages) as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > 2 and (";" in f[1] or ";" in f[2]):
                    mixed.add(f[0].split(".")[0])

    seen_bs, keep, drop = {}, [], []
    for acc in sorted(man):
        org = man[acc]; st = strain.get(acc, "")
        c = col.get(acc)
        why = None
        if acc in manual:
            why = "manual_qc:" + manual[acc].split(":")[0]
        elif c and float(c["ref_covered"]) < SCRAMBLED_COV and int(c["breakpoints"]) > 0:
            why = "scrambled"
        elif acc in mixed:
            why = "mixed_lineage"
        elif acc in ENGINEERED:
            why = "engineered:" + ENGINEERED[acc]
        elif acc in bcg_barcode:
            why = f"BCG_by_barcode:{call[acc]}"
        elif re.search(r"BCG", org + st):
            why = "BCG_by_name"
        elif re.search(r"H37Ra", org + st):
            why = "H37Ra_attenuated"
        elif re.search(r"CDC1551", org + st):
            why = "CDC1551_dropped"
        else:
            bs = biosample.get(acc, "")
            if bs and bs in seen_bs:
                why = f"duplicate_biosample_of:{seen_bs[bs]}"
            elif bs:
                seen_bs[bs] = acc
        (drop if why else keep).append((acc, org, st, why or ""))

    with open(a.out, "w") as fh:
        fh.write("accession\torganism\tstrain\n")
        for acc, org, st, _ in keep:
            fh.write(f"{acc}\t{org}\t{st}\n")
    with open(a.excluded, "w") as fh:
        fh.write("accession\torganism\tstrain\treason\n")
        for acc, org, st, why in drop:
            fh.write(f"{acc}\t{org}\t{st}\t{why}\n")

    import collections
    print(f"[build_panel] {len(man)} selected -> {len(keep)} kept, {len(drop)} excluded")
    for r, n in collections.Counter(w.split(":")[0] for _, _, _, w in drop).most_common():
        print(f"    {r:<26s} {n}")
    inv_kept = sum(1 for acc, _, _, _ in keep if acc in INVERSION_CLADE)
    print(f"  inversion clade retained: {inv_kept}/8")
    if bcg_barcode:
        nb = [x for x in bcg_barcode if not any(x == d[0] for d in drop)]
        print(f"  BCG by barcode: {len(bcg_barcode)} found, "
              f"{len(bcg_barcode)-len(nb)} excluded")
    h = [acc for acc, o, s, _ in keep if "H37Rv" in o + s]
    print(f"  H37Rv assemblies retained: {len(h)}  {' '.join(h)}")
    e = [acc for acc, o, s, _ in keep if "Erdman" in o + s]
    print(f"  Erdman retained: {len(e)}  {' '.join(e)}")


if __name__ == "__main__":
    main()
