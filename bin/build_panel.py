#!/usr/bin/env python3
"""Apply the recorded exclusions to the selected assemblies and emit the panel list.

Every decision between "selected from RefSeq" and "a path in the graph" is
made here, from recorded inputs, so the panel is reproducible from the
manifest plus QC tables. CX333 (514 selected -> 333) is reproduced by:

    build_panel.py --manifest data/ncbi/panel_manifest.tsv \\
        --collinearity data/qc/collinearity_qc.tsv \\
        --summary data/ncbi/assembly_summary_refseq.txt \\
        --barcode data/qc/lineages.all.tsv \\
        --manual-exclusions data/qc/manual_exclusions.tsv \\
        --exclusions data/qc/engineered_exclusions.tsv \\
        --provenance data/qc/assembly_provenance.tsv \\
        --retained data/qc/panel_retained_exceptions.tsv \\
        --out panel.tsv --excluded excluded.tsv

The output is a TSV list (accession, organism, strain), the input of
make_pansn_fasta.sh --panel. It is never a FASTA: the script refuses an --out
or --excluded that looks like a FASTA or is one of its own inputs (the
runbook once passed the panel FASTA as --out, which would have overwritten it
with this list; audit PGB-3).

Decisions recorded here:

  EXCLUDE recorded decisions   --manual-exclusions and --exclusions: TSVs of
                               accession<TAB>reason, one row per genome
                               removed by review (SNP excess within
                               sublineage; the two attB-vector genomes the
                               foreign screen found).
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
  EXCLUDE provenance           assembly_provenance_screen.py flag other than
                               "ok" (REFERENCE_STRUCTURE, SHORT_READ_SUSPECT).
                               Every genome that reaches this rule must have a
                               row: a missing row is "not screened", never
                               "passed".
  REVIEW --review tables       any other screen table with an `accession` and a
                               `flag` column (e.g. snp_outlier_screen.py,
                               assembly_qc_stats.py's base_run_flag). A
                               non-empty flag on a genome that would otherwise
                               be kept stops the build until the genome is
                               given a recorded decision: an --exclusions row,
                               or a --retained row.

  KEEP --retained              genomes kept against a flag, with the reason:
                               H37Rv (GCF_000195955) is flagged by the
                               provenance screen and kept as the reference path.
  KEEP all H37Rv assemblies    they differ structurally and the differences are
                               the biology: three of five lack the esxN.2 /
                               esxJ.3 paralogs that Chitale 2022 added, two
                               carry them.
  KEEP all Erdman              used in the lab.
  KEEP the lineage-2.2.1       eight genomes sharing a real ~20% inverted
  inversion clade              segment at 0.9997 reference coverage. The 2025
                               block-count QC would have failed these; they are
                               structural variation, not misassembly.

A retained exception overrides the provenance flag and the --review flags only;
it never overrides a recorded exclusion or the first-round rules.

Reference stays NC_000962.3 (GCF_000195955) so the snpEff database, RD BED and
IS6110 GFF continue to apply, even though that assembly is the one MISSING the
esx paralogs.
"""
import argparse, collections, csv, hashlib, os, re, sys

SCRAMBLED_COV = 0.97
INVERSION_CLADE = set("""GCF_000193185 GCF_003265005 GCF_014899425 GCF_014900195
GCF_016917775 GCF_026167585 GCF_027912555 GCF_014899745""".split())
ENGINEERED = {"GCF_045345655": "mc2 6030", "GCF_005155785": "BCG Danish dsapM"}
FASTA_LIKE = (".fa", ".fasta", ".fna", ".gz", ".fa.gz", ".fasta.gz", ".bgz")


def read_manifest(path):
    """accession -> organism. Headerless panel_manifest.tsv / selected.*.tsv
    (accession, versioned, organism, ftp), or a table with a header naming an
    `organism` column. A header row is never a genome."""
    with open(path) as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    if not lines:
        sys.exit(f"empty manifest {path}")
    first = lines[0]
    if first[0] in ("accession", "#accession") or first[0].startswith("#"):
        if "organism" not in first:
            sys.exit(f"{path}: has a header but no 'organism' column; pass the "
                     "selection manifest (data/ncbi/panel_manifest.tsv)")
        oc = first.index("organism")
        return {f[0]: f[oc] for f in lines[1:] if len(f) > oc}
    for f in lines:
        if not re.match(r"GC[AF]_\d+", f[0]):
            sys.exit(f"{path}: row '{f[0]}' is not an assembly accession")
    return {f[0]: f[2] for f in lines if len(f) >= 3}


def read_decisions(path):
    """accession -> reason, from an accession<TAB>reason[<TAB>...] TSV"""
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if not r.get("accession") or not r.get("reason"):
            sys.exit(f"{path}: every row needs an accession and a reason")
        out[r["accession"]] = r["reason"]
    return out


def check_outputs(outs, inputs):
    """refuse to write a FASTA-named list, or over any input"""
    real_in = {os.path.realpath(p) for p in inputs if p}
    seen = set()
    for p in outs:
        rp = os.path.realpath(p)
        if p.lower().endswith(FASTA_LIKE):
            sys.exit(f"refusing --out/--excluded {p}: this script writes a TSV "
                     "list, not a FASTA (make_pansn_fasta.sh writes the FASTA)")
        if rp in real_in:
            sys.exit(f"refusing to overwrite an input: {p}")
        if rp in seen:
            sys.exit(f"--out and --excluded are the same file: {p}")
        seen.add(rp)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True,
                    help="selection manifest from refresh_assemblies.sh "
                         "(data/ncbi/panel_manifest.tsv)")
    ap.add_argument("--collinearity", required=True)
    ap.add_argument("--summary", required=True, help="RefSeq assembly_summary")
    ap.add_argument("--barcode", required=True,
                    help="lineages.all.tsv from tbprofiler_collect.py: "
                         "authoritative source for mixed-lineage and BCG calls, "
                         "because it is derived from the sequence rather than "
                         "from metadata")
    ap.add_argument("--lineages", help="legacy tb-profiler table (semicolon calls)")
    ap.add_argument("--manual-exclusions", help="TSV of accession<TAB>reason for "
                    "genomes removed on QC grounds by review -- currently SNP "
                    "excess within sublineage; written as manual_qc:<reason>")
    ap.add_argument("--exclusions", action="append", default=[],
                    help="further TSV(s) of accession<TAB>reason, e.g. "
                         "engineered_exclusions.tsv; the reason is written as given")
    ap.add_argument("--provenance",
                    help="assembly_provenance.tsv; `flag` other than ok excludes. "
                         "Required unless --no-provenance")
    ap.add_argument("--no-provenance", action="store_true",
                    help="first round only (reproduces the 484 of "
                         "panel.rebuild.tsv); NOT a panel")
    ap.add_argument("--retained", help="TSV of accession<TAB>reason: kept "
                    "against a provenance or --review flag")
    ap.add_argument("--review", action="append", default=[],
                    help="screen table(s) with accession and flag columns; a "
                         "flagged genome needs a recorded decision")
    ap.add_argument("--out", required=True, help="panel list (TSV)")
    ap.add_argument("--excluded", required=True, help="excluded list (TSV)")
    a = ap.parse_args()
    if not a.provenance and not a.no_provenance:
        ap.error("--provenance is required (or --no-provenance for the first "
                 "round only)")

    inputs = [a.manifest, a.collinearity, a.summary, a.barcode, a.lineages,
              a.manual_exclusions, a.provenance, a.retained] + a.exclusions + a.review
    check_outputs([a.out, a.excluded], inputs)

    man = read_manifest(a.manifest)
    col = {r["accession"]: r for r in csv.DictReader(open(a.collinearity), delimiter="\t")}

    manual = read_decisions(a.manual_exclusions) if a.manual_exclusions else {}
    recorded = {}
    for p in a.exclusions:
        recorded.update(read_decisions(p))
    retained = read_decisions(a.retained) if a.retained else {}
    prov = {}
    if a.provenance:
        prov = {r["accession"]: r for r in
                csv.DictReader(open(a.provenance), delimiter="\t")}
        if prov and "flag" not in next(iter(prov.values())):
            sys.exit(f"{a.provenance}: no `flag` column")
    review = collections.defaultdict(list)
    for p in a.review:
        for r in csv.DictReader(open(p), delimiter="\t"):
            if "flag" not in r:
                sys.exit(f"{p}: no `flag` column")
            if (r["flag"] or "").strip():
                review[r["accession"]].append(
                    f"{os.path.basename(p)}:{r['flag'].strip()}")

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
    unscreened, unresolved = [], []
    for acc in sorted(man):
        org = man[acc]; st = strain.get(acc, "")
        c = col.get(acc)
        why = src = None
        if acc in manual:
            why, src = "manual_qc:" + manual[acc].split(":")[0], a.manual_exclusions
        elif acc in recorded:
            why, src = recorded[acc], "exclusions"
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
        if why is None and a.provenance:
            p = prov.get(acc)
            if p is None:
                unscreened.append(acc)
            elif p["flag"] not in ("", "ok") and acc not in retained:
                why, src = f"provenance:{p['flag']}", a.provenance
        if why is None and review.get(acc) and acc not in retained:
            unresolved.append((acc, review[acc]))
        (drop if why else keep).append((acc, org, st, why or "", src or ("rule" if why else "")))

    if unscreened:
        sys.exit(f"{len(unscreened)} genome(s) reach the provenance rule with no "
                 f"row in {a.provenance} (not screened is not passed): "
                 f"{' '.join(unscreened[:10])}")
    if unresolved:
        print(f"[build_panel] {len(unresolved)} genome(s) carry a --review flag "
              "and no recorded decision:", file=sys.stderr)
        for acc, fl in unresolved:
            print(f"    {acc}  {'; '.join(fl)}", file=sys.stderr)
        sys.exit("add each to an --exclusions or --retained file, with the "
                 "reason; nothing written")

    with open(a.out, "w") as fh:
        fh.write("accession\torganism\tstrain\n")
        for acc, org, st, _, _ in keep:
            fh.write(f"{acc}\t{org}\t{st}\n")
    with open(a.excluded, "w") as fh:
        fh.write("accession\torganism\tstrain\treason\tdecided_by\n")
        for acc, org, st, why, src in drop:
            fh.write(f"{acc}\t{org}\t{st}\t{why}\t{src}\n")
    with open(a.out + ".inputs.tsv", "w") as fh:
        fh.write("role\tpath\tsha256\n")
        for role, p in [("manifest", a.manifest), ("collinearity", a.collinearity),
                        ("summary", a.summary), ("barcode", a.barcode),
                        ("lineages", a.lineages),
                        ("manual_exclusions", a.manual_exclusions),
                        ("provenance", a.provenance), ("retained", a.retained)] + \
                [("exclusions", p) for p in a.exclusions] + \
                [("review", p) for p in a.review]:
            if p:
                fh.write(f"{role}\t{os.path.abspath(p)}\t{sha256(p)}\n")

    print(f"[build_panel] {len(man)} selected -> {len(keep)} kept, {len(drop)} excluded")
    for r, n in collections.Counter(w.split(":")[0] for _, _, _, w, _ in drop).most_common():
        print(f"    {r:<26s} {n}")
    inv_kept = sum(1 for acc, *_ in keep if acc in INVERSION_CLADE)
    print(f"  inversion clade retained: {inv_kept}/8")
    nb = [x for x in bcg_barcode if not any(x == d[0] for d in drop)]
    print(f"  BCG by barcode: {len(bcg_barcode)} found, "
          f"{len(bcg_barcode)-len(nb)} excluded")
    kept_ids = {k[0] for k in keep}
    for acc in sorted(retained):
        if acc in kept_ids:
            print(f"  retained against a flag: {acc} -- {retained[acc]}")
    h = [acc for acc, o, s, _, _ in keep if "H37Rv" in o + s]
    print(f"  H37Rv assemblies retained: {len(h)}  {' '.join(h)}")
    e = [acc for acc, o, s, _, _ in keep if "Erdman" in o + s]
    print(f"  Erdman retained: {len(e)}  {' '.join(e)}")
    print(f"  -> {a.out}  (inputs and sha256: {a.out}.inputs.tsv)")


if __name__ == "__main__":
    main()
