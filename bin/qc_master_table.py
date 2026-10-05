#!/usr/bin/env python3
"""Join every per-genome QC measurement into one auditable table.

One row per assembly in the 514-genome selection, so that a panel decision can be
traced to the numbers behind it rather than to a narrative. Sources:

  assembly_summary_refseq   organism, strain, biosample, submitter, date
  assembly_qc_stats.tsv     length, contigs, N content, IUPAC ambiguity codes
  collinearity_qc.tsv       large blocks, order correlation, breakpoints,
                            inverted fraction, reference coverage
  rotation_verify.tsv       achieved dnaA offset and strand (panel genomes only)
  lineages.all.tsv          tb-profiler 6.7.0 barcode call and mixed-lineage flag
  snp_outliers.tsv          within-sublineage robust-z on SNP count
  mash_nearest.tsv          nearest-neighbour mash distance and identity
  assembly_provenance.tsv   reference-guided / short-read-only flag
  panel / excluded          first-round membership and, if excluded, the reason

MEMBERSHIP IS READ FROM THE BUILT PANEL, NOT FROM THE FIRST-ROUND LIST
`--panel-fasta` makes the FASTA that was actually built the authority for
`in_panel`, and the exclusion reason is then resolved from whichever screen
accounts for the genome.  Without it the script falls back to
`panel.rebuild.tsv`, which is the first-round list of 484 and predates the
provenance screen -- so the table then claimed 484 members against a
333-path panel, with both engineered genomes still marked `yes` and a blank
reason.  That mismatch is the reason this option exists.

The resolution order is: an explicit first-round reason, then an engineered
exclusion, then the provenance flag, then `not_in_built_panel` for anything
left, which should be nothing and is reported loudly if it is not.

ONE GENOME IS RETAINED AGAINST ITS OWN FLAG
H37Rv carries a provenance flag because every other genome's structure is
described relative to it, and it is kept anyway as the panel's reference path.
`--retained` lists such cases so they read as a recorded decision rather than
as a screen that failed.

COVERAGE IS UNEVEN AND THE TABLE SAYS SO
SNP counts come from the existing 412-sample complex graph, so the 79 genomes
added in this refresh have no SNP screen yet -- they cannot until the graph is
rebuilt. Rotation and mash cover the 490 panel genomes only, since excluded
genomes were never rotated. Empty cells mean "not measured", never "zero".
"""
import argparse, collections, csv, gzip, os, re, sys


def load(path, key="accession"):
    if not path or not os.path.exists(path):
        return {}
    with open(path) as fh:
        return {r[key]: r for r in csv.DictReader(fh, delimiter="\t")}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qc-dir", required=True)
    ap.add_argument("--ncbi-dir", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--panel-fasta", default="",
                    help="the panel actually built; makes it the authority "
                         "for in_panel")
    ap.add_argument("--provenance", default="",
                    help="assembly_provenance.tsv, for exclusion reasons; "
                         "defaults to <qc-dir>/assembly_provenance.tsv")
    ap.add_argument("--engineered", default="",
                    help="defaults to <qc-dir>/engineered_exclusions.tsv")
    ap.add_argument("--retained", default="",
                    help="genomes kept despite a flag; defaults to "
                         "<qc-dir>/panel_retained_exceptions.tsv")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    Q, N = a.qc_dir, a.ncbi_dir

    panel = {}
    for i, l in enumerate(open(f"{N}/panel.rebuild.tsv")):
        if i:
            f = l.rstrip("\n").split("\t"); panel[f[0]] = f
    excl = {}
    for i, l in enumerate(open(f"{N}/panel.excluded.tsv")):
        if i:
            f = l.rstrip("\n").split("\t"); excl[f[0]] = f[3] if len(f) > 3 else "excluded"

    meta = {}
    for line in open(a.summary):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        acc = f[0].split(".")[0]
        if acc in panel or acc in excl:
            meta[acc] = {"organism": f[7], "strain": f[8].replace("strain=", ""),
                         "biosample": f[2], "date": f[14],
                         "submitter": f[16] if len(f) > 16 else ""}

    # membership from the built panel, when one is given
    built = set()
    if a.panel_fasta:
        op = gzip.open if a.panel_fasta.endswith(".gz") else open
        with op(a.panel_fasta, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    built.add(line[1:].split("#")[0].split()[0])
        if not built:
            sys.exit(f"no sequence names in {a.panel_fasta}")

    prov = load(a.provenance or f"{Q}/assembly_provenance.tsv")
    eng = load(a.engineered or f"{Q}/engineered_exclusions.tsv")
    retained = load(a.retained or f"{Q}/panel_retained_exceptions.tsv")

    stats = load(f"{Q}/assembly_qc_stats.tsv")
    coll = load(f"{Q}/collinearity_qc.tsv")
    linq = load(f"{Q}/lineages.all.tsv")
    snp = load(f"{Q}/snp_outliers.tsv")
    mash = load(f"{Q}/mash_nearest.tsv")
    rot = {}
    rp = f"{Q}/rotation_verify.tsv"
    if os.path.exists(rp):
        for l in open(rp):
            f = l.rstrip("\n").split("\t")
            if len(f) >= 6:
                rot[f[0]] = {"dnaA_offset": f[1], "dnaA_strand": f[2],
                             "dnaA_alen": f[3], "rot_status": f[5]}

    cols = ["accession", "in_panel", "exclusion_reason",
            "organism", "strain", "biosample", "submitter", "date",
            "lineage", "sublineage", "mixed_lineage",
            "length", "contigs", "N_count", "N_pct", "iupac_ambiguous",
            "large_blocks", "order_rho", "breakpoints", "inverted_frac", "ref_covered",
            "dnaA_offset", "dnaA_strand", "rot_status",
            "snps", "snp_group_median", "snp_excess", "snp_fold", "snp_robust_z", "snp_flag",
            "nearest_mash_dist", "nearest_neighbour"]
    def member(acc):
        return ("yes" if acc in built else "no") if built else \
               ("yes" if acc in panel else "no")

    def reason(acc):
        if member(acc) == "yes":
            return ""
        if acc in excl:
            return excl[acc]
        if acc in eng:
            return eng[acc].get("reason", "engineered")
        f = prov.get(acc, {}).get("flag", "")
        if f and f != "ok":
            return f"provenance:{f}"
        return "not_in_built_panel"

    rows = []
    for acc in sorted(set(panel) | set(excl)):
        m = meta.get(acc, {}); s = stats.get(acc, {}); c = coll.get(acc, {})
        L = linq.get(acc, {}); p = snp.get(acc, {}); mn = mash.get(acc, {})
        r = rot.get(acc, {})
        sl = (L.get("strain", "") or "")
        mm = re.match(r"(lineage\d+|La\d+|M\.\w+)", sl)
        sn = int(p["snps"]) if p.get("snps") else None
        gm = int(p["group_median"]) if p.get("group_median") else None
        rows.append({
            "accession": acc,
            "in_panel": member(acc),
            "exclusion_reason": reason(acc),
            "organism": m.get("organism", ""), "strain": m.get("strain", ""),
            "biosample": m.get("biosample", ""), "submitter": m.get("submitter", ""),
            "date": m.get("date", ""),
            "lineage": mm.group(1) if mm else ("no_call" if not sl else sl),
            "sublineage": sl.split(":")[0], "mixed_lineage": L.get("mixed", ""),
            "length": s.get("length", ""), "contigs": s.get("contigs", ""),
            "N_count": s.get("N_count", ""), "N_pct": s.get("N_pct", ""),
            "iupac_ambiguous": s.get("other_ambiguous", ""),
            "large_blocks": c.get("n_large_blocks", ""), "order_rho": c.get("order_rho", ""),
            "breakpoints": c.get("breakpoints", ""), "inverted_frac": c.get("inverted_frac", ""),
            "ref_covered": c.get("ref_covered", ""),
            "dnaA_offset": r.get("dnaA_offset", ""), "dnaA_strand": r.get("dnaA_strand", ""),
            "rot_status": r.get("rot_status", ""),
            "snps": sn if sn is not None else "", "snp_group_median": gm if gm else "",
            "snp_excess": (sn - gm) if (sn is not None and gm) else "",
            "snp_fold": round(sn / gm, 3) if (sn is not None and gm) else "",
            "snp_robust_z": p.get("robust_z", ""), "snp_flag": p.get("flag", ""),
            "nearest_mash_dist": mn.get("nearest_mash_dist", ""),
            "nearest_neighbour": mn.get("nearest_neighbour", ""),
        })
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader(); w.writerows(rows)

    print(f"[qc_master] {len(rows)} rows x {len(cols)} columns -> {a.out}")
    n_in = sum(1 for r in rows if r["in_panel"] == "yes")
    print(f"  in panel: {n_in}   excluded: {len(rows) - n_in}"
          + (f"   (membership read from {a.panel_fasta})" if built else
             "   (membership from panel.rebuild.tsv; pass --panel-fasta to "
             "use the built panel)"))
    if built:
        miss = built - {r["accession"] for r in rows}
        if miss:
            print(f"  WARNING: {len(miss)} paths in the panel have no QC row: "
                  f"{sorted(miss)[:5]}")
        if n_in != len(built):
            print(f"  WARNING: {n_in} marked in_panel against {len(built)} "
                  f"paths in the panel FASTA")
        print("\n  why each excluded genome is out:")
        c = collections.Counter(r["exclusion_reason"] for r in rows
                                if r["in_panel"] == "no")
        for k, v in c.most_common():
            print(f"    {k or '(blank)':<44s} {v:>4d}")
        if c.get("not_in_built_panel"):
            print("    ^ 'not_in_built_panel' means NO screen accounts for it; "
                  "that is a gap, not a decision")
        for acc in retained:
            if acc in built:
                print(f"\n  retained against its own flag: {acc} -- "
                      f"{retained[acc].get('reason','')}")
    print("\n  measurement coverage (blank = not measured, never zero):")
    for c in ("length", "large_blocks", "dnaA_offset", "sublineage", "snps", "nearest_mash_dist"):
        n = sum(1 for r in rows if str(r[c]) != "")
        print(f"    {c:<20s} {n:>4d} / {len(rows)}")


if __name__ == "__main__":
    main()
