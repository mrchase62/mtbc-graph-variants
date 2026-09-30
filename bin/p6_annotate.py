#!/usr/bin/env python3
"""P6: annotate the merged matrix, two-tier.

  on the H37Rv path    H37Rv's own RefSeq gene models. This is where H37Rv earns
                       its place: not as the coordinate system every genotype
                       lives in, but as the annotation layer, which is what
                       WORKFLOW_PROPOSAL section 2c concluded.
  off the H37Rv path   the MATCHED REFERENCE's own RefSeq GFF, at that sample's
                       own coordinate. The off-path sequence IS that reference's
                       sequence, so its GFF annotates it directly with no
                       projection step. These are the records H37Rv annotation
                       structurally cannot reach.

                       The header used to say ~1-2% of records need this route,
                       citing section 11. Measured on the cohort tables it is
                       7.1% of pilot sites and 12.0% at scale100, so the route
                       is not the rounding error that figure implies.

A site is annotated by gene and product where a feature contains it, and left
explicitly unannotated where none does -- intergenic is a finding, not a gap to
be filled by attaching the nearest gene.
"""
import argparse, bisect, collections, csv, gzip, os, sys


def op(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def load_genes(path, want=("gene", "pseudogene", "CDS")):
    """[(start, end, locus_tag, gene, product)] sorted by start, per contig."""
    by_contig = collections.defaultdict(list)
    prod = {}
    for line in op(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in want:
            continue
        attrs = {}
        for kv in f[8].split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                attrs[k] = v
        s, e = int(f[3]), int(f[4])
        lt = attrs.get("locus_tag") or attrs.get("ID", "").replace("gene-", "")
        if f[2] == "CDS":
            # products live on the CDS; genes carry the name
            if lt:
                prod[lt] = attrs.get("product", "")
            continue
        by_contig[f[0]].append((s, e, lt, attrs.get("gene", ""), ""))
    for c in by_contig:
        by_contig[c].sort()
    # attach products collected from CDS rows
    for c, rows in by_contig.items():
        by_contig[c] = [(s, e, lt, g, prod.get(lt, "")) for s, e, lt, g, _ in rows]
    return by_contig


def make_lookup(rows):
    """The gene containing p, preferring the latest-starting one.

    The backward scan used to stop at the first gene that ENDS before p, so a
    longer, earlier gene still containing p was never reached whenever a short
    gene sat inside it -- 9,418 genic H37Rv bases came back intergenic. The
    running maximum of end coordinates bounds the scan correctly: once no
    earlier gene can reach p, stop.
    """
    starts = [r[0] for r in rows]
    reach, m = [], 0
    for r in rows:
        m = max(m, r[1]); reach.append(m)

    def hit(p):
        i = bisect.bisect_right(starts, p) - 1
        while i >= 0 and reach[i] >= p:
            if rows[i][0] <= p <= rows[i][1]:
                return rows[i]
            i -= 1
        return None
    return hit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", default="refbias/p5/matrix.tsv",
                    help="legacy input; --sites is read instead when given")
    ap.add_argument("--sites", default="",
                    help="sites.tsv from p5_matrix.py: the same per-site "
                         "columns as the dense matrix, without the samples")
    ap.add_argument("--refmap", default="refbias/p1/refmap.tsv")
    ap.add_argument("--p4-dir", default="refbias/p4")
    ap.add_argument("--build", required=True)
    ap.add_argument("--h37rv-acc", default="GCF_000195955")
    ap.add_argument("--out", default="refbias/p6/annotated.tsv")
    a = ap.parse_args()

    gff = os.path.join(a.build, "annotation", f"{a.h37rv_acc}.gff.gz")
    if not os.path.exists(gff):
        print(f"no H37Rv GFF at {gff}", file=sys.stderr); return 1
    h37 = load_genes(gff)
    contig = "NC_000962.3"
    if contig not in h37:
        print(f"H37Rv GFF has no {contig}: {list(h37)[:3]}", file=sys.stderr)
        return 1
    hit_h37 = make_lookup(h37[contig])
    print(f"  H37Rv annotation: {len(h37[contig])} genes on {contig}")

    # off-path records need their carrier's reference coordinate, which the
    # matrix does not carry: it is in P4's per-sample output
    offpath_src = {}
    for r in csv.DictReader(open(a.refmap, newline=""), delimiter="\t"):
        p = os.path.join(a.p4_dir, f"{r['sample']}.placed.tsv")
        if not os.path.exists(p):
            continue
        for x in csv.DictReader(open(p, newline=""), delimiter="\t"):
            if x["frame"] == "node" and x["key"] not in offpath_src:
                offpath_src[x["key"]] = (r["reference"], x)
    print(f"  off-path keys with a carrier to annotate from: {len(offpath_src)}")

    # THE TWO KEYS ARE NOT THE SAME STRING. P4 writes `node:<id>:<offset>`,
    # while p5_keys.py writes the COHORT key through mtb_norm.node_key, which
    # appends the normalised allele: `node:<id>:<offset>:<REF>><ALT>`. Looking
    # the matrix key up in P4's index therefore never matched, and every
    # off-path row fell through to offpath_no_carrier -- 7,290 of 60,802 rows
    # at scale100, which is the whole of the second annotation tier. Match on
    # the node prefix, which is what both agree on.
    def carrier_of(key):
        hit = offpath_src.get(key)
        if hit is not None:
            return hit
        return offpath_src.get(":".join(key.split(":")[:3]), (None, None))

    ref_gff = {}

    def annotate_offpath(refid, rec):
        """Annotate from the matched reference's own GFF, at its own coordinate."""
        if refid not in ref_gff:
            p = os.path.join(a.build, "annotation", f"{refid}.gff.gz")
            if os.path.exists(p):
                g = load_genes(p)
                ref_gff[refid] = {c: make_lookup(rows) for c, rows in g.items()}
            else:
                ref_gff[refid] = {}
        lut = ref_gff[refid]
        if not lut:
            return None, "no_gff"
        rp = rec.get("r_pos") or ""
        if not rp.isdigit():
            return None, "no_r_coordinate"
        rp = int(rp)
        # the reference may have several contigs; try each, since the record does
        # not carry its contig name
        for c, hit in lut.items():
            g = hit(rp)
            if g:
                return g, "annotated"
        return None, "intergenic"

    rows = []
    counts = collections.Counter()
    # Only the per-site columns are read, so sites.tsv serves exactly as the
    # dense matrix did, at a fraction of the size.
    with open(a.sites or a.matrix, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        fixed = hdr.index("n_nocall") + 1
        idx = {k: hdr.index(k) for k in
               ("key", "frame", "region", "kind", "h37rv_pos", "node",
                "acc_locus", "canonical_ref", "n_alt", "n_ref", "n_absent",
                "n_nocall")}
        for row in rd:
            frame = row[idx["frame"]]
            gene = locus = product = ""
            src = ""
            if frame == "h37rv":
                pos = int(row[idx["h37rv_pos"]] or 0)
                g = hit_h37(pos) if pos else None
                if g:
                    locus, gene, product = g[2], g[3], g[4]
                    src = "h37rv_gff"
                    counts["h37rv_annotated"] += 1
                else:
                    src = "h37rv_intergenic"
                    counts["h37rv_intergenic"] += 1
            else:
                refid, rec = carrier_of(row[idx["key"]])
                if refid is None:
                    src = "offpath_no_carrier"
                else:
                    g, why = annotate_offpath(refid, rec)
                    src = f"offpath_{why}"
                    if g:
                        locus, gene, product = g[2], g[3], g[4]
                    else:
                        locus = row[idx["acc_locus"]]
                counts[src] += 1
            rows.append([row[idx[k]] for k in
                         ("key", "frame", "region", "kind", "h37rv_pos", "node",
                          "acc_locus", "canonical_ref", "n_alt", "n_ref",
                          "n_absent", "n_nocall")]
                        + [locus, gene, product, src])

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["key", "frame", "region", "kind", "h37rv_pos", "node",
                    "acc_locus", "canonical_ref", "n_alt", "n_ref", "n_absent",
                    "n_nocall", "locus_tag", "gene", "product", "annotation_source"])
        w.writerows(rows)

    print(f"\n  {len(rows)} sites annotated")
    for k, n in counts.most_common():
        print(f"    {k:<26s}{n:>7d}  {100*n/len(rows):>5.1f}%")
    genes = collections.Counter(r[13] for r in rows if r[13])
    print(f"\n  distinct genes carrying a variant: {len(genes)}")
    for g, n in genes.most_common(8):
        print(f"    {g:<12s}{n:>5d} sites")
    print(f"\n  written: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
