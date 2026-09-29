#!/usr/bin/env python3
"""Write IS6110 insertion sites as VCF, per is6110/docs/P1I_VCF_DESIGN.md.

FOUR DECISIONS, TAKEN AND NOT RE-OPENED HERE
  filtered but visible    a one-sided call is written with FILTER=LowSupport
                          rather than dropped, so a tool reading only PASS is
                          unaffected and the rows stop being invisible
  measured SVLEN          the element length is taken from what stage 1 actually
                          found in that sample's own matched reference, not from
                          the canonical 1355
  matched ref authoritative   the per-sample VCF is written in the frame the
                          evidence was observed in. The H37Rv file is derived
  omit GT=0               a copy the matched reference already carries is a
                          reference allele in that frame and gets no record

WHY MOST SITES PRODUCE NO RECORD IN THE PER-SAMPLE FILE
132 of 179 two-sided sites sit at a removed-interval join, meaning the matched
reference carries that copy too -- `site_class = ref_shared`. Against that reference they are REF, not insertions, and
writing them as ALT would assert 132 insertions that are not there. They are not
lost: they are facts about the isolate that belong in the cohort frame, where
P5's four states already distinguish REF from ABSENT. This script emits the
cohort key table for that; it does not modify P5.

WHY <INS:ME:IS6110> AND NOT A BREAKEND PAIR
dysgu already writes <INS> into refbias/p2/<sample>.dysgu.vcf and delly writes
<DEL>/<DUP>/<INV>, so symbolic ALT is the house style and these records sit
beside those. A breakend is for a junction whose far side is unknown or is a
rearrangement; here the far side is a known 1355 bp element. BND would also need
two records and a mate id per site, giving the cohort merge two keys per locus.

SVLEN FOR A NOVEL INSERTION IS NOT DIRECTLY MEASURABLE, and the field says so.
Junction evidence gives the position, not the length of the element inserted
there. What IS measured is the length of that reference's own copies, from the
stage 1 GFF, so the median of those is used and IS6110_SVLEN_SRC records that it
is a per-reference median rather than a per-site measurement. A site at a join
does have a measured length, but those get no ALT record, so in practice the
median is what is written.
"""
import argparse, collections, csv, datetime, gzip, os, statistics, sys


def read_fasta_one(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                break
            name = line[1:].split()[0]
        else:
            buf.append(line.strip())
    return name, "".join(buf)


def load_excisions(path):
    """Both bases adjacent to each excised span; a junction sits AT a seam, so a
    stack reports orig_start-1 or orig_end+1, never orig_start."""
    left, right = {}, {}
    if os.path.exists(path):
        for r in csv.DictReader(open(path), delimiter="\t"):
            s0, e0 = int(r["orig_start"]), int(r["orig_end"])
            left[s0 - 1] = (s0, e0)
            right[e0 + 1] = (s0, e0)
    return left, right


def gff_lengths(path):
    out = []
    if os.path.exists(path):
        for line in open(path):
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) > 4:
                out.append(int(f[4]) - int(f[3]) + 1)
    return out


HDR_COMMON = [
    '##ALT=<ID=INS:ME:IS6110,Description="IS6110 mobile element insertion">',
    '##FILTER=<ID=LowSupport,Description="IS6110_EVIDENCE=one_sided: only one flank of the junction was seen. Promoted rather than dropped; contradicted by an assembly 8.8% of the time against 1.7% for two_sided">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the inserted element">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
    '##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info: name,start,end,polarity">',
    '##INFO=<ID=IS6110_EVIDENCE,Number=1,Type=String,Description="Quality of support, an ordering. two_sided: both element termini seen at the site. one_sided: one flank, promoted after collapsing flank pairs">',
    '##INFO=<ID=IS6110_CLASS,Number=1,Type=String,Description="Relationship to the reference this record is written against, NOT a quality ordering. ref_lacking: the reference does not carry this insertion. ref_shared: it does, and the build removed the interval. In the per-sample file this is always ref_lacking, because a ref_shared copy is a reference allele in that frame and gets no record; both classes appear in the cohort key table">',
    '##INFO=<ID=IS6110_CHROM_SIDE,Number=1,Type=String,Description="Outcome of the independent chromosome-side scan at this site. Corroboration, reported separately from IS6110_EVIDENCE rather than folded into it">',
    '##INFO=<ID=IS6110_READS,Number=1,Type=Integer,Description="Uniquely-placed element-side reads supporting the junction">',
    '##INFO=<ID=IS6110_SPAN,Number=1,Type=Integer,Description="Distance in bp between the two flank junction coordinates">',
    '##INFO=<ID=IS6110_SVLEN_SRC,Number=1,Type=String,Description="How SVLEN was obtained: ref_median (median measured element length in this sample matched reference) or measured_locus">',
    '##INFO=<ID=IS6110_H37RV,Number=1,Type=Integer,Description="H37Rv coordinate from flank placement, both flanks required to agree">',
    '##INFO=<ID=IS6110_H37RV_STATE,Number=1,Type=String,Description="empty: H37Rv carries no copy at this locus. occupied: it does. unplaced: the two flanks did not agree or did not place uniquely">',
    '##INFO=<ID=IS6110_NODE,Number=1,Type=String,Description="Pangenome graph node key, always recorded. It is the only identity a site has when IS6110_H37RV_STATE is unplaced, and is build-scoped: see IS6110_BUILD">',
    '##INFO=<ID=IS6110_BUILD,Number=1,Type=String,Description="Build id; graph node ids are build-scoped">',
    '##INFO=<ID=IS6110_ISM,Number=1,Type=Integer,Description="1 if an ISMapper region lies within 50 bp of the H37Rv coordinate">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--reconcile", default="is6110/results/p1i_reconcile.tsv")
    ap.add_argument("--flank", default="is6110/results/p1i_sites_flank_all.tsv")
    ap.add_argument("--refmap", default="refbias/p1/refmap.tsv")
    ap.add_argument("--refs", default="refbias/build/7713a8d71d8e/refs")
    ap.add_argument("--crossmap-dir", default="is6110/assets/isclean_matched")
    ap.add_argument("--gff-dir", default="is6110/assets/matched_gff")
    ap.add_argument("--h37rv", default="refbias/build/7713a8d71d8e/refs/GCF_000195955.fasta")
    ap.add_argument("--h37rv-contig", default="NC_000962.3")
    ap.add_argument("--build-id", default="7713a8d71d8e")
    ap.add_argument("--outdir", default="refbias/p1i/vcf")
    ap.add_argument("--keys-out", default="is6110/results/p1i_cohort_keys.tsv")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    ref_of = {r["sample"]: r["reference"] for r in
              csv.DictReader(open(a.refmap), delimiter="\t")}
    flank = {}
    for r in csv.DictReader(open(a.flank), delimiter="\t"):
        flank[(r["sample"], int(r["r_pos"]))] = r

    sites = collections.defaultdict(list)
    for r in csv.DictReader(open(a.reconcile), delimiter="\t"):
        if not r["geometry"]:
            continue                       # chrom_side_only rows have no stack
        # Two axes, per is6110/docs/PROMOTED_TIER.md.  Evidence is whether
        # both element termini were seen, nothing else.  The old tier B also
        # required chromosome-side agreement, which is corroboration and not
        # evidence, so it is now carried in its own INFO field instead.  PASS
        # therefore means two_sided, and a two-sided call the chromosome-side
        # scan did not confirm is PASS with IS6110_CHROM_SIDE recording that.
        evidence = ("two_sided"
                    if r["geometry"] in ("tsd", "two_sided_wide")
                    else "one_sided")
        sites[r["sample"]].append(dict(r, evidence=evidence))

    seqs, exc, medlen, hcontig = {}, {}, {}, None
    _, h37seq = read_fasta_one(a.h37rv)
    keys_rows = []
    today = datetime.date.today().strftime("%Y%m%d")
    n_alt = n_ref = n_filt = 0

    # EVERY isolate in the refmap gets a file, including one with no surviving
    # site. Iterating `sites` alone silently omitted 3 isolates from scale200
    # and 2 from scale100, and an absent file cannot be told apart from an
    # isolate that was never processed -- the same conflation of "no data" with
    # "nothing there" that this pipeline has been careful about everywhere else.
    #
    # Those five are all lineage 1 with about one element copy, and all their
    # junction stacks fall below --min-reads-q 10, so having no ALT record is
    # the correct answer. Having no file is not. (An earlier note in the
    # handoff blamed their reference carrying zero annotated intervals; that
    # was a correlation -- a lineage-1 reference paired with lineage-1 isolates
    # -- and not the cause.)
    for sample in sorted(ref_of):
        ref = ref_of[sample]
        if ref not in seqs:
            seqs[ref] = read_fasta_one(os.path.join(a.refs, f"{ref}.fasta"))
            exc[ref] = load_excisions(os.path.join(a.crossmap_dir, f"{ref}.crossmap.tsv"))
            L = gff_lengths(os.path.join(a.gff_dir, f"{ref}.is6110.gff"))
            medlen[ref] = int(statistics.median(L)) if L else 1355
        contig, gseq = seqs[ref]
        exl, exr = exc[ref]

        recs_r, recs_h = [], []
        for s in sorted(sites.get(sample, []), key=lambda x: int(x["orig_pos"])):
            p = int(s["orig_pos"])
            shared = (p in exl) or (p in exr)
            fr = flank.get((sample, p), {})
            v = fr.get("verdict", "")
            hstate = ("empty" if v == "placed_h37rv_empty" else
                      "occupied" if v == "placed_h37rv_occupied" else "unplaced")
            hpos = fr.get("h37rv_pos", "") if hstate != "unplaced" else ""
            node = f"{fr.get('node','')}:{fr.get('node_offset','')}" if fr.get("node") else ""
            info = [f"SVTYPE=INS", f"SVLEN={medlen[ref]}",
                    f"MEINFO=IS6110,1,{medlen[ref]},+",
                    f"IS6110_EVIDENCE={s['evidence']}",
                    f"IS6110_CLASS={'ref_shared' if shared else 'ref_lacking'}",
                    f"IS6110_CHROM_SIDE={s.get('chrom_side','') or 'none'}",
                    f"IS6110_READS={s['reads_q']}",
                    f"IS6110_SPAN={s['span']}",
                    f"IS6110_SVLEN_SRC=ref_median",
                    f"IS6110_H37RV_STATE={hstate}", f"IS6110_BUILD={a.build_id}"]
            if hpos:
                info.append(f"IS6110_H37RV={hpos}")
                if fr.get("ismapper") != "":
                    info.append(f"IS6110_ISM={fr.get('ismapper')}")
            if node:
                info.append(f"IS6110_NODE={node}")
            filt = "PASS" if s["evidence"] == "two_sided" else "LowSupport"

            keys_rows.append(dict(
                sample=sample, reference=ref, build_id=a.build_id,
                evidence=s["evidence"],
                site_class="ref_shared" if shared else "ref_lacking", r_pos=p,
                state="REF" if shared else "ALT",
                frame="h37rv" if hpos else "node",
                key=f"h37rv:{hpos}" if hpos else f"node:{node}",
                h37rv_pos=hpos, h37rv_state=hstate, node=node,
                reads=s["reads_q"], ismapper=fr.get("ismapper", "")))

            if shared:
                n_ref += 1
                continue                   # REF in this frame; no record, by decision
            n_alt += 1
            n_filt += (filt != "PASS")
            base = gseq[p - 1] if 0 < p <= len(gseq) else "N"
            recs_r.append((p, f"{contig}\t{p}\t.\t{base}\t<INS:ME:IS6110>\t.\t{filt}\t"
                              f"{';'.join(info)};END={p}\tGT\t1"))
            if hpos:
                hb = h37seq[int(hpos) - 1] if 0 < int(hpos) <= len(h37seq) else "N"
                recs_h.append((int(hpos),
                    f"{a.h37rv_contig}\t{hpos}\t.\t{hb}\t<INS:ME:IS6110>\t.\t{filt}\t"
                    f"{';'.join(info)};END={hpos}\tGT\t1"))

        def write(path, ctg, clen, rows):
            with open(path, "w") as fh:
                fh.write("##fileformat=VCFv4.2\n")
                fh.write(f"##fileDate={today}\n")
                fh.write(f"##source=is6110_write_vcf.py\n")
                fh.write(f"##reference={ctg}\n")
                fh.write(f"##contig=<ID={ctg},length={clen}>\n")
                for h in HDR_COMMON:
                    fh.write(h + "\n")
                fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                         f"{sample}\n")
                for _, line in sorted(rows):
                    fh.write(line + "\n")

        write(os.path.join(a.outdir, f"{sample}.is6110.vcf"), contig, len(gseq), recs_r)
        write(os.path.join(a.outdir, f"{sample}.is6110.h37rv.vcf"),
              a.h37rv_contig, len(h37seq), recs_h)

    with open(a.keys_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(keys_rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(keys_rows)

    unplaced = sum(1 for k in keys_rows if k["frame"] == "node")
    print(f"  {len(keys_rows)} sites over {len(sites)} isolates")
    print(f"    ALT records written (matched-reference frame) : {n_alt}"
          f"   of which LowSupport {n_filt}")
    print(f"    REF in that frame, no record written          : {n_ref}")
    print(f"    carried on a graph node key, no H37Rv position: {unplaced}")
    print(f"\n  per sample: {a.outdir}/<sample>.is6110.vcf        (authoritative)")
    print(f"              {a.outdir}/<sample>.is6110.h37rv.vcf  (derived)")
    print(f"  cohort keys: {a.keys_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
