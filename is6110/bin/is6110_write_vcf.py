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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from is6110_seam import Seams


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
    """The reference's removed spans. A site within is6110_seam.SEAM_SLOP of
    one is its own copy (ref_shared); the rule is shared with promotion and
    flank placement. It was an exact match on the seam base here, which made
    a site 1-3 bp off the seam, inside the target-site duplication, an ALT."""
    return Seams(path)


def cluster_keys(obs, window):
    """{(sample, r_pos): canonical position} for each coordinate axis.

    One insertion seen in several isolates can be placed a few bases apart --
    the stack peak sits on either copy of the target-site duplication, or a
    base of wobble -- and each placement became its own key (review 4.5).
    Distances between neighbouring distinct H37Rv keys in gwas1000 (2026-09-30):

        1 bp 204   2 bp 137   3 bp 74   4 bp 46   5 bp 40   6 bp 29   7 bp 16
        15-24 bp: about 12 per bp, the background of distinct nearby insertions

    Merging at a distance is right when the excess over background is larger
    than the background, which holds to 6 bp (29 against 12) and fails at 7
    (16 against 12), hence --key-window 6. No two keys within 11 bp ever
    occur in the same isolate.

    The rule: sort the distinct positions; a cluster starts at its first
    position, and each next position joins the earliest cluster whose START is
    within `window` (so a chain of 1 bp steps cannot grow without bound) and
    holds none of the isolates carrying it -- two sites in one isolate are two
    insertions. Otherwise it starts a new cluster. The canonical position is the one carried by the most
    isolates, ties to the smallest. Clustering is across the cohort, so a key
    can differ between cohorts by a few bases; within a cohort it is fixed.

    obs: iterable of (axis, pos, sample, r_pos); positions are clustered within
    an axis only (one axis for H37Rv, one per graph node).
    """
    by = collections.defaultdict(lambda: collections.defaultdict(set))
    site = collections.defaultdict(list)
    for axis, pos, sample, r_pos in obs:
        by[axis][pos].add(sample)
        site[(axis, pos)].append((sample, r_pos))
    out = {}
    for axis, carriers in by.items():
        clusters = []
        for pos in sorted(carriers):
            # the earliest cluster still open (start within window) that no
            # isolate carrying this position is already in
            for c in clusters[-(window + 2):] if window > 0 else []:
                if pos - c["start"] <= window and not (carriers[pos] & c["samples"]):
                    c["pos"].append(pos); c["samples"] |= carriers[pos]
                    break
            else:
                clusters.append(dict(start=pos, pos=[pos], samples=set(carriers[pos])))
        for c in clusters:
            canon = min(c["pos"], key=lambda q: (-len(carriers[q]), q))
            for q in c["pos"]:
                for k in site[(axis, q)]:
                    out[k] = canon
    return out


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
    '##INFO=<ID=IS6110_H37RV,Number=1,Type=Integer,Description="H37Rv coordinate of this insertion: flank placement, both flanks required to agree, then the canonical position of placements within --key-window bp across the cohort">',
    '##INFO=<ID=IS6110_H37RV_PLACED,Number=1,Type=Integer,Description="This isolate\'s own flank-placed H37Rv coordinate, present only where it differs from IS6110_H37RV">',
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
    # NO BUILD DEFAULTS: --refs, --h37rv and --build-id were build
    # 7713a8d71d8e's (CX333), and bin/p1i_vcf.sh did not pass them, so every
    # build's IS6110 records were stamped 7713a8d71d8e. p1i_vcf.sh passes the
    # resolved build's.
    ap.add_argument("--refs", required=True, help="<build>/refs")
    ap.add_argument("--crossmap-dir", default="is6110/assets/isclean_matched")
    ap.add_argument("--gff-dir", default="is6110/assets/matched_gff")
    ap.add_argument("--h37rv", required=True,
                    help="<build>/refs/GCF_000195955.fasta")
    ap.add_argument("--h37rv-contig", default="NC_000962.3")
    ap.add_argument("--build-id", required=True)
    ap.add_argument("--outdir", default="refbias/p1i/vcf")
    ap.add_argument("--keys-out", default="is6110/results/p1i_cohort_keys.tsv")
    ap.add_argument("--key-window", type=int, default=6,
                    help="merge keys of one insertion placed up to this many "
                         "bases apart across isolates; 0 disables. See "
                         "cluster_keys for how 6 was measured")
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

    seqs, exc, medlen = {}, {}, {}
    hcontig, h37seq = read_fasta_one(a.h37rv)
    # The derived VCF's CHROM is --h37rv-contig, but its bases come from this
    # FASTA. They agree today (NC_000962.3); make a mismatch an error rather
    # than a VCF whose CHROM names a different sequence from its REF bases.
    if hcontig.split()[0] != a.h37rv_contig:
        sys.exit(f"FATAL: --h37rv-contig {a.h37rv_contig!r} does not match the "
                 f"H37Rv FASTA's contig {hcontig.split()[0]!r}")
    # Each site's H37Rv position or graph node, exactly as placed, then one
    # canonical position per insertion across the cohort (cluster_keys).
    # A node is an identity only where the carrier's own path visits it once
    # (audit P3IS-3): node 46966, a 1 bp node every path walks 85-453 times,
    # gathered insertions over 1 Mb apart into one key. node_occ is that visit
    # count, from is6110_project_sites.py; a table without it predates the
    # check and cannot say which node keys are safe.
    if flank and "node_occ" not in next(iter(flank.values())):
        sys.exit(f"FATAL: {a.flank} has no node_occ column; rerun "
                 f"is6110_project_sites.py and is6110_place_by_flank.py")

    def node_ok(fr):
        return bool(fr.get("node")) and fr.get("node_occ") == "1"

    def placement(sample, p):
        fr = flank.get((sample, p), {})
        v = fr.get("verdict", "")
        hstate = ("empty" if v == "placed_h37rv_empty" else
                  "occupied" if v == "placed_h37rv_occupied" else "unplaced")
        hpos = fr.get("h37rv_pos", "") if hstate != "unplaced" else ""
        return fr, hstate, hpos
    obs_h, obs_n = [], []
    for sample in ref_of:
        for s in sites.get(sample, []):
            p = int(s["orig_pos"])
            fr, _, hpos = placement(sample, p)
            if hpos:
                obs_h.append(("h37rv", int(hpos), sample, p))
            elif node_ok(fr):
                obs_n.append((fr["node"], int(fr.get("node_offset") or 0), sample, p))
    canon_h = cluster_keys(obs_h, a.key_window)
    canon_n = cluster_keys(obs_n, a.key_window)
    n_moved = (sum(1 for o in obs_h if canon_h[(o[2], o[3])] != o[1])
               + sum(1 for o in obs_n if canon_n[(o[2], o[3])] != o[1]))

    keys_rows = []
    today = datetime.date.today().strftime("%Y%m%d")
    n_alt = n_ref = n_filt = n_repeat = 0

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

        recs_r, recs_h = [], []
        for s in sorted(sites.get(sample, []), key=lambda x: int(x["orig_pos"])):
            p = int(s["orig_pos"])
            shared = exc[ref].shared(p)
            fr, hstate, hpos_raw = placement(sample, p)
            node_raw = (f"{fr.get('node','')}:{fr.get('node_offset','')}"
                        if fr.get("node") else "")
            hpos, node = hpos_raw, node_raw
            if hpos_raw:
                hpos = str(canon_h[(sample, p)])
            elif node_ok(fr):
                node = f"{fr['node']}:{canon_n[(sample, p)]}"
            # MEINFO's polarity is the element's strand, which no step of the
            # arm determines (review 4.10), so it is written unknown rather
            # than the '+' it used to claim for every site.
            info = [f"SVTYPE=INS", f"SVLEN={medlen[ref]}",
                    f"MEINFO=IS6110,1,{medlen[ref]},.",
                    f"IS6110_EVIDENCE={s['evidence']}",
                    f"IS6110_CLASS={'ref_shared' if shared else 'ref_lacking'}",
                    f"IS6110_CHROM_SIDE={s.get('chrom_side','') or 'none'}",
                    f"IS6110_READS={s['reads_q']}",
                    f"IS6110_SPAN={s['span']}",
                    f"IS6110_SVLEN_SRC=ref_median",
                    f"IS6110_H37RV_STATE={hstate}", f"IS6110_BUILD={a.build_id}"]
            if hpos:
                info.append(f"IS6110_H37RV={hpos}")
                if hpos != hpos_raw:
                    info.append(f"IS6110_H37RV_PLACED={hpos_raw}")
                if fr.get("ismapper") != "":
                    info.append(f"IS6110_ISM={fr.get('ismapper')}")
            if node:
                info.append(f"IS6110_NODE={node}")
            filt = "PASS" if s["evidence"] == "two_sided" else "LowSupport"

            # THE COHORT KEY IS IN THE H37Rv OR NODE FRAME, NOT THE MATCHED ONE
            # (audit P3IS-1). site_class says whether the isolate's matched
            # reference also holds this copy, which makes it REF in the matched
            # frame; carried into the cohort as REF it became GT 0 on an
            # H37Rv-frame <INS> record, indistinguishable from a non-carrier.
            # The state is now what the isolate carries against what the key's
            # frame holds:
            #   H37Rv empty     carrier ALT, whatever its matched reference holds
            #   H37Rv occupied  carrier REF: it has the copy H37Rv has. A
            #                   stage-2 non-carrier here lacks H37Rv's copy, which
            #                   an <INS> allele cannot express, so stage 2 leaves
            #                   it NOCALL (is6110_p5_stage2.py)
            #   graph node      carrier ALT: the node is flank sequence and holds
            #                   no element; a stage-2 REF means no junction there
            # site_class stays in the table as information about the matched
            # frame only.
            # A node the carrier's path visits more than once is no identity
            # (P3IS-3), nor one whose count is unknown: the site gets frame
            # repeat_node and no key, and P5 and the merged VCF leave it out,
            # counted.
            if hpos:
                frame, key = "h37rv", f"h37rv:{hpos}"
                state = "REF" if hstate == "occupied" else "ALT"
            elif node_ok(fr):
                frame, key, state = "node", f"node:{node}", "ALT"
            elif node_raw:
                frame, key, state = "repeat_node", "", "ALT"
                n_repeat += 1
            else:
                frame, key, state = "node", "node:", "ALT"   # refused below
            keys_rows.append(dict(
                sample=sample, reference=ref, build_id=a.build_id,
                evidence=s["evidence"],
                site_class="ref_shared" if shared else "ref_lacking", r_pos=p,
                state=state, frame=frame, key=key,
                h37rv_pos=hpos, h37rv_state=hstate,
                node=node if frame == "node" else "",
                reads=s["reads_q"], ismapper=fr.get("ismapper", ""),
                h37rv_pos_placed=hpos_raw, node_placed=node_raw,
                node_occ=fr.get("node_occ", "")))

            # The derived H37Rv-frame file follows the same rule: a record for
            # every carrier at an H37Rv-empty locus, shared or not; none at an
            # occupied one, where the carrier matches H37Rv (GT 0, omitted).
            if hpos and hstate != "occupied":
                hb = h37seq[int(hpos) - 1] if 0 < int(hpos) <= len(h37seq) else "N"
                recs_h.append((int(hpos),
                    f"{a.h37rv_contig}\t{hpos}\t.\t{hb}\t<INS:ME:IS6110>\t.\t{filt}\t"
                    f"{';'.join(info)};END={hpos}\tGT\t1"))
            if shared:
                n_ref += 1
                continue                   # REF in this frame; no record, by decision
            n_alt += 1
            n_filt += (filt != "PASS")
            base = gseq[p - 1] if 0 < p <= len(gseq) else "N"
            recs_r.append((p, f"{contig}\t{p}\t.\t{base}\t<INS:ME:IS6110>\t.\t{filt}\t"
                              f"{';'.join(info)};END={p}\tGT\t1"))

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

    # A site with neither an H37Rv position nor a graph node has no identity
    # in P5's key space. It used to be written as the key "node:", which
    # is6110_p5_merge.py then skipped without a word: 1,800 gwas1000 rows, 816
    # of them ALT, because p1i_vcf.sh ran the projection without --all-stacks.
    # Every stack gets a flank row under --all-stacks, so an empty key now
    # means an upstream step was skipped -- stop rather than drop sites.
    unkeyed = [k for k in keys_rows if not k["h37rv_pos"] and not k["node"]
               and k["frame"] != "repeat_node"]
    if unkeyed:
        ex = ", ".join(f'{k["sample"]}@{k["r_pos"]}' for k in unkeyed[:5])
        sys.exit(f"FATAL: {len(unkeyed)} sites have no H37Rv position and no "
                 f"graph node (e.g. {ex}); was is6110_project_sites.py run with "
                 f"--all-stacks?")
    # Two sites of one isolate on one key are two insertions merged into one
    # record, and the readers kept whichever row came last (audit P3IS-5).
    # cluster_keys never puts two sites of an isolate in one cluster, and a
    # node visited once cannot hold two, so this is a defect upstream: stop.
    seen = collections.Counter((k["sample"], k["key"]) for k in keys_rows if k["key"])
    dup = [sk for sk, n in seen.items() if n > 1]
    if dup:
        ex = ", ".join(f"{s}@{k}" for s, k in dup[:5])
        sys.exit(f"FATAL: {len(dup)} (sample, key) pairs hold more than one "
                 f"site (e.g. {ex})")

    KEY_FIELDS = ["sample", "reference", "build_id", "evidence", "site_class",
                  "r_pos", "state", "frame", "key", "h37rv_pos", "h37rv_state",
                  "node", "reads", "ismapper", "h37rv_pos_placed", "node_placed",
                  "node_occ"]
    with open(a.keys_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=KEY_FIELDS, delimiter="\t")
        w.writeheader(); w.writerows(keys_rows)

    unplaced = sum(1 for k in keys_rows if k["frame"] == "node")
    print(f"  {len(keys_rows)} sites over {len(sites)} isolates")
    print(f"    ALT records written (matched-reference frame) : {n_alt}"
          f"   of which LowSupport {n_filt}")
    print(f"    REF in that frame, no record written          : {n_ref}")
    print(f"    carried on a graph node key, no H37Rv position: {unplaced}")
    print(f"    on a node not visited once by its path, NOT keyed: {n_repeat}")
    print(f"    moved to their insertion's canonical key       : {n_moved}"
          f"   (--key-window {a.key_window})")
    print(f"\n  per sample: {a.outdir}/<sample>.is6110.vcf        (authoritative)")
    print(f"              {a.outdir}/<sample>.is6110.h37rv.vcf  (derived)")
    print(f"  cohort keys: {a.keys_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
