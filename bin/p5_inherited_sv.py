#!/usr/bin/env python3
"""P5: read each SV a sample inherits from its matched reference, in BOTH frames.

WHY (HANDOFF fix list ABS-2, with ABS-1 for the SV block). An inherited SV is
a difference between the sample's matched reference R and H37Rv. The sample's
reads are aligned to R, so if the sample shares it no caller sees it, and the
pipeline reported it on R's word:

  deletions   re-genotyped per sample by bin/p5_sv_genotype.py, but from the
              H37Rv frame only. The matched frame contributed "probes do not
              project and the flanks are covered", an inference: R's junction
              was never read. scale200_fix: 5,049 ALT confirmed, 3,119 REF
              contradicted, 4,479 NOCALL.
  insertions  never checked. The interval catalogue is deletions only, and the
              merged VCF wrote them from the caller matrix, presence-only, so
              every inherited insertion was present in every sample matched to
              that R.

THIS SCRIPT READS R's JUNCTION. For one sample it writes one row per catalogue
interval and one per inherited insertion (P4b `component=inherited`), with:

  matched frame  R's allele from the graph VCF is placed on R's own sequence
                 (minimap2 asm10 of R's local haplotype) and the sample's P2
                 reads are read there (bin/inherited_check.py): a deletion's
                 junction spanned cleanly ("carries") or disrupted ("lacks"),
                 an insertion's depth across Y and its clean junctions.
  H37Rv frame    for the catalogue intervals, R's deletion over the interval
                 checked in the sample's P1 alignment to H37Rv (ABS-1): reads
                 across both of its ends, or joining them. p5_sv_genotype.py
                 uses it where the interval does not project into R at all and
                 used to be stated ABSENT unread (the ctpV case: SAMEA2297133's
                 872 and 82 bp intervals, nested in R's 1,296 bp deletion).
  masked         the event's span, or either end within 150 bp, in the build's
                 repeat_mask.bed (the end rule kept, user 2026-10-14). The
                 HANDOFF asked for masked inherited events to be NOCALL until
                 adjudicated (ABS-2 item 5). Adjudicated on ctrl25_p200's 25
                 finished assemblies (2026-10-14): inside the mask the
                 combined two-frame rule was right at 432 deletion cells and
                 wrong at 16, so the user decided ("yes, accept all three")
                 that masked inherited DELETIONS are genotyped by that rule;
                 masked inherited INSERTIONS stay NOCALL "masked" (the matched
                 frame was 157 right / 32 wrong there).
  tandem         the changed sequence is a VNTR unit or duplication of its
                 flanks: undecidable from short reads, NOCALL "tandem"
                 (ABS-2 item 3; parked until a repeat-length genotyper is
                 chosen).

THE STATE COLUMN. For an insertion it is final, and the merged VCF writes it
instead of the caller matrix's presence-only ALT: "carries" -> ALT, anything
else NOCALL with its reason (masked, tandem, repeat, low_depth, unclear,
unplaced, no_allele, ins_lacks_unvalidated). "lacks" is NOT REF: on
ctrl25_p200 it was wrong at 11 of 12 cells (the isolate carried the insert),
and truth7 never saw a truly absent insertion called "lacks"; the user
decided 2026-10-14 ("yes, accept all three") to state it NOCALL
"ins_lacks_unvalidated" until it is validated. Expect many inherited insertions to
move from present to NOCALL. For a deletion interval it is the matched frame's
call alone; p5_sv_genotype.py combines it with the H37Rv frame (both agree ->
that state; one decides -> that state; they conflict -> NOCALL
"frames_conflict").

ONE GRAPH VCF PER RUN (ABS-2 item 6). R's alleles are read from the build's
graph_collapsed.vcf.gz, the panel table must have been made from that file
(its sha256 is in the table's header), and every P4b row must carry this
build's id; anything else is FATAL rather than matched by tolerance.

    p5_inherited_sv.py --sample S --reference R --intervals sv_intervals.tsv \\
        --sv-placed p4b/S.sv_placed.tsv --build <build> --bam p2/S.bam \\
        --h37rv-bam p1/S.h37rv.bam --minimap2 minimap2 --out S.inherited_sv.tsv
"""
import argparse, collections, csv, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import inherited_check as ic  # noqa: E402

COLUMNS = ["sample", "reference", "kind", "id", "start", "end", "svlen",
           "r_carries", "r_del_start", "r_del_end", "masked", "tandem",
           "placed", "contig", "r_start", "r_end", "span", "disrupted",
           "flank", "y20", "y0", "clips", "call", "h_junction",
           "h_bridge_start", "h_bridge_end", "h_joined", "state", "evidence"]


def read_fasta(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq:
                    break
                continue
            seq.append(line.strip())
    s = "".join(seq).upper()
    if len(s) < 4_000_000:
        sys.exit(f"FATAL: {path} is {len(s):,} bp; not an H37Rv FASTA")
    return s


def overlap_frac(a0, a1, b0, b1):
    shared = min(a1, b1) - max(a0, b0) + 1
    return max(0, shared) / max(1, max(a1 - a0 + 1, b1 - b0 + 1))


def ins_state(masked, tandem, call):
    """(state, evidence) for an inherited insertion."""
    if masked:
        return "NOCALL", "masked"
    if tandem:
        return "NOCALL", "tandem"
    if call == "carries":
        return "ALT", "matched_carries"
    if call == "lacks":
        # unvalidated (wrong 11 of 12 on ctrl25_p200): NOCALL, user 2026-10-14
        return "NOCALL", "ins_lacks_unvalidated"
    return "NOCALL", call


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--intervals", required=True,
                    help="the cohort's catalogue, bin/sv_intervals.py")
    ap.add_argument("--sv-placed", required=True,
                    help="P4b's <sample>.sv_placed.tsv, for the inherited insertions")
    ap.add_argument("--build", required=True,
                    help="the build directory: graph VCF, panel table, mask, refs")
    ap.add_argument("--bam", required=True,
                    help="the sample's P2 alignment to its matched reference")
    ap.add_argument("--h37rv-bam", required=True,
                    help="the sample's P1 alignment to H37Rv")
    ap.add_argument("--minimap2", required=True)
    ap.add_argument("--samtools", default=os.environ.get("MTB_SAMTOOLS", "samtools"))
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    import pysam
    B = a.build
    gvcf = os.path.join(B, "assets", "graph_collapsed.vcf.gz")
    bid = ""
    for line in open(os.path.join(B, "build_info.tsv")):
        k, _, v = line.rstrip("\n").partition("\t")
        if k == "build_id":
            bid = v
    if not bid:
        sys.exit(f"FATAL: no build_id in {B}/build_info.tsv")
    for p in (gvcf, a.intervals, a.sv_placed, a.bam):
        if not os.path.exists(p):
            sys.exit(f"FATAL: {p} does not exist")
    diffs = ic.PanelDiffs(os.path.join(B, "assets", "panel_differences.tsv.gz"),
                          a.reference, gvcf)
    mask = ic.Mask.from_bed(os.path.join(B, "assets", "repeat_mask.bed"))
    h37 = read_fasta(os.path.join(B, "refs", f"{ic.H37RV_ACC}.fasta"))
    h37f = ic.H37rvFrame(a.h37rv_bam, a.samtools)

    # ---- the catalogue: R's deletion over each interval ---------------------
    ivs = list(csv.DictReader(open(a.intervals, newline=""), delimiter="\t"))
    if not ivs:
        sys.exit(f"FATAL: no intervals in {a.intervals}")
    rows, events = [], {}
    for iv in ivs:
        s0, e0 = int(iv["start"]), int(iv["end"])
        carriers = set((iv.get("ref_carriers") or "").split(",")) - {""}
        hits = diffs.deletions_over(s0, e0)
        r_car = a.reference in carriers
        if r_car and hits:
            # the deletion of R's that this interval clusters: the best
            # reciprocal overlap, not merely the largest
            d = max(hits, key=lambda x: (overlap_frac(s0, e0, x[0], x[1]),
                                         -abs((x[1] - x[0]) - (e0 - s0))))
        else:
            d = diffs.deletion_over(s0, e0)
        msk = mask.site(s0, e0, (s0, e0)) or (d is not None and mask.site(s0, e0, d))
        row = dict(sample=a.sample, reference=a.reference, kind="DEL",
                   id=iv["interval"], start=s0, end=e0, svlen=iv.get("svlen", ""),
                   r_carries="yes" if r_car else "no", masked="yes" if msk else "no")
        if d is None:
            row["h_junction"] = ("masked" if msk else "r_no_genotype"
                                 if diffs.no_genotype_over(s0, e0) else "no_r_deletion")
            row["call"] = "no_r_deletion"
        else:
            row.update(r_del_start=d[0], r_del_end=d[1])
            if msk:
                row["h_junction"] = "masked"
            else:
                v, ev = h37f.junction(d[0], d[1])
                row.update(h_junction=v, h_bridge_start=ev["bridge_start"],
                           h_bridge_end=ev["bridge_end"], h_joined=ev["joined"])
            events[("DEL", iv["interval"])] = (d[2], d[3])
        rows.append(row)

    # ---- the inherited insertions, P4b's rows -------------------------------
    n_ins = collections.Counter()
    with open(a.sv_placed, newline="") as fh:
        for x in csv.DictReader(fh, delimiter="\t"):
            if x.get("component") != "inherited":
                continue
            if x.get("build_id") != bid:
                sys.exit(f"FATAL: {a.sv_placed} row {x.get('key')} is from build "
                         f"{x.get('build_id')!r}, not {bid}: P4b and P5 must read "
                         f"one graph VCF. Rerun P4b on this build.")
            if x["svtype"] != "INS":
                continue
            pos, ln = int(x["h37rv_pos"]), int(x["svlen"])
            if x.get("callable") in ("yes", "no"):
                msk = x["callable"] == "no"
            else:
                msk = mask.site(pos, pos + 1, (pos, pos + 1))
            row = dict(sample=a.sample, reference=a.reference, kind="INS",
                       id=x["key"], start=pos, end=pos, svlen=ln, r_carries="yes",
                       masked="yes" if msk else "no")
            al = diffs.ins.get((pos, ln))
            if al is None:
                row.update(call="no_allele", state="NOCALL", evidence="no_allele")
                n_ins["no allele in the panel table"] += 1
            else:
                events[("INS", x["key"])] = (pos, al)
            rows.append(row)

    # ---- place R's alleles on R, read the sample's P2 reads there -----------
    alleles = ic.fetch_alleles(gvcf, set(events.values()))
    todo = {k: (v[0],) + alleles[v] for k, v in events.items() if v in alleles}
    rfa = os.path.join(B, "refs", f"{a.reference}.fasta")
    placed = ic.place_on_r(todo, h37, rfa, a.minimap2, a.threads) if todo else {}
    bam = pysam.AlignmentFile(a.bam)
    lens = dict(zip(bam.references, bam.lengths))
    for row in rows:
        k = (row["kind"], row["id"])
        if k not in events:
            if row["kind"] == "DEL":
                row.update(state="NOCALL", evidence=row["call"])
            continue
        pl = placed.get(k)
        if pl is None:
            row.update(placed="no_allele", call="unplaced")
        else:
            row["tandem"] = pl["tandem"]
            if pl["status"] != "placed":
                row.update(placed="no", call="unplaced")
            else:
                if pl["contig"] not in lens:
                    sys.exit(f"FATAL: {a.bam} has no contig {pl['contig']}; it is "
                             f"not aligned to {rfa}")
                L = lens[pl["contig"]]
                row.update(placed="yes", contig=pl["contig"], r_start=pl["start"],
                           r_end=pl["end"])
                if row["kind"] == "DEL":
                    m = ic.matched_junction(bam, pl["contig"], pl["start"], L)
                    row.update(span=m["span"], disrupted=m["disrupted"], call=m["call"])
                else:
                    m = ic.matched_insertion(bam, pl["contig"], pl["start"], pl["end"], L)
                    row.update(flank=f"{m['flank']:.1f}", y20=f"{m['y20']:.1f}",
                               y0=f"{m['y0']:.1f}", clips=m["clips"],
                               span=m["span"], disrupted=m["disrupted"],
                               call=m["call"])
        if row["kind"] == "INS":
            st, why = ins_state(row["masked"] == "yes", row.get("tandem") == "yes",
                                row["call"])
            row.update(state=st, evidence=why)
            n_ins[f"{st} {why}"] += 1
        else:
            row.update(state=ic.MATCHED_STATE.get(row["call"], "NOCALL"),
                       evidence=row["call"])

    with open(a.out + ".tmp", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                           lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)
    os.replace(a.out + ".tmp", a.out)
    nd = collections.Counter((r["r_carries"], r.get("call", ""))
                             for r in rows if r["kind"] == "DEL" and r.get("r_del_start"))
    print(f"  {a.sample} (ref {a.reference}): {len(ivs):,} intervals, "
          f"{sum(nd.values()):,} with an R deletion over them; "
          f"{sum(1 for r in rows if r['kind'] == 'INS'):,} inherited insertions")
    print(f"    H37Rv-frame median depth {h37f.baseline:.0f}, bridging floor "
          f"{h37f.min_bridge:.0f} reads")
    for (rc, c), n in sorted(nd.items()):
        print(f"    DEL, R {'carries' if rc == 'yes' else 'overlaps'} the interval: "
              f"matched frame {c:<10} {n:>6,}")
    for k, n in sorted(n_ins.items()):
        print(f"    INS {k:<28} {n:>6,}")
    print(f"  -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
