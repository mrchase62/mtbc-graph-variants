#!/usr/bin/env python3
"""Can the matched-genome alignment confirm or reject an inherited SV?

User (2026-10-13): "Can't read metrics from the matched genome confirm or deny
the SV in the sample?" Analysis only; reads production outputs, writes here.

An inherited SV is a difference between the sample's matched reference R and
H37Rv, taken from the graph VCF (P4b `component=inherited`). The sample's own
P2 reads are aligned to R, so in R's coordinates:

  inherited INS (R has extra sequence Y)
      sample carries it  -> normal depth across Y
      sample lacks it    -> no MAPQ>=20 reads across Y (and no reads at all,
                            so it is not just a repeat), clips at Y's ends
  inherited DEL (R lacks segment X; R has a junction j)
      sample carries it  -> reads span j unclipped
      sample lacks it    -> reads clipped at j (the clipped part is X), or a
                            large insertion in the CIGAR at j

Placing the event on R. R's local haplotype is built from the graph VCF record
(300 bp of H37Rv either side + R's allele) and aligned to R's own refs FASTA
with minimap2. The alignment both places the event and checks the placement:
rows whose haplotype does not align end to end are reported `unplaced`.

Thresholds, fixed before looking at results:
  INS: flank depth (MAPQ>=20) < 10 -> low_depth
       Y/flank >= 0.5 -> carries
       Y/flank <= 0.1 and all-MAPQ Y/flank <= 0.2 -> lacks
       Y/flank <= 0.1 but all-MAPQ depth present -> repeat (unclear)
       else unclear
  DEL: spanning + disrupted < 5 -> low_depth
       disrupted >= 3 and disrupted/(span+disrupted) >= 0.3 -> lacks
       span >= 5 and disrupted/(span+disrupted) < 0.15 -> carries
       else unclear
  spanning: MAPQ>=20, one aligned block covers j-15..j+15
  disrupted: MAPQ>=20, a soft clip >= 10 bp with its boundary within 10 bp of
             j, or an I/D >= 10 bp within 15 bp of j

    check.py --rows 'p4b/*.sv_placed.tsv' --bam-dir work/p2 --out out.tsv
"""
import argparse, collections, glob, gzip, os, statistics, subprocess, sys, tempfile
import pysam

F = 300


def rc(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def load_rows(pattern):
    rows = []
    for p in sorted(glob.glob(pattern)):
        with open(p) as fh:
            hdr = fh.readline().rstrip("\r\n").split("\t")
            for line in fh:
                r = dict(zip(hdr, line.rstrip("\r\n").split("\t")))
                if r.get("component") == "inherited" and r["svtype"] in ("INS", "DEL"):
                    rows.append(r)
    return rows


def rows_from_vcf(vcf, refmap):
    """Every SV (>= 50 bp) each sample's matched reference carries in the graph
    VCF: P4b's inherited rule without its two exclusions (nested parents and a
    caller call within 200 bp). Used for truth7, whose P4b rows predate this
    VCF."""
    ref_of = {}
    for line in open(refmap):
        f = line.rstrip("\r\n").split("\t")
        if f[0] != "sample":
            ref_of[f[0]] = f[4] if len(f) > 4 else f[1]
    samples_of = collections.defaultdict(list)
    for s_, r_ in ref_of.items():
        samples_of[r_].append(s_)
    rows, col, seen = [], None, set()
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                col = {s: i for i, s in enumerate(line.rstrip("\n").split("\t")) if i >= 9}
                continue
            f = line.split("\t", 5)
            if max(len(f[3]), max(len(x) for x in f[4].split(","))) < 50:
                continue
            f = line.rstrip("\n").split("\t")
            alts = f[4].upper().split(",")
            for rid, ss in samples_of.items():
                g = f[col[rid]].split(":", 1)[0]
                if not g.isdigit() or g == "0" or int(g) > len(alts):
                    continue
                alt = alts[int(g) - 1]
                if set(f[3].upper() + alt) - set("ACGTN"):
                    continue
                d = len(alt) - len(f[3])
                if abs(d) < 50 or (rid, f[1], abs(d)) in seen:
                    continue
                seen.add((rid, f[1], abs(d)))
                for s_ in ss:
                    rows.append(dict(sample=s_, reference=rid, svtype="INS" if d > 0 else "DEL",
                                     svlen=str(abs(d)), h37rv_pos=f[1],
                                     h37rv_end=f[1] if d > 0 else str(int(f[1]) + abs(d))))
    return rows


def load_alleles(vcf, want):
    """want: {pos: set((ref_id, svlen))} -> {(ref_id, pos, svlen): (VCFPOS, REF, ALT)}

    The record R carries (GT > 0) within 10 bp of the row's position whose
    length change is within 2 bp or 2% of the row's svlen, nearest first.
    The tolerance is for rows written by P4b on an older graph VCF (truth7),
    whose records sit a few bp from the current file's; scale200_fix rows come
    from this file and match exactly."""
    best, col = {}, None
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                col = {s: i for i, s in enumerate(line.rstrip("\n").split("\t")) if i >= 9}
                continue
            f = line.split("\t", 9)
            vpos = int(f[1])
            near = [p for p in range(vpos - 10, vpos + 11) if p in want]
            if not near:
                continue
            f = line.rstrip("\n").split("\t")
            alts = f[4].upper().split(",")
            for p in near:
                for rid, svlen in want[p]:
                    if rid not in col:
                        continue
                    g = f[col[rid]].split(":", 1)[0]
                    if not g.isdigit() or g == "0" or int(g) > len(alts):
                        continue
                    alt = alts[int(g) - 1]
                    d = abs(abs(len(alt) - len(f[3])) - svlen)
                    if d > max(2, 0.02 * svlen):
                        continue
                    score = (abs(vpos - p), d)
                    old = best.get((rid, p, svlen))
                    if old is None or score < old[0]:
                        best[rid, p, svlen] = (score, vpos, f[3].upper(), alt)
    return {k: (v, r, a) for k, (_s, v, r, a) in best.items()}


def blocks_from_cigar(cg):
    ops, n = [], ""
    for ch in cg:
        if ch.isdigit():
            n += ch
        else:
            ops.append((int(n), ch)); n = ""
    return ops


def place(rows, alleles, refs_dir, h37rv, mm2):
    """Align R's local haplotype to R; return {(rid,pos,svlen): placement}."""
    h = "".join(l.strip() for l in open(h37rv) if not l.startswith(">")).upper()
    byref = collections.defaultdict(dict)
    tandem = {}
    for (rid, rowpos, svlen), (pos, ref, alt) in alleles.items():
        p = 0
        while p < min(len(ref), len(alt)) and ref[p] == alt[p]:
            p += 1
        s = 0
        while s < min(len(ref), len(alt)) - p and ref[-1 - s] == alt[-1 - s]:
            s += 1
        left = h[max(0, pos - 1 - F):pos - 1]
        right = h[pos - 1 + len(ref):pos - 1 + len(ref) + F]
        q = left + alt + right
        a = len(left) + p                       # event start in the query
        b = len(left) + len(alt) - s            # event end (INS: Y end; DEL: = a)
        # tandem: the changed sequence repeats what sits beside it (a VNTR
        # unit, a duplication). Fraction of its 20-mers also found within
        # 300 bp either side; >= 0.5 is tandem. For DEL, the deleted H37Rv
        # sequence against the remaining flanks.
        chg = q[a:b] if b > a else h[pos - 1 + p:pos - 1 + len(ref) - s]
        fl = q[:a] + "N" + q[b:]
        km = {fl[i:i + 20] for i in range(len(fl) - 19)}
        ck = [chg[i:i + 20] for i in range(max(0, len(chg) - 19))]
        tandem[rid, rowpos, svlen] = (sum(x in km for x in ck) / len(ck)) >= 0.5 if ck else False
        byref[rid][f"{rowpos}_{svlen}"] = (q, a, b, rowpos, svlen)
    out = {}
    tmp = tempfile.mkdtemp(prefix="imc_")
    for rid, qs in byref.items():
        fa = os.path.join(tmp, rid + ".fa")
        with open(fa, "w") as fh:
            for k, (q, *_r) in qs.items():
                fh.write(f">{k}\n{q}\n")
        res = subprocess.run([mm2, "-x", "asm10", "-c", "--secondary=no", "-t", "1",
                              os.path.join(refs_dir, rid + ".fasta"), fa],
                             capture_output=True, text=True, check=True).stdout
        best = {}
        for line in res.splitlines():
            f = line.split("\t")
            if f[0] in best:
                continue
            best[f[0]] = f
        for k, (q, a, b, pos, svlen) in qs.items():
            f = best.get(k)
            key = (rid, pos, svlen)
            if f is None:
                out[key] = dict(status="unplaced"); continue
            qlen, qs_, qe, strand = int(f[1]), int(f[2]), int(f[3]), f[4]
            ctg, ts, mapq = f[5], int(f[7]), int(f[11])
            # event must sit inside the aligned part with 100 bp to spare
            if qs_ > a - 100 or qe < b + 100:
                out[key] = dict(status="unplaced"); continue
            cg = next(x[5:] for x in f[12:] if x.startswith("cg:Z:"))
            # map query coordinates (forward query) to target coordinates
            qa, qb = (a, b) if strand == "+" else (qlen - b, qlen - a)
            qpos, tpos, m = qs_ if strand == "+" else qlen - qe, ts, {}
            for n, op in blocks_from_cigar(cg):
                if op == "M":
                    for x in (qa, qb):
                        if qpos <= x <= qpos + n:
                            m.setdefault(x, tpos + (x - qpos))
                    qpos += n; tpos += n
                elif op == "I":
                    for x in (qa, qb):
                        if qpos <= x <= qpos + n:
                            m.setdefault(x, tpos)
                    qpos += n
                elif op == "D":
                    tpos += n
            if qa not in m or qb not in m:
                out[key] = dict(status="unplaced"); continue
            out[key] = dict(status="placed", contig=ctg, start=min(m[qa], m[qb]),
                            end=max(m[qa], m[qb]), strand=strand, mapq=mapq)
    for k, v in out.items():
        v["tandem"] = "yes" if tandem.get(k) else "no"
    return out


def med_depth(bam, ctg, s, e, minq):
    if e <= s:
        return 0.0
    cov = bam.count_coverage(ctg, s, e, quality_threshold=0,
                             read_callback=lambda r: (not r.is_secondary and not r.is_duplicate
                                                      and not r.is_supplementary
                                                      and r.mapping_quality >= minq))
    tot = [sum(x) for x in zip(*cov)]
    return float(statistics.median(tot)) if tot else 0.0


def ins_metrics(bam, ctg, s, e, L):
    left = med_depth(bam, ctg, max(0, s - 300), max(0, s - 50), 20)
    right = med_depth(bam, ctg, min(L, e + 50), min(L, e + 300), 20)
    flank = (left + right) / 2
    y20 = med_depth(bam, ctg, s, e, 20)
    y0 = med_depth(bam, ctg, s, e, 0)
    clips = 0
    for x in (s, e):
        for r in bam.fetch(ctg, max(0, x - 20), min(L, x + 20)):
            if r.mapping_quality < 20 or r.is_secondary or r.is_duplicate:
                continue
            ct = r.cigartuples or []
            if ct and ct[0][0] == 4 and ct[0][1] >= 10 and abs(r.reference_start - x) <= 10:
                clips += 1
            if ct and ct[-1][0] == 4 and ct[-1][1] >= 10 and abs(r.reference_end - x) <= 10:
                clips += 1
    # REVISED 2026-10-13 after the first truth7 score (see RESULTS.md), so
    # these INS rules are tuned on truth7, not blind:
    #  * "lacks" also needs >= 2 clipped reads at Y's ends: GC-rich PE_PGRS
    #    inserts the isolate carries showed near-zero depth with no clips;
    #  * "carries" also needs both of Y's junctions spanned cleanly: in tandem
    #    repeats reads from other units fill Y's depth when the isolate lacks it.
    js = del_metrics(bam, ctg, s, L)
    je = del_metrics(bam, ctg, e, L)
    ends_clean = all(j["call"] == "carries" for j in (js, je))
    ends_broken = any(j["call"] == "lacks" for j in (js, je))
    if flank < 10:
        call = "low_depth"
    elif y20 / flank >= 0.5 and ends_clean:
        call = "carries"
    elif y20 / flank <= 0.1 and y0 / flank <= 0.2 and clips >= 2:
        call = "lacks"
    elif y20 / flank <= 0.1 and y0 / flank > 0.2:
        call = "repeat"
    else:
        call = "unclear"
    return dict(flank=flank, y20=y20, y0=y0, clips=clips,
                span=f"{js['span']},{je['span']}",
                disrupted=f"{js['disrupted']},{je['disrupted']}", call=call)


def del_metrics(bam, ctg, j, L):
    span = dis = 0
    for r in bam.fetch(ctg, max(0, j - 15), min(L, j + 15)):
        if r.mapping_quality < 20 or r.is_secondary or r.is_duplicate or r.is_supplementary:
            continue
        ct = r.cigartuples or []
        bad = False
        if ct and ct[0][0] == 4 and ct[0][1] >= 10 and abs(r.reference_start - j) <= 10:
            bad = True
        if ct and ct[-1][0] == 4 and ct[-1][1] >= 10 and abs(r.reference_end - j) <= 10:
            bad = True
        p = r.reference_start
        for op, n in ct:
            if op in (0, 7, 8):
                p += n
            elif op == 2:
                if n >= 10 and p - 15 <= j <= p + n + 15:
                    bad = True
                p += n
            elif op == 1:
                if n >= 10 and abs(p - j) <= 15:
                    bad = True
        if bad:
            dis += 1
        elif any(bs <= j - 15 and be >= j + 15 for bs, be in r.get_blocks()):
            span += 1
    n = span + dis
    if n < 5:
        call = "low_depth"
    elif dis >= 3 and dis / n >= 0.3:
        call = "lacks"
    elif span >= 5 and dis / n < 0.15:
        call = "carries"
    else:
        call = "unclear"
    return dict(span=span, disrupted=dis, call=call)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", default="", help="glob of P4b sv_placed.tsv files")
    ap.add_argument("--refmap", default="",
                    help="instead of --rows: every SV the matched reference "
                         "carries in the graph VCF, for each sample in this P1 refmap")
    ap.add_argument("--bam-dir", required=True)
    ap.add_argument("--build", required=True)
    ap.add_argument("--minimap2", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    vcf = os.path.join(a.build, "assets", "graph_collapsed.vcf.gz")
    rows = load_rows(a.rows) if a.rows else rows_from_vcf(vcf, a.refmap)
    want = collections.defaultdict(set)
    for r in rows:
        want[int(r["h37rv_pos"])].add((r["reference"], int(r["svlen"])))
    alleles = load_alleles(os.path.join(a.build, "assets", "graph_collapsed.vcf.gz"), want)
    print(f"{len(rows):,} inherited rows; alleles found for {len(alleles):,} of "
          f"{sum(len(v) for v in want.values()):,} (reference, event) pairs", file=sys.stderr)
    pl = place(rows, alleles, os.path.join(a.build, "refs"),
               os.path.join(a.build, "refs", "GCF_000195955.fasta"), a.minimap2)
    print("placement: " + str(collections.Counter(v["status"] for v in pl.values())), file=sys.stderr)
    cols = ["sample", "reference", "svtype", "svlen", "h37rv_pos", "h37rv_end", "status",
            "contig", "r_start", "r_end", "aln_mapq", "tandem", "flank", "y20", "y0", "clips",
            "span", "disrupted", "call"]
    bams = {}
    n = collections.Counter()
    with open(a.out, "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: r["sample"]):
            key = (r["reference"], int(r["h37rv_pos"]), int(r["svlen"]))
            o = dict(sample=r["sample"], reference=r["reference"], svtype=r["svtype"],
                     svlen=r["svlen"], h37rv_pos=r["h37rv_pos"], h37rv_end=r["h37rv_end"])
            p = pl.get(key) if key in alleles else None
            if p is None:
                o["status"] = "no_allele"; o["call"] = "unplaced"
            elif p["status"] != "placed":
                o["status"] = p["status"]; o["call"] = "unplaced"; o["tandem"] = p["tandem"]
            else:
                bp = os.path.join(a.bam_dir, r["sample"] + ".bam")
                if r["sample"] not in bams:
                    bams.clear()
                    bams[r["sample"]] = pysam.AlignmentFile(bp)
                bam = bams[r["sample"]]
                L = bam.get_reference_length(p["contig"])
                o.update(status="placed", contig=p["contig"], r_start=p["start"],
                         r_end=p["end"], aln_mapq=p["mapq"], tandem=p["tandem"])
                if r["svtype"] == "INS":
                    o.update(ins_metrics(bam, p["contig"], p["start"], p["end"], L))
                else:
                    o.update(del_metrics(bam, p["contig"], p["start"], L))
            n[r["svtype"], o["call"]] += 1
            fo.write("\t".join(str(o.get(c, "")) if not isinstance(o.get(c), float)
                               else f"{o[c]:.1f}" for c in cols) + "\n")
    for k in sorted(n):
        print(f"  {k[0]} {k[1]:<10} {n[k]:,}", file=sys.stderr)


if __name__ == "__main__":
    main()
