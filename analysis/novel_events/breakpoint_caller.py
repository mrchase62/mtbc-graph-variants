#!/usr/bin/env python3
"""Prototype breakpoint-evidence caller (Phase A): structural events between
a sample and its matched reference R, from the P2 alignment of the sample's
reads to R.

Evidence (MAPQ >= 20, primary alignments):
  clipped reads   soft clip >= 10 bp; the clip side and position, and the
                  split-read partner from the SA tag
  CIGAR events    D or I of 50 bp or more inside a read
  pairs           wrong orientation, or more than --max-insert apart, or the
                  mate unmapped (the read anchors one side of new sequence)
Clustering: clip positions within 5 bp on the same side form a cluster. It
needs max(4, 10% of median depth) reads, and must lie more than 200 bp from
either contig end (the circular origin; v2).
Clip filter (v3, for real reads; off with --no-clip-filter):
  - each clip is quality-trimmed: its length is counted outward from the
    alignment end up to the first base below Q20, and must still be >= 10 bp;
  - at least half of the cluster's reads clip within 1 bp of its main
    position (a real breakpoint clips every read at the same base; untrimmed
    read ends and library chimeras clip at scattered positions);
  - the clipped sequences agree: at least 60% of the reads at the main
    position differ from the consensus of their first 20 clipped bases at no
    more than 10% of those bases.
  Every cluster, kept or not, is listed in <out>.clusters.tsv with these
  measures.
Events:
  DEL   right-clip cluster at A, left-clip cluster at B > A + 50, joined by
        >= 2 split reads or >= 3 pairs spanning A..B; or >= 2 reads with a
        CIGAR deletion of the same span
  INS   right-clip and left-clip clusters within 20 bp of each other and not
        joined to each other; or >= 2 reads with a CIGAR insertion there.
        The size is unknown here, and comes from local assembly.
  INV   split reads joining two places on opposite strands
  INS   (tandem duplication) reads clipped at A continue at a left-clip
        cluster B < A on the same strand: the segment B..A is repeated;
        scored as an insertion of A - B bp
  DEL   (depth-joined, v2) a right-clip cluster at A and a left-clip cluster
        at B within 5 kb, with no split reads but depth between them below
        30% of the median: a deletion in sequence too repetitive for SA tags
  REARR split reads joining two places on the same strand out of order
        (B < A), or more than 100 kb apart, not explained above
  BND   one clip cluster with no partner: unresolved breakpoint
Depth: for DEL the mean depth between A and B, against the sample median,
is reported, not used as a filter.

Output: one row per event, in R coordinates (1-based).
"""
import argparse
import collections
import statistics
import sys

import pysam

MINQ = 20
MINCLIP = 10
TOL = 5
JOIN_TOL = 20


def clip_out(r, n, right, filt):
    """The clipped bases read outward from the alignment end, quality-trimmed
    at the first base below Q20 (v3). None if fewer than MINCLIP remain."""
    seq, q = r.query_sequence, r.query_qualities
    if right:
        bases, quals = seq[-n:], (q[-n:] if q is not None else None)
    else:
        bases, quals = seq[:n][::-1], (q[:n][::-1] if q is not None else None)
    if filt and quals is not None:
        k = next((i for i, x in enumerate(quals) if x < MINQ), len(quals))
        bases = bases[:k]
    return bases if len(bases) >= MINCLIP else None


def agreement(seqs):
    """Fraction of clipped sequences within 10% mismatches of the consensus
    of their first 20 bases."""
    if not seqs:
        return 0.0
    cons = []
    for i in range(20):
        col = collections.Counter(x[i] for x in seqs if len(x) > i)
        if not col:
            break
        cons.append(col.most_common(1)[0][0])
    ok = 0
    for x in seqs:
        m = min(len(x), len(cons))
        ok += sum(1 for i in range(m) if x[i] != cons[i]) <= 0.1 * m
    return ok / len(seqs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--max-insert", type=int, default=1000)
    ap.add_argument("--no-clip-filter", action="store_true", help="v2 behaviour")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    filt = not a.no_clip_filter
    bam = pysam.AlignmentFile(a.bam)
    ctg = bam.references[0]
    L = bam.lengths[0]

    rclip = collections.defaultdict(list)   # pos -> [sa partners]
    lclip = collections.defaultdict(list)
    cig_del = collections.Counter()           # (start, end)
    cig_ins = collections.Counter()           # pos
    pairs = []                                # (left pos, right pos) for long pairs
    one_end = collections.Counter()           # pos // 50 for mate-unmapped anchors
    for r in bam.fetch(ctg):
        if (r.is_secondary or r.is_supplementary or r.is_unmapped
                or r.mapping_quality < MINQ):
            continue
        ct = r.cigartuples
        sa = []
        if r.has_tag("SA"):
            for x in r.get_tag("SA").rstrip(";").split(";"):
                c, p, st, cg, mq, _ = x.split(",")
                if int(mq) >= MINQ and c == ctg:
                    sa.append((int(p), st))
        st = "+" if not r.is_reverse else "-"
        if ct[-1][0] == 4 and ct[-1][1] >= MINCLIP:
            k = clip_out(r, ct[-1][1], right=True, filt=filt)
            if k is not None:
                rclip[r.reference_end].append((sa, st, k))
        if ct[0][0] == 4 and ct[0][1] >= MINCLIP:
            k = clip_out(r, ct[0][1], right=False, filt=filt)
            if k is not None:
                lclip[r.reference_start + 1].append((sa, st, k))
        p = r.reference_start
        for op, n in ct:
            if op == 2 and n >= 50:
                cig_del[(p + 1, p + n)] += 1
            if op == 1 and n >= 50:
                cig_ins[p] += 1
            if op in (0, 2, 3, 7, 8):
                p += n
        if r.is_paired and r.is_read1:
            if r.mate_is_unmapped:
                one_end[r.reference_start // 50] += 1
            elif (r.next_reference_name == ctg and not r.is_reverse
                    and r.mate_is_reverse and r.template_length > a.max_insert):
                pairs.append((r.reference_end, r.next_reference_start + 1))
    # median depth, sampled
    dep = []
    for x in range(1000, L - 1000, max(1, L // 2000)):
        dep.append(bam.count(ctg, x, x + 1))
    med = statistics.median(dep) if dep else 0
    minr = max(4, 0.1 * med)

    def clusters(d, side):
        out, cur = [], None
        for p in sorted(d):
            if cur and p - cur["end"] <= TOL:
                cur["end"] = p
                cur["reads"] += d[p]
            else:
                if cur:
                    out.append(cur)
                cur = dict(start=p, end=p, reads=list(d[p]))
        if cur:
            out.append(cur)
        for c in out:
            c["n"] = len(c["reads"])
            c["pos"] = c["start"] if c["end"] - c["start"] < 1 else \
                max(range(c["start"], c["end"] + 1), key=lambda q: len(d.get(q, [])))
            near_main = [k for q in range(c["pos"] - 1, c["pos"] + 2) for _, _, k in d.get(q, [])]
            c["main_frac"] = len(near_main) / c["n"]
            c["agree"] = agreement(near_main)
            c["keep"] = (c["n"] >= minr and 200 < c["pos"] < L - 200
                         and (not filt or (c["main_frac"] >= 0.5 and c["agree"] >= 0.6)))
        diag.extend((side, c) for c in out if c["n"] >= minr)
        # the circular origin: every read that crosses it is clipped at the
        # contig ends, so clusters within 200 bp of either end are dropped
        return [c for c in out if c["keep"]]

    diag = []
    R, Lc = clusters(rclip, "right"), clusters(lclip, "left")
    with open(a.out + ".clusters.tsv", "w") as fo:
        fo.write("side\tpos\treads\tmain_frac\tagree\tkept\n")
        for sd, c in sorted(diag, key=lambda x: x[1]["pos"]):
            fo.write(f"{sd}\t{c['pos']}\t{c['n']}\t{c['main_frac']:.2f}\t{c['agree']:.2f}\t{int(c['keep'])}\n")

    def depth_ratio(s, e):
        if not med or e < s:
            return None
        pts = [(s + e) // 2] if e - s < 40 else \
            range(s + 10, e - 10, max(1, (e - s - 20) // 20))
        return statistics.mean(bam.count(ctg, x - 1, x) for x in pts) / med
    used_r, used_l = set(), set()
    events = []

    def partners(c):
        return [(p, st, rst) for sa, rst, _ in c["reads"] for p, st in sa]

    # joins from right clips (A) via SA to a partner position
    for i, c in enumerate(R):
        A = c["pos"]
        ps = partners(c)
        same = [p for p, st, rst in ps if st == rst]
        opp = [p for p, st, rst in ps if st != rst]
        # DEL: partner at a left-clip cluster B > A + 50
        best = None
        for j, d in enumerate(Lc):
            B = d["pos"]
            if B <= A + 50:
                continue
            k = sum(1 for p in same if abs(p - B) <= JOIN_TOL + 150)
            span = sum(1 for x, y in pairs if x <= A + 20 and y >= B - 20 and x >= A - 600)
            if k >= 2 or span >= 3:
                if best is None or k + span > best[1]:
                    best = (j, k + span)
        if best is not None:
            j = best[0]
            B = Lc[j]["pos"]
            if B - A <= 100000:
                events.append(dict(type="DEL", start=A + 1, end=B - 1, size=B - A - 1,
                                   support=best[1], left=c["n"], right=Lc[j]["n"],
                                   source="clip"))
                used_r.add(i)
                used_l.add(j)
                continue
        if len(opp) >= 2:
            q = int(statistics.median(opp))
            events.append(dict(type="INV", start=min(A, q), end=max(A, q), size=abs(q - A),
                               support=len(opp), left=c["n"], right=0, source="clip"))
            used_r.add(i)
            continue
        # tandem duplication: reads clipped at the copy's end (A) continue at
        # its start (a left-clip cluster at B < A); scored as an insertion of
        # A - B bp at A
        dup = None
        for j, d in enumerate(Lc):
            B = d["pos"]
            if j in used_l or not (A - 100000 < B < A - 30):
                continue
            k = sum(1 for p in same if abs(p - B) <= JOIN_TOL + 150)
            if k >= 2 and (dup is None or k > dup[1]):
                dup = (j, k)
        if dup is not None:
            j = dup[0]
            B = Lc[j]["pos"]
            events.append(dict(type="INS", start=A, end=A, size=A - B, support=dup[1],
                               left=c["n"], right=Lc[j]["n"], source="clip;tandem_dup"))
            used_r.add(i)
            used_l.add(j)
            continue
        # deletion in a repeat: no split reads, but a left-clip cluster within
        # 5 kb after A with the depth between them below 30% of the median
        dd = None
        for j, d in enumerate(Lc):
            B = d["pos"]
            if j in used_l or not (A + 50 <= B <= A + 5000):
                continue
            dr = depth_ratio(A + 1, B - 1)
            if dr is not None and dr < 0.3:
                dd = (j, dr)
                break
        if dd is not None:
            j = dd[0]
            B = Lc[j]["pos"]
            events.append(dict(type="DEL", start=A + 1, end=B - 1, size=B - A - 1,
                               support=0, left=c["n"], right=Lc[j]["n"],
                               source="clip;depth_joined"))
            used_r.add(i)
            used_l.add(j)
            continue
        back = [p for p in same if p < A - 50 or p > A + 100000]
        if len(back) >= 2:
            q = int(statistics.median(back))
            events.append(dict(type="REARR", start=min(A, q), end=max(A, q), size=abs(q - A),
                               support=len(back), left=c["n"], right=0, source="clip"))
            used_r.add(i)
    # INS: right and left clusters close together, not joined
    for i, c in enumerate(R):
        if i in used_r:
            continue
        for j, d in enumerate(Lc):
            if j in used_l:
                continue
            if abs(d["pos"] - c["pos"]) <= JOIN_TOL:
                anch = sum(one_end.get(x, 0) for x in range(c["pos"] // 50 - 10, c["pos"] // 50 + 11))
                events.append(dict(type="INS", start=min(c["pos"], d["pos"]),
                                   end=max(c["pos"], d["pos"]), size=0,
                                   support=c["n"] + d["n"], left=c["n"], right=d["n"],
                                   source=f"clip;one_end={anch}"))
                used_r.add(i)
                used_l.add(j)
                break
    # CIGAR events
    for (s, e), n in cig_del.items():
        if n >= 2:
            events.append(dict(type="DEL", start=s, end=e, size=e - s + 1, support=n,
                               left=0, right=0, source="cigar"))
    for p, n in cig_ins.items():
        if n >= 2:
            events.append(dict(type="INS", start=p, end=p, size=0, support=n,
                               left=0, right=0, source="cigar"))
    # unpaired clusters
    for i, c in enumerate(R):
        if i not in used_r:
            events.append(dict(type="BND", start=c["pos"], end=c["pos"], size=0,
                               support=c["n"], left=c["n"], right=0, source="rclip"))
    for j, d in enumerate(Lc):
        if j not in used_l:
            events.append(dict(type="BND", start=d["pos"], end=d["pos"], size=0,
                               support=d["n"], left=0, right=d["n"], source="lclip"))
    # depth inside deletions
    for e in events:
        e["depth_ratio"] = ""
        if e["type"] == "DEL" and med:
            mid = [(e["start"] + e["end"]) // 2]
            if e["size"] > 200:
                mid = range(e["start"] + 50, e["end"] - 50, max(1, (e["size"] - 100) // 20))
            d = [bam.count(ctg, x - 1, x) for x in mid]
            e["depth_ratio"] = f"{statistics.mean(d) / med:.2f}"
    cols = ["sample", "type", "start", "end", "size", "support", "left", "right",
            "depth_ratio", "source"]
    with open(a.out, "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for e in sorted(events, key=lambda x: x["start"]):
            fo.write("\t".join(str(e.get(k, a.sample if k == "sample" else ""))
                               for k in cols) + "\n")
    c = collections.Counter(e["type"] for e in events)
    print(f"{a.sample}: median depth {med}, min cluster {minr:.0f}, "
          f"{len(R)} right / {len(Lc)} left clip clusters, events {dict(c)}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
