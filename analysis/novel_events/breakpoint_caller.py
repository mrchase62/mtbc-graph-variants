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
needs max(4, 10% of median depth) reads.
Events:
  DEL   right-clip cluster at A, left-clip cluster at B > A + 50, joined by
        >= 2 split reads or >= 3 pairs spanning A..B; or >= 2 reads with a
        CIGAR deletion of the same span
  INS   right-clip and left-clip clusters within 20 bp of each other and not
        joined to each other; or >= 2 reads with a CIGAR insertion there.
        The size is unknown here, and comes from local assembly.
  INV   split reads joining two places on opposite strands
  REARR split reads joining two places on the same strand out of order
        (B < A), or more than 100 kb apart: duplication or relocation
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--max-insert", type=int, default=1000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
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
        if ct[-1][0] == 4 and ct[-1][1] >= MINCLIP:
            rclip[r.reference_end].append((sa, "+" if not r.is_reverse else "-"))
        if ct[0][0] == 4 and ct[0][1] >= MINCLIP:
            lclip[r.reference_start + 1].append((sa, "+" if not r.is_reverse else "-"))
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

    def clusters(d):
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
        return [c for c in out if c["n"] >= minr]

    R, Lc = clusters(rclip), clusters(lclip)
    used_r, used_l = set(), set()
    events = []

    def partners(c):
        return [(p, st, rst) for sa, rst in c["reads"] for p, st in sa]

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
