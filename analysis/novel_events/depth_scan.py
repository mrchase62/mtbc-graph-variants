#!/usr/bin/env python3
"""Step 2 of the Phase A follow-up: a read-depth copy-number scan in the
matched reference's frame, for tandem-repeat expansions and contractions.
These leave no clipped reads, because junction reads align to the next
repeat unit.

1. samtools depth -a, all mapping qualities (tandem-repeat reads are often
   MAPQ 0), over the P2 BAM.
2. Mean depth in 100 bp windows, divided by the sample's median window depth.
3. Runs of 2 or more consecutive windows at ratio >= 1.4 are a gain (scored
   as INS); runs at <= 0.6 are a loss (scored as DEL). Windows within 1 kb
   of either contig end are skipped (the circular origin).
4. The size is estimated as |ratio - 1| x run length.
No GC correction: the reads are simulated without GC bias. Real reads
(Phase B) will need one.
"""
import argparse
import statistics
import subprocess
import sys

W = 100


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--samtools", required=True)
    ap.add_argument("--gain", type=float, default=1.4)
    ap.add_argument("--loss", type=float, default=0.6)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    p = subprocess.Popen([a.samtools, "depth", "-a", "-Q", "0", a.bam],
                         stdout=subprocess.PIPE, text=True)
    sums, L = {}, 0
    for line in p.stdout:
        c, pos, d = line.split("\t")
        pos = int(pos)
        sums[(pos - 1) // W] = sums.get((pos - 1) // W, 0) + int(d)
        L = max(L, pos)
    p.wait()
    nwin = (L + W - 1) // W
    win = [sums.get(i, 0) / W for i in range(nwin)]
    med = statistics.median(win)
    ratio = [x / med if med else 0 for x in win]
    edge = 1000 // W
    events = []
    i = edge
    while i < nwin - edge:
        r = ratio[i]
        kind = "INS" if r >= a.gain else "DEL" if r <= a.loss else None
        if not kind:
            i += 1
            continue
        j = i
        while j + 1 < nwin - edge and ((kind == "INS" and ratio[j + 1] >= a.gain)
                                      or (kind == "DEL" and ratio[j + 1] <= a.loss)):
            j += 1
        if j - i + 1 >= 2:
            mr = statistics.mean(ratio[i:j + 1])
            span = (j - i + 1) * W
            events.append(dict(type=kind, start=i * W + 1, end=(j + 1) * W, size=round(abs(mr - 1) * span),
                               depth_ratio=f"{mr:.2f}"))
        i = j + 1
    cols = ["sample", "type", "start", "end", "size", "support", "left", "right",
            "depth_ratio", "source"]
    with open(a.out, "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for e in events:
            fo.write("\t".join([a.sample, e["type"], str(e["start"]), str(e["end"]), str(e["size"]),
                                "", "", "", e["depth_ratio"], "depth"]) + "\n")
    print(f"{a.sample}: median window depth {med:.1f}, {len(events)} depth events", file=sys.stderr)


if __name__ == "__main__":
    main()
