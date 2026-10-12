#!/usr/bin/env python3
"""Homopolymer tracts in the graph's panel genomes (user, 2026-10-12: "we
need to redo the graph because of homopolymeric tracts in some of the
genomes"). Report only. Run from runroot with $MTB_PY.

For each of the build's 333 references (refbias/build/<build>/refs):
  pure runs   every run of one base of MIN_RUN (15) bp or more: natural MTBC
              homopolymers are short, so the distribution shows where
              artifacts begin (bin/assembly_qc_stats.py flags >= 100 bp)
  tracts      stretches where one base is >= PURITY (90%) of every WIN (50)
              bp window, windows merged: catches A-rich stretches that a pure
              run interrupts (GCF_045348265: 484 bp pure inside 1,243 bp)
Each pure run is placed on H37Rv by its 40 bp flanks (exact, unique match on
either strand), with H37Rv's own run length there: a locus where many
genomes carry runs of varying length is a natural variable tract; a run
found in one genome only, far longer than any other, is an artifact.
Then which cohort samples use a genome with a tract or a run >= 50 bp as
matched reference (every cohort in refbias/cohorts.tsv, P1 refmap).
Writes out/runs.tsv (every pure run >= 15 bp), out/tracts.tsv (every tract),
out/per_genome.tsv, out/loci.tsv, and prints a summary.
"""
import collections
import csv
import glob
import os
import re

import numpy as np

BUILD = "refbias/build/7713a8d71d8e-fix1"
O = "../analysis/homopolymer_scan/out"
MIN_RUN, WIN, PURITY = 15, 50, 0.9
RUN_RE = re.compile("|".join(f"{b}{{{MIN_RUN},}}" for b in "ACGTN"))


def read_fasta(p):
    seqs, name, buf = [], None, []
    for line in open(p):
        if line.startswith(">"):
            if name:
                seqs.append((name, "".join(buf).upper()))
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs.append((name, "".join(buf).upper()))
    return seqs


def tracts(seq):
    """Maximal stretches covered by WIN-bp windows in which one base (the
    same base throughout) is >= PURITY of the window."""
    out = []
    need = int(PURITY * WIN + 0.999)
    arr = np.frombuffer(seq.encode(), dtype=np.uint8)
    for b in "ACGT":
        c = np.concatenate(([0], np.cumsum(arr == ord(b))))
        ok = (c[WIN:] - c[:-WIN]) >= need          # window starting at i
        if not ok.any():
            continue
        d = np.diff(np.concatenate(([0], ok.view(np.int8), [0])))
        for i, j in zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)):
            out.append((int(i), int(j - 1 + WIN), b))   # last window ends at j-1+WIN
    return out


COMP = str.maketrans("ACGTN", "TGCAN")
FL = 40


def rc(x):
    return x.translate(COMP)[::-1]


def place(h37, seq, a, e, base):
    """(H37Rv 1-based position of the run's first base, H37Rv run length
    there, strand) from the run's flanks; None if neither flank is unique."""
    left, right = seq[max(0, a - FL):a], seq[e:e + FL]
    for fl, side in ((left, "L"), (right, "R")):
        if len(fl) < FL:
            continue
        for strand, q in (("+", fl), ("-", rc(fl))):
            i = h37.find(q)
            if i < 0 or h37.find(q, i + 1) >= 0:
                continue
            b = base if strand == "+" else base.translate(COMP)
            # the run sits after a forward left flank or a reverse right flank
            if (side, strand) in (("L", "+"), ("R", "-")):
                j = i + FL
                k = j
                while k < len(h37) and h37[k] == b:
                    k += 1
                return j + 1, k - j, strand
            j = i
            k = j
            while k > 0 and h37[k - 1] == b:
                k -= 1
            return k + 1, j - k, strand
    return None


def main():
    os.makedirs(O, exist_ok=True)
    h37 = read_fasta(os.environ["MTB_H37RV"])[0][1]
    runs, trs, per = [], [], []
    for p in sorted(glob.glob(f"{BUILD}/refs/*.fasta")):
        g = os.path.basename(p)[:-6]
        longest = (0, "", "", 0)
        for ctg, seq in read_fasta(p):
            for m in RUN_RE.finditer(seq):
                n = m.end() - m.start()
                pl = place(h37, seq, m.start(), m.end(), seq[m.start()])
                runs.append(dict(genome=g, contig=ctg, start=m.start() + 1, length=n, base=seq[m.start()],
                                 h37rv_pos=pl[0] if pl else "", h37rv_run=pl[1] if pl else "",
                                 strand=pl[2] if pl else ""))
                if n > longest[0]:
                    longest = (n, seq[m.start()], ctg, m.start() + 1)
            for a, e, b in tracts(seq):
                sub = seq[a:e]
                trs.append(dict(genome=g, contig=ctg, start=a + 1, end=e, length=e - a, base=b,
                                purity=round(sub.count(b) / len(sub), 3),
                                longest_pure=max(len(x) for x in re.findall(f"{b}+", sub))))
        gt = [t for t in trs if t["genome"] == g]
        per.append(dict(genome=g, longest_run=longest[0], longest_run_base=longest[1],
                        longest_run_pos=f"{longest[2]}:{longest[3]}" if longest[0] else "",
                        tracts=len(gt), longest_tract=max((t["length"] for t in gt), default=0)))
    for name, rows in (("runs", runs), ("tracts", trs), ("per_genome", per)):
        with open(f"{O}/{name}.tsv", "w") as fo:
            w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    print(f"genomes scanned: {len(per)}")
    dist = collections.Counter(r["longest_run"] for r in per)
    print("longest pure run per genome (bp: genomes): " +
          ", ".join(f"{k}: {v}" for k, v in sorted(dist.items())))
    print(f"\ntracts (one base >= {PURITY:.0%} of every {WIN} bp window): {len(trs)} in "
          f"{len({t['genome'] for t in trs})} genomes")
    print("genome\tcontig\tstart\tlength\tbase\tpurity\tlongest_pure")
    for t in sorted(trs, key=lambda t: -t["length"]):
        print(f"{t['genome']}\t{t['contig']}\t{t['start']}\t{t['length']}\t{t['base']}\t{t['purity']}\t{t['longest_pure']}")

    hit = {t["genome"] for t in trs} | {r["genome"] for r in per if r["longest_run"] >= 50}
    # loci: runs within 20 bp of each other on H37Rv
    placed = sorted((r for r in runs if r["h37rv_pos"] != ""), key=lambda r: r["h37rv_pos"])
    loci = []
    for r in placed:
        if loci and r["h37rv_pos"] - loci[-1]["last"] <= 20:
            loci[-1]["runs"].append(r)
            loci[-1]["last"] = r["h37rv_pos"]
        else:
            loci.append(dict(last=r["h37rv_pos"], runs=[r]))
    lrows = []
    for l in loci:
        rs = l["runs"]
        lens = sorted(x["length"] for x in rs)
        lrows.append(dict(h37rv_pos=rs[0]["h37rv_pos"], h37rv_run=rs[0]["h37rv_run"],
                          base=rs[0]["base"] if rs[0]["strand"] == "+" else rs[0]["base"].translate(COMP),
                          genomes=len({x["genome"] for x in rs}), min_len=lens[0],
                          median_len=lens[len(lens) // 2], max_len=lens[-1],
                          max_genome=max(rs, key=lambda x: x["length"])["genome"]))
    with open(f"{O}/loci.tsv", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(lrows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(lrows)
    print(f"\nruns >= {MIN_RUN} bp: {len(runs)}, placed on H37Rv {len(placed)}, at {len(lrows)} loci")
    print("H37Rv_pos\tH37Rv_run\tbase\tgenomes\tmin\tmedian\tmax\tlongest in")
    for l in sorted(lrows, key=lambda l: -l["genomes"]):
        print("\t".join(str(v) for v in l.values()))
    unplaced = [r for r in runs if r["h37rv_pos"] == ""]
    print(f"not placed (flanks absent or repeated in H37Rv): {len(unplaced)}; "
          f"longest {max((r['length'] for r in unplaced), default=0)}")

    reg = [l.rstrip("\n").split("\t") for l in open("refbias/cohorts.tsv") if not l.startswith("#")][1:]
    print("\ncohort samples whose matched reference has a tract or a run >= 50 bp")
    for f in reg:
        coh, root = f[0], f[3]
        rm = next((p for p in (f"{root}/p1/refmap.tsv", f"{root}/refmap.tsv") if os.path.exists(p)), None)
        if not rm:
            print(f"  {coh}: no refmap under {root}")
            continue
        use = collections.Counter(r["reference"] for r in csv.DictReader(open(rm), delimiter="\t")
                                  if r["reference"] in hit)
        if use:
            print(f"  {coh}: " + ", ".join(f"{g} {n}" for g, n in sorted(use.items())))


if __name__ == "__main__":
    main()
