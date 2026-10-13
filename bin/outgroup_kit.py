#!/usr/bin/env python3
"""Outgroup kit, route A: the event states of an outgroup genome that is not in
the graph, read from its assembly aligned directly to H37Rv. No graph is used.

Decided 2026-10-13 (HANDOFF OUT-1): ET1291 (GCF_035581225) and L8
(GCF_012923765). Their SNP states come from the existing per-position kit
(analysis/panel_checks/bin/genome_states.py); this script adds the rest.

    outgroup_kit.py --genome X.fna.gz --name X --h37rv H37Rv.fasta \
        --is6110 IS6110.fasta [--accessory catalogue.fasta catalogue.tsv] --outdir D

Steps, all against H37Rv (1-based coordinates in the outputs):
1. minimap2 -cx asm5 --cs, the assembly as query (as truth.py, against R).
   Writes X.paf.
2. paftools.js call -L 10000 -l 1000: X.var. `R` lines are the H37Rv
   stretches covered by exactly one alignment (an absence there is observed);
   `V` lines are SNPs and indels inside alignments.
3. X.sv.tsv: every difference of >= 50 bp. Indels inside alignments (from
   X.var) plus events between consecutive alignments, by truth.py's rules
   (DEL, INS, REPL, INV, REARR; the circular-origin wrap is dropped).
4. X.is6110.tsv: each IS6110 copy in the assembly (minimap2 of the element,
   >= 90% of its length), with 300 bp of flank on each side placed on H37Rv
   (minimap2 -x sr, MAPQ >= 20). Both flanks must place on the same strand;
   the H37Rv gap between them gives the state:
     |gap| <= 20           insertion site H37Rv lacks  -> site_absent_in_h37rv
     gap within 1355+-60   H37Rv has a copy here       -> shared_with_h37rv
     otherwise             unresolved
5. X.accessory.tsv (optional): each accessory catalogue locus aligned to the
   assembly (minimap2 -cx asm10); `cov` is the share of the locus covered.
   present >= 0.8, absent < 0.2, otherwise partial. A locus is a position,
   so a present hit must also sit at it: the assembly bases just outside the
   hit are projected onto H37Rv through X.paf, and one of them must land
   within 5 kb of the locus `pos` (the graph anchor of a locus amid other
   insertions can sit a few kb off); otherwise elsewhere (the sequence is in
   the genome, at another site). For copy_number loci (sequence H37Rv already
   carries) alignment cannot tell carriers apart, so their call is
   not_by_alignment.
"""
import argparse, collections, csv, gzip, os, subprocess, sys

MIN = 50
IS_LEN = 1355


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, check=True, **kw).stdout


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    seqs, name, buf = {}, None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def svs(paf, var):
    rows = []
    for line in open(var):
        c = line.rstrip("\n").split("\t")
        if c[0] != "V":
            continue
        ref = c[6] if c[6] != "-" else ""
        alt = c[7] if c[7] != "-" else ""
        if abs(len(ref) - len(alt)) < MIN and min(len(ref), len(alt)) < MIN:
            continue
        kind = "DEL" if len(ref) > len(alt) else "INS"
        st = int(c[2]) + 1
        rows.append((kind, st, max(int(c[3]), st), len(ref), len(alt), "cs"))
    rlen, alns = None, []
    for line in open(paf):
        if "tp:A:P" not in line:
            continue
        c = line.split("\t")
        rlen = int(c[6])
        alns.append(dict(q=c[0], qs=int(c[2]), qe=int(c[3]), strand=c[4],
                         ts=int(c[7]), te=int(c[8])))
    alns.sort(key=lambda x: (x["q"], x["qs"]))
    for x, y in zip(alns, alns[1:]):
        if x["q"] != y["q"]:
            continue
        qgap = y["qs"] - x["qe"]
        if x["strand"] != y["strand"]:
            kind, st, en, rg = "INV", min(x["te"], y["ts"]) + 1, max(x["te"], y["ts"]), 0
        else:
            if x["strand"] == "+":
                rg, st, en = y["ts"] - x["te"], x["te"] + 1, y["ts"]
            else:
                rg, st, en = x["ts"] - y["te"], y["te"] + 1, x["ts"]
            if rlen and abs(rg) > rlen / 2:
                continue
            if rg < -MIN:
                kind = "REARR"
            elif rg - qgap >= MIN and qgap < MIN:
                kind = "DEL"
            elif qgap - rg >= MIN and rg < MIN:
                kind = "INS"
            elif qgap >= MIN and rg >= MIN:
                kind = "REPL"
            else:
                continue
        rows.append((kind, st, max(st, en), max(rg, 0), max(qgap, 0), "break"))
    return sorted(rows, key=lambda r: r[1])


def is6110(a, asm, D):
    hits = run([a.minimap2, "-c", "-x", "asm10", "-N", "50", "--secondary=yes",
                os.path.join(D, "asm.fa"), a.is6110])
    copies = []
    for line in hits.splitlines():
        c = line.split("\t")
        if int(c[3]) - int(c[2]) >= 0.9 * int(c[1]):
            copies.append((c[5], int(c[7]), int(c[8]), c[4]))
    copies.sort()
    G = 300
    fq = []
    for i, (ctg, s, e, strand) in enumerate(copies):
        seq = asm[ctg]
        if s >= G:
            fq.append((f"c{i}_L", seq[s - G:s]))
        if e + G <= len(seq):
            fq.append((f"c{i}_R", seq[e:e + G]))
    fa = os.path.join(D, "is6110_flanks.fa")
    with open(fa, "w") as fo:
        for n, s in fq:
            fo.write(f">{n}\n{s}\n")
    place = {}
    for line in run([a.minimap2, "-c", "-x", "sr", "--secondary=no", a.h37rv, fa]).splitlines():
        c = line.split("\t")
        if int(c[11]) >= 20 and "tp:A:P" in line:
            place[c[0]] = (c[4], int(c[7]), int(c[8]))
    out = []
    for i, (ctg, s, e, strand) in enumerate(copies):
        L, R = place.get(f"c{i}_L"), place.get(f"c{i}_R")
        state, pos, gap = "unresolved", "", ""
        if L and R and L[0] == R[0]:
            # the base where the element sits, on H37Rv
            if L[0] == "+":
                gap, pos = R[1] - L[2], L[2]
            else:
                gap, pos = L[1] - R[2], R[2]
            if abs(gap) <= 20:
                state = "site_absent_in_h37rv"
            elif abs(gap - IS_LEN) <= 60:
                state = "shared_with_h37rv"
        elif not (L and R):
            state = "flank_unplaced"
        out.append((ctg, s + 1, e, strand, pos, gap, state))
    return out


def project(alns, ctg, q):
    """H37Rv coordinate of assembly base q, from the nearest primary alignment
    within 2 kb (linear within it)."""
    best = None
    for x in alns.get(ctg, ()):
        d = 0 if x["qs"] <= q < x["qe"] else min(abs(q - x["qs"]), abs(q - x["qe"]))
        if d <= 2000 and (best is None or d < best[0]):
            qq = min(max(q, x["qs"]), x["qe"] - 1)
            t = x["ts"] + (qq - x["qs"]) if x["strand"] == "+" else x["te"] - (qq - x["qs"])
            best = (d, t)
    return best[1] if best else None


def accessory(a, D, paf):
    alns = collections.defaultdict(list)
    for line in open(paf):
        if "tp:A:P" in line:
            c = line.split("\t")
            alns[c[0]].append(dict(qs=int(c[2]), qe=int(c[3]), strand=c[4],
                                   ts=int(c[7]), te=int(c[8])))
    fa, tsv = a.accessory
    info = {r["locus_id"]: r for r in csv.DictReader(open(tsv), delimiter="\t")}
    lens = {k: len(v) for k, v in read_fasta(fa).items()}
    cov, where = collections.defaultdict(list), collections.defaultdict(list)
    for line in run([a.minimap2, "-c", "-x", "asm10", "--secondary=no",
                     os.path.join(D, "asm.fa"), fa]).splitlines():
        c = line.split("\t")
        cov[c[0]].append((int(c[2]), int(c[3])))
        where[c[0]] += [project(alns, c[5], int(c[7]) - 1), project(alns, c[5], int(c[8]))]
    out = []
    for k, L in lens.items():
        iv, n, last = sorted(cov[k]), 0, 0
        for s, e in iv:
            s = max(s, last)
            if e > s:
                n += e - s
            last = max(last, e)
        f = n / L if L else 0
        r = info.get(k, {})
        nov = r.get("novelty", "")
        call = ("not_by_alignment" if nov == "copy_number" else
                "present" if f >= 0.8 else "absent" if f < 0.2 else "partial")
        if call == "present" and not any(t is not None and abs(t - int(r["pos"])) <= 5000
                                         for t in where[k]):
            call = "elsewhere"
        out.append((k, nov, L, round(f, 4), call))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--h37rv", required=True)
    ap.add_argument("--is6110", required=True)
    ap.add_argument("--accessory", nargs=2, metavar=("FASTA", "TSV"))
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--minimap2", default=os.environ.get("MTB_SV_MINIMAP2", "minimap2"))
    ap.add_argument("--k8", default=os.environ.get("MTB_K8", "k8"))
    ap.add_argument("--paftools", default=os.environ.get("MTB_PAFTOOLS", "paftools.js"))
    a = ap.parse_args()
    D = os.path.join(a.outdir, a.name)
    os.makedirs(D, exist_ok=True)
    asm = read_fasta(a.genome)
    with open(os.path.join(D, "asm.fa"), "w") as fo:
        for k, v in asm.items():
            fo.write(f">{k}\n{v}\n")

    paf = os.path.join(a.outdir, f"{a.name}.paf")
    with open(paf, "w") as fo:
        fo.write(run([a.minimap2, "-c", "-x", "asm5", "--cs", "-t", "2",
                      a.h37rv, os.path.join(D, "asm.fa")]))
    var = os.path.join(a.outdir, f"{a.name}.var")
    with open(var, "w") as fo:
        fo.write(run([a.k8, a.paftools, "call", "-L", "10000", "-l", "1000", "-"],
                     input=run(["sort", "-k6,6", "-k8,8n", paf])))

    rows = svs(paf, var)
    with open(os.path.join(a.outdir, f"{a.name}.sv.tsv"), "w") as fo:
        fo.write("type\th37rv_start\th37rv_end\tref_len\talt_len\tsource\n")
        for r in rows:
            fo.write("\t".join(map(str, r)) + "\n")

    isr = is6110(a, asm, D)
    with open(os.path.join(a.outdir, f"{a.name}.is6110.tsv"), "w") as fo:
        fo.write("contig\tcopy_start\tcopy_end\tstrand\th37rv_pos\th37rv_gap\tstate\n")
        for r in isr:
            fo.write("\t".join(map(str, r)) + "\n")

    acc = accessory(a, D, paf) if a.accessory else []
    if acc:
        with open(os.path.join(a.outdir, f"{a.name}.accessory.tsv"), "w") as fo:
            fo.write("locus\tnovelty\tlength\tcov\tcall\n")
            for r in acc:
                fo.write("\t".join(map(str, r)) + "\n")

    os.remove(os.path.join(D, "asm.fa"))
    nR = sum(int(l.split("\t")[3]) - int(l.split("\t")[2])
             for l in open(var) if l.startswith("R\t"))
    v = collections.Counter()
    for l in open(var):
        c = l.split("\t")
        if c[0] == "V":
            d = len(c[7].strip("-")) - len(c[6].strip("-"))
            v["snp" if len(c[6]) == 1 and len(c[7]) == 1 and "-" not in c[6] + c[7]
              else "indel<50" if abs(d) < MIN else "indel>=50"] += 1
    print(f"{a.name}: callable {nR:,} bp of H37Rv; {dict(v)}; "
          f"SV rows {dict(collections.Counter(r[0] for r in rows))}; "
          f"IS6110 {dict(collections.Counter(r[-1] for r in isr))}"
          + (f"; accessory {dict(collections.Counter(r[-1] for r in acc))}" if acc else ""))


if __name__ == "__main__":
    main()
