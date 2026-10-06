#!/usr/bin/env python3
"""SNP and indel counts against H37Rv, and homopolymer-indel privacy, for every
candidate genome: the producer of the SNP-outlier screen's input and of the
long-read error checks (audit PGB-7, PGB-8).

    $MTB_PY bin/variant_counts.py --accessions data/qc/accessions.txt \\
        --assembly-dir data/rotated --h37rv data/ref/H37Rv.fasta \\
        --minimap2 "$MTB_MINIMAP2" --k8 "$MTB_K8" --paftools "$MTB_PAFTOOLS" \\
        --out data/qc/variants --threads 8

Ported from analysis/external_assemblies/bin/variant_counts.py (the 2026-10-04
external-assembly QC, which ran it on the 332 CX333 genomes as the baseline);
the counting is unchanged.

Per genome: minimap2 -cx asm5 --cs to H37Rv, then paftools call (-l 1000
-L 10000), keeping variants in uniquely aligned regions (depth 1). Counted:

  snps            single-base substitutions
  indel1          1 bp insertions and deletions
  indel1_hp       1 bp indels inside an H37Rv homopolymer run of >= 4 of that base
                  (insertion: the inserted base matches a run of >= 4 at the
                  site; deletion: the deleted base is part of a run of >= 4)
  indel2_49       2-49 bp indels
  aligned_bp      H37Rv bases covered (R lines)

PRIVACY. A homopolymer indel is "private" when no other genome in the run
carries the same (position, ref, alt). Real indels are mostly shared down a
clade; private homopolymer indels are the uncorrected long-read error
signature (CP010333: 417; a clean genome: 0 to 3). Privacy depends on the set
run together, so run every candidate in one call.

Outputs in --out:
  variant_counts.tsv   one row per genome, every count
  counts.snp.tsv       accession<TAB>snps         -> snp_outlier_screen.py
  counts.indel.tsv     accession<TAB>indel1       -> snp_outlier_screen.py --column indel1
  counts.hpindel.tsv   accession<TAB>indel1_hp_private
                                                   -> snp_outlier_screen.py --column indel1_hp_private
  lineages.tsv         barcode lineage (with --barcode), in the layout
                       snp_outlier_screen.py reads; tb-profiler's
                       lineages.all.tsv is the panel's authority

LINEAGE (optional, --barcode tbdb.barcode.bed). tb-profiler's own barcode applied
to the genome's calls: a marker counts as seen when its position is aligned, and
as carried when the genome has the marker allele. Lineage = the deepest label
with >= 50% of its markers carried whose ancestors are also supported.
"""
import argparse, bisect, collections, os, subprocess, sys
from multiprocessing import Pool

G = {}


def init(h37, mm2, k8, paft):
    G.update(h37=h37, mm2=mm2, k8=k8, paft=paft)


def run_one(job):
    sid, path = job
    p1 = subprocess.Popen([G["mm2"], "-cx", "asm5", "--cs", "-t", "1", G["h37"], path],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(["sort", "-k6,6", "-k8,8n"], stdin=p1.stdout, stdout=subprocess.PIPE)
    p3 = subprocess.run([G["k8"], G["paft"], "call", "-l", "1000", "-L", "10000", "-"],
                        stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                        universal_newlines=True)
    return (sid,) + parse_call(p3.stdout)


def parse_call(text):
    """paftools call output -> (regions, variants at depth 1)"""
    regions, var = [], []
    for l in text.splitlines():
        f = l.split("\t")
        if f[0] == "R":
            regions.append((int(f[2]), int(f[3])))
        elif f[0] == "V" and f[4] == "1":
            var.append((int(f[2]), int(f[3]), f[6].upper(), f[7].upper()))
    return regions, var


def hp_run(seq, pos0, base):
    """length of the run of `base` touching 0-based position pos0 in seq"""
    i = pos0
    while i > 0 and seq[i - 1] == base:
        i -= 1
    j = pos0
    while j < len(seq) and seq[j] == base:
        j += 1
    return j - i


def count_variants(var, h37):
    """(snps, indel1, indel1_hp, indel2_49, hp keys, {pos: alt} of SNPs)"""
    snps = indel1 = hp = i2 = 0
    keys, snp_at = set(), {}
    for s, e, ref, alt in var:
        if ref != "-" and alt != "-" and len(ref) == 1 and len(alt) == 1:
            snps += 1; snp_at[e] = alt
            continue
        L = max(len(ref.replace("-", "")), len(alt.replace("-", "")))
        if L == 1:
            indel1 += 1
            base = (alt if ref == "-" else ref).replace("-", "")
            run = hp_run(h37, s, base) if base else 0
            if ref == "-":   # insertion before 0-based s: run on either side
                run = max(hp_run(h37, s, base) if s < len(h37) and h37[s] == base else 0,
                          hp_run(h37, s - 1, base) if s > 0 and h37[s - 1] == base else 0)
            if run >= 4:
                hp += 1; keys.add((s, ref, alt))
        elif L < 50:
            i2 += 1
    return snps, indel1, hp, i2, keys, snp_at


def private_counts(hpkeys):
    """{sid: number of its hp keys carried by no other genome}"""
    cnt = collections.Counter(k for ks in hpkeys.values() for k in ks)
    return {sid: sum(1 for k in ks if cnt[k] == 1) for sid, ks in hpkeys.items()}


def barcode_lineage(regions, snp_at, h37, markers, by_lin):
    reg = sorted(regions)
    starts = [r[0] for r in reg]

    def covered(p):
        i = bisect.bisect_right(starts, p - 1) - 1
        return i >= 0 and reg[i][0] < p <= reg[i][1]
    seen, hit = collections.Counter(), collections.Counter()
    for pos, L, allele in markers:
        if not covered(pos):
            continue
        gt = snp_at.get(pos, h37[pos - 1])
        seen[L] += 1; hit[L] += (gt == allele)
    sup = {L for L in by_lin if seen[L] and hit[L] / seen[L] >= 0.5}
    ok = [L for L in sup if all(".".join(L.split(".")[:i]) in sup or ".".join(L.split(".")[:i]) not in by_lin
                                for i in range(1, L.count(".") + 1))]
    return max(ok, key=lambda L: (L.count("."), len(L))) if ok else ""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accessions", required=True)
    ap.add_argument("--assembly-dir", required=True)
    ap.add_argument("--assembly-suffix", default=".dnaA_rotated.fasta")
    ap.add_argument("--h37rv", required=True)
    ap.add_argument("--barcode", help="tb-profiler tbdb.barcode.bed, for a barcode lineage")
    ap.add_argument("--minimap2", default=os.environ.get("MTB_MINIMAP2"))
    ap.add_argument("--k8", default=os.environ.get("MTB_K8"))
    ap.add_argument("--paftools", default=os.environ.get("MTB_PAFTOOLS"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    for k in ("minimap2", "k8", "paftools"):
        if not getattr(a, k):
            ap.error(f"--{k} is required (or set MTB_{k.upper()})")
    os.makedirs(a.out, exist_ok=True)
    h37 = "".join(l.strip() for l in open(a.h37rv) if not l.startswith(">")).upper()
    accs = [l.strip() for l in open(a.accessions) if l.strip()]
    jobs, missing = [], []
    for acc in accs:
        p = os.path.join(a.assembly_dir, acc + a.assembly_suffix)
        (jobs.append((acc, p)) if os.path.exists(p) else missing.append(acc))
    if missing:
        sys.exit(f"no assembly for {len(missing)} accession(s): {' '.join(missing[:10])}")
    markers = []
    if a.barcode:
        for l in open(a.barcode):
            f = l.rstrip("\n").split("\t")
            if len(f) >= 5 and f[1].isdigit():
                markers.append((int(f[2]), f[3], f[4].upper()))
    by_lin = collections.Counter(L for _, L, _ in markers)

    rows, hpkeys = {}, {}
    with Pool(a.threads, initializer=init, initargs=(a.h37rv, a.minimap2, a.k8, a.paftools)) as pool:
        for n, (sid, regions, var) in enumerate(pool.imap_unordered(run_one, jobs), 1):
            snps, indel1, hp, i2, keys, snp_at = count_variants(var, h37)
            aligned = sum(e - s for s, e in regions)
            if not aligned:
                sys.exit(f"{sid}: nothing aligned to H37Rv; check the assembly")
            lin = barcode_lineage(regions, snp_at, h37, markers, by_lin) if markers else ""
            rows[sid] = dict(accession=sid, barcode_lineage=lin, aligned_bp=aligned, snps=snps,
                             indel1=indel1, indel1_hp=hp, indel2_49=i2)
            hpkeys[sid] = keys
            if n % 50 == 0:
                print(f"  {n}/{len(jobs)}", flush=True)
    for sid, n in private_counts(hpkeys).items():
        rows[sid]["indel1_hp_private"] = n
    cols = ["accession", "barcode_lineage", "aligned_bp", "snps", "indel1", "indel1_hp",
            "indel1_hp_private", "indel2_49"]
    with open(os.path.join(a.out, "variant_counts.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for sid in sorted(rows):
            fh.write("\t".join(str(rows[sid][c]) for c in cols) + "\n")
    for key, col in (("snp", "snps"), ("indel", "indel1"), ("hpindel", "indel1_hp_private")):
        with open(os.path.join(a.out, f"counts.{key}.tsv"), "w") as fh:
            for sid in sorted(rows):
                fh.write(f"{sid}\t{rows[sid][col]}\n")
    if markers:
        with open(os.path.join(a.out, "lineages.tsv"), "w") as fh:
            fh.write("accession\tstrain\n")
            for sid in sorted(rows):
                fh.write(f"{sid}\t{rows[sid]['barcode_lineage']}\n")
    print(f"  {len(rows)} genomes -> {a.out}")


if __name__ == "__main__":
    main()
