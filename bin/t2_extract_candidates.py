#!/usr/bin/env python3
"""Extract pure-insertion alleles from the graph VCF as accessory-panel candidates.

REF/ALT are reduced on both ends first, so what is emitted is the inserted
segment itself and not an arbitrary slice; only alleles that reduce to a pure
insertion (empty REF side) are kept. Complex records are excluded because their
inserted sequence is not well defined -- the same reduction the direct-repeat
screen uses (see DR_ELEMENTS.md).
"""
import argparse, subprocess, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--bcftools", required=True)
    ap.add_argument("--min-len", type=int, default=500)
    ap.add_argument("--out-fasta", default="refbias/t2/candidates.fasta")
    ap.add_argument("--out-meta", default="refbias/t2/candidates.tsv")
    a = ap.parse_args()

    hdr = subprocess.run([a.bcftools, "view", "-h", a.vcf],
                         capture_output=True, text=True, check=True)
    samples = hdr.stdout.strip().split("\n")[-1].split("\t")[9:]
    proc = subprocess.Popen([a.bcftools, "view", "-H", a.vcf],
                            stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    fa, meta, n = open(a.out_fasta, "w"), open(a.out_meta, "w"), 0
    meta.write("allele_id\tpos\tlen\tAC\tcarriers\n")
    for line in proc.stdout:
        f = line.rstrip("\n").split("\t")
        ref, alt = f[3].upper(), f[4].upper()
        if "," in alt or not set(ref + alt) <= set("ACGTN"):
            continue
        if len(alt) - len(ref) < a.min_len:
            continue
        p = 0
        while p < len(ref) and p < len(alt) and ref[p] == alt[p]:
            p += 1
        s = 0
        while (s < len(ref) - p and s < len(alt) - p
               and ref[len(ref) - 1 - s] == alt[len(alt) - 1 - s]):
            s += 1
        r, ins = ref[p:len(ref) - s], alt[p:len(alt) - s]
        if r or len(ins) < a.min_len:
            continue
        carriers = [x for x, g in zip(samples, f[9:]) if g == "1"]
        if not carriers:
            continue
        n += 1
        aid = f"A{n:05d}_{f[1]}_{len(ins)}"
        fa.write(f">{aid}\n{ins}\n")
        meta.write(f"{aid}\t{f[1]}\t{len(ins)}\t{len(carriers)}\t{','.join(carriers)}\n")
    proc.wait(); fa.close(); meta.close()
    print(f"  {n} pure-insertion alleles >= {a.min_len} bp -> {a.out_fasta}")


if __name__ == "__main__":
    sys.exit(main())
