#!/usr/bin/env python3
"""Quick cohort tree for the NEW 1K: biallelic H37Rv-frame SNPs of the NEW
merged VCF (CLASS=small, outside the build repeat mask), variable among the
1,000 isolates; GT 0 -> REF base, 1 -> ALT base, anything else -> N. Writes
out/snps.fasta; FastTree -nt -gtr then builds out/new1k.fasttree.nwk.
"""
import subprocess, os, sys
import numpy as np
R = "/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants"
B = f"{R}/runroot/refbias/build/3893a3762ed1"
V = f"{R}/runroot/refbias/new1k_p200/p5/merged.vcf.gz"
BCF = os.environ["MTB_BCFTOOLS"]
mask = np.zeros(4411600, bool)
for l in open(f"{B}/assets/repeat_mask.bed"):
    c, s, e = l.split("\t")[:3]; mask[int(s):int(e)] = True
samples = subprocess.run([BCF, "query", "-l", V], capture_output=True, text=True).stdout.split()
p = subprocess.Popen([BCF, "query", "-r", "NC_000962.3", "-i", 'INFO/CLASS="small" && TYPE="snp" && N_ALT=1',
                      "-f", "%POS\t%REF\t%ALT[\t%GT]\n", V], stdout=subprocess.PIPE, text=True)
cols = []
for l in p.stdout:
    f = l.rstrip("\n").split("\t")
    pos = int(f[0])
    if mask[pos - 1] or len(f[1]) != 1 or len(f[2]) != 1: continue
    g = f[3:]
    na = sum(x == "1" for x in g); nr = sum(x == "0" for x in g)
    if na < 1 or nr < 1: continue
    m = {"0": f[1], "1": f[2]}
    cols.append("".join(m.get(x, "N") for x in g))
print(len(cols), "SNP columns", file=sys.stderr)
A = np.array([list(c) for c in cols]).T
with open(f"{R}/analysis/sv_eval_1k/out/snps.fasta", "w") as fo:
    for s, row in zip(samples, A): fo.write(f">{s}\n{''.join(row)}\n")
