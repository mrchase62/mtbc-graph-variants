#!/bin/bash
# odgi position timing on the test graphs and CX333 (user, 2026-10-12): how
# much does graph composition change the cost of the lookups the pipeline
# makes? Same queries on every graph, one thread, one run at a time.
#   graphs   A L1-4,7 (50), B A+L5+L6 (52), C A+canettii (51), D B sparse
#            (52), CX333 (333)
#   refs     GCF_009730215 (L1), GCF_000193185 (L2), GCF_013010385 (L4),
#            all in every graph
#   queries  r2h: 20,000 random positions on the reference -> H37Rv
#            (P4/P4b form); h2r: 20,000 random H37Rv positions -> the
#            reference (SV projection / P3 form); load: 10 positions, which
#            is graph loading alone
#   two repeats each; /usr/bin/time: elapsed, user+sys, max RSS
# Submit from the repository root:
#   sbatch analysis/graph_tests/odgi_timing.sbatch
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 06:00:00
#SBATCH -J odgi_timing
#SBATCH -o analysis/graph_tests/logs/odgi_timing_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source config/project_env.sh
T=analysis/graph_tests
Q=$T/timing/queries
mkdir -p $Q
CX=/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome/graphs/CX333.s10k.k23.K15/mtb.complex333.fasta.gz.f4f5ee2.11fba48.36e68b6.smooth.final.og
H=GCF_000195955#1#NC_000962.3
declare -A OG=([A]=$(ls $T/arm_A_L14/pggb/*.og) [B]=$(ls $T/arm_B_L14_L56/pggb/*.og)
               [C]=$(ls $T/arm_C_L14_canettii/pggb/*.og) [D]=$(ls $T/arm_D_L14_L56_sparse/pggb/*.og) [CX333]=$CX)
REFS="GCF_009730215 GCF_000193185 GCF_013010385"
# queries from the A graph's input lengths (the same assemblies in every graph)
"$MTB_PY" - "$T/arm_A_L14/in.fa.gz.fai" "$Q" $REFS <<'PY'
import random, sys
fai, q, refs = sys.argv[1], sys.argv[2], sys.argv[3:]
L = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}
h = "GCF_000195955#1#NC_000962.3"
for g in refs:
    p = next(k for k in L if k.startswith(g + "#"))
    rng = random.Random(g)
    for name, path in (("r2h", p), ("h2r", h)):
        pos = sorted(rng.randrange(1000, L[path] - 1000) for _ in range(20000))
        with open(f"{q}/{g}.{name}.txt", "w") as fo:
            fo.writelines(f"{path},{x},+\n" for x in pos)
        with open(f"{q}/{g}.{name}.load.txt", "w") as fo:
            fo.writelines(f"{path},{x},+\n" for x in pos[:10])
    open(f"{q}/{g}.path", "w").write(p + "\n")
PY
OUT=$T/timing/timing.tsv
printf 'graph\tref\tquery\trep\tn_out\telapsed_s\tcpu_s\tmax_rss_kb\n' > $OUT
for rep in 1 2; do for gr in A B C D CX333; do for g in $REFS; do
  P=$(cat $Q/$g.path)
  for qn in r2h h2r r2h.load h2r.load; do
    case $qn in r2h*) tgt=$H ;; h2r*) tgt=$P ;; esac
    /usr/bin/time -f "%e\t%U\t%S\t%M" -o $T/timing/t.tmp \
      "$MTB_ODGI" position -i "${OG[$gr]}" -F $Q/$g.$qn.txt -r "$tgt" -t 1 > $T/timing/o.tmp 2>/dev/null
    n=$(grep -vc '^#' $T/timing/o.tmp || true)
    awk -F'\t' -v g=$gr -v r=$g -v q=$qn -v rep=$rep -v n=$n '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%s\n", g, r, q, rep, n, $1, $2+$3, $4}' $T/timing/t.tmp >> $OUT
  done
done; done; done
rm -f $T/timing/t.tmp $T/timing/o.tmp
echo "done: $(($(wc -l < $OUT) - 1)) runs"
