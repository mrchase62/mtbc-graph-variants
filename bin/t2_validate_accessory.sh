#!/usr/bin/env bash
#SBATCH --job-name=T2_acc
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 0-02:00
#SBATCH -p shared
#SBATCH --mem=16G
#SBATCH --output=slurm/T2_%A_%a.out
#SBATCH --error=slurm/T2_%A_%a.err
#
# T2 of REFERENCE_BIAS_TESTS.md: validate the accessory catalogue against the
# assemblies, independently of the graph that produced it.
#
#   sbatch --array=1-332%40 bin/t2_validate_accessory.sh
#
# For every candidate accessory allele, blast it against one panel assembly and
# record the best hit. Aggregated across all 332 genomes this asks whether the
# graph's carrier set corresponds to genomes that actually contain the sequence.
#
# The test is NOT "is the sequence present" -- at a multiallelic locus every
# genome carries some version, so non-carriers hit a similar sequence too. The
# test is whether identity and coverage in graph-carriers are materially higher
# than in graph-non-carriers. If they are not, the allele call distinguishes
# nothing. This is the check that was skipped when the ~5 kb alleles at 2,268,725
# were wrongly dismissed as artefacts.
set -euo pipefail
_mtb_env=""
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh" \
          "${MTB_WORK:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
source "$_mtb_env"

# NO DEFAULTS. They were the pilot's: CX333's 332 genomes
# (loci/hely_tatc/samples.txt), its candidates and its hits folder, where a
# sample's existing hits file is kept -- so a new graph's run into it would
# have reused CX333's results. docs/PANEL_TREE.md section 4 gives the commands.
#   SAMPLES  the graph VCF's samples, one per line (bcftools query -l)
#   QUERY    t2_extract_candidates.py's --out-fasta
#   OUT      a hits folder for this graph alone
SAMPLES="${SAMPLES:?SAMPLES is unset: the graph VCF's samples (docs/PANEL_TREE.md section 4)}"
QUERY="${QUERY:?QUERY is unset: t2_extract_candidates.py's --out-fasta}"
OUT="${OUT:?OUT is unset: a hits folder for this graph}"
mkdir -p "$OUT"

# An explicit argument WINS over SLURM_ARRAY_TASK_ID. The other order
# silently ignores a caller-supplied index inside an array job, which
# made a targeted six-genome probe process samples 1-6 instead.
IDX="${1:-${SLURM_ARRAY_TASK_ID:-1}}"
SAMPLE=$(sed -n "${IDX}p" "$SAMPLES")
[[ -n "$SAMPLE" ]] || { echo "no sample at line ${IDX}" >&2; exit 1; }
ASM="data/rotated/${SAMPLE}.dnaA_rotated.fasta"
DEST="${OUT}/${SAMPLE}.hits.tsv.gz"
[[ -s "$DEST" ]] && { echo "[T2] ${SAMPLE}: already done"; exit 0; }
mtb_require_file "$ASM" "$QUERY"

# -max_hsps 1 keeps one alignment per query/subject pair; the question is whether
# the allele is there at full length, not how many times it occurs.
"${MTB_QC_BIN}/blastn" -query "$QUERY" -subject "$ASM" \
    -outfmt "6 qseqid pident length qlen mismatch gapopen" \
    -evalue 1e-20 -max_hsps 1 -dust no -num_threads "${SLURM_CPUS_PER_TASK:-4}" \
  | awk -v S="$SAMPLE" 'BEGIN{OFS="\t"}
      { if (!($1 in best) || $3 > bl[$1]) { best[$1]=$2; bl[$1]=$3; ql[$1]=$4 } }
      END { print "sample","allele_id","pident","length","qlen","cov"
            for (q in best)
              print S, q, best[q], bl[q], ql[q], bl[q]/ql[q] }' \
  | gzip -c > "${DEST}.tmp" && mv -f "${DEST}.tmp" "$DEST"

echo "[T2] ${SAMPLE}: done"
