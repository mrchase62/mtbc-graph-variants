#!/usr/bin/env bash
#SBATCH --job-name=snp_tree
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 0-24:00
#SBATCH -p sapphire
#SBATCH --mem=64G
#SBATCH --output=slurm/snptree_%j.out
#SBATCH --error=slurm/snptree_%j.err
#
# Build a maximum-likelihood tree from a SNP-only alignment.
#
#   sbatch bin/build_snp_tree.sh <alignment.fasta> [outgroup] [prefix]
#
# The alignment holds only variable sites, so an ascertainment-bias correction
# (+ASC) is required -- without it branch lengths are inflated because the model
# expects to have seen the invariant majority of the genome.
set -euo pipefail

_mtb_env=""
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh" "${MTB_WORK:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
[[ -n "$_mtb_env" ]] || { echo "FATAL: cannot locate config/project_env.sh" >&2; exit 1; }
source "$_mtb_env"

: "${MTB_IQTREE:=/n/boslfs02/LABS/sfortune_lab/Lab/software/iqtree-3.0.1-Linux-intel/bin/iqtree3}"
[[ $# -ge 1 ]] || { mtb_usage "${BASH_SOURCE[0]}"; exit 2; }
ALN="$1"
OUTGROUP="${2:-GCF_035581225}"          # M. canettii
PREFIX="${3:-${ALN%.fasta}}"
mtb_require_file "$MTB_IQTREE" "$ALN"

echo "### alignment : $ALN"
echo "### outgroup  : $OUTGROUP"
"$MTB_IQTREE" -s "$ALN" \
    -m GTR+F+ASC+G4 \
    -B 1000 -alrt 1000 \
    -o "$OUTGROUP" \
    -T "${SLURM_CPUS_PER_TASK:-24}" \
    --prefix "$PREFIX" \
    -redo

# IQ-TREE's -o places the outgroup but still writes a trifurcating root, and the
# ASR and loss scripts want a rooted, bifurcating newick. Rooting was previously
# a manual step, which is how the 412 tree ended up rooted on a different taxon
# than the -o argument it was built with. Do it here instead.
"$MTB_PY" - "${PREFIX}.treefile" "$OUTGROUP" "${PREFIX}.rooted.nwk" <<'PY'
import sys, dendropy
tf, og, out = sys.argv[1], sys.argv[2], sys.argv[3]
t = dendropy.Tree.get(path=tf, schema="newick", preserve_underscores=True)
node = t.find_node_with_taxon_label(og)
if node is None:
    sys.exit(f"FATAL: outgroup {og} is not a tip in {tf}")
t.reroot_at_edge(node.edge, update_bipartitions=True,
                 suppress_unifurcations=True)
t.write(path=out, schema="newick", suppress_rooting=True,
        unquoted_underscores=True)
n = len(t.leaf_nodes()); k = len(t.seed_node.child_nodes())
print(f"### rooted on {og}: {n} tips, root children={k} -> {out}")
PY

echo "### done -> ${PREFIX}.treefile"
