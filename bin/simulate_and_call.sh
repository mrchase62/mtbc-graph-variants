#!/usr/bin/env bash
#
# Simulate Illumina reads from a source genome, map to a reference, call variants.
#
#   bin/simulate_and_call.sh <source.fasta> <reference.fasta> <outdir> <prefix> [depth] [readlen]
#
# The instrument for the reference-bias test program (REFERENCE_BIAS_TESTS.md).
# Reference bias is a property of the reference and the aligner, not of
# sequencing chemistry, so simulated reads isolate exactly the effect under study
# and give perfect per-base truth. Absolute sensitivities from this harness are
# optimistic; only RELATIVE comparisons between arms on identical reads are to be
# trusted, and the absolute offset is calibrated separately against real data.
#
# Calls haploid, matching the project convention and the MTBC literature.
set -euo pipefail
_mtb_env=""
for _c in "${MTB_ENV_FILE:-}" \
          "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)/config/project_env.sh" \
          "${SLURM_SUBMIT_DIR:-}/config/project_env.sh" \
          "${MTB_WORK:-}/config/project_env.sh"; do
    [[ -n "$_c" && -r "$_c" ]] && { _mtb_env="$_c"; break; }
done
[[ -n "$_mtb_env" ]] || { echo "FATAL: cannot locate config/project_env.sh" >&2; exit 1; }
source "$_mtb_env"

[[ $# -ge 4 ]] || { mtb_usage "${BASH_SOURCE[0]}"; exit 2; }
SRC="$1"; REF="$2"; OUTDIR="$3"; PREFIX="$4"
DEPTH="${5:-60}"; READLEN="${6:-150}"; FRAG="${7:-350}"

BWA="${MTB_BWA:-${MTB_QC_BIN}/bwa}"
WGSIM="${MTB_WGSIM:-${MTB_QC_BIN}/wgsim}"
SAMTOOLS="${MTB_SAMTOOLS}"
# GATK is run from the container the Liftover pipeline already uses. The conda
# gatk wrapper fails on a compute node with "/usr/bin/env: 'python': No such
# file or directory" -- it needs a python on PATH that is not there under sbatch.
GATK_SIF="${MTB_GATK_SIF:?MTB_GATK_SIF is unset; see config/project_env.sh}"
mtb_require_file "$SRC" "$REF" "$BWA" "$WGSIM" "$SAMTOOLS" "$GATK_SIF"
gatk_run() { singularity exec -B "${MTB_WORK}:${MTB_WORK}" "$GATK_SIF" gatk "$@"; }
mkdir -p "$OUTDIR"

# --- reads --------------------------------------------------------------------
# wgsim -e/-r/-R/-X at 0 gives an error-free, variant-free read set: the harness
# must not inject differences of its own, or the null controls in T0 are
# meaningless. Sequencing error is not what this program measures.
# SIM_FQ_PREFIX lets a caller supply the reads instead of simulating them. Two
# uses: comparing several references on BYTE-IDENTICAL reads (refbias_t5.sh, where
# re-simulating per arm would make a per-sample arm difference partly noise), and
# aligning a pre-selected subset of reads (arm D's extract-then-align). If the
# pair already exists it is reused verbatim; otherwise it is simulated to there.
FQP="${SIM_FQ_PREFIX:-${OUTDIR}/${PREFIX}}"
FQ1="${FQP}_1.fq"; FQ2="${FQP}_2.fq"
SRCPLAIN="${OUTDIR}/${PREFIX}.src.fa"
# SIM_REQUIRE_FQ=1 forbids the fallback. P1 and P2 pass an isolate's REAL reads
# this way, with $SRC set to the reference: if either FASTQ were missing, the
# fallback below would simulate reads FROM THE REFERENCE and call them as the
# isolate -- a perfect, silent substitution.
if [[ "${SIM_REQUIRE_FQ:-0}" == 1 && ! ( -s "$FQ1" && -s "$FQ2" ) ]]; then
    echo "FATAL: SIM_REQUIRE_FQ=1 but ${FQ1} / ${FQ2} are missing or empty;" \
         "refusing to simulate reads in place of real ones" >&2
    exit 1
fi
if [[ -e "$FQ1" && -e "$FQ2" ]]; then
    echo "[sim] ${PREFIX}: reusing reads ${FQ1} ($(( $(wc -l < "$FQ1") / 4 )) pairs)"
else
    case "$SRC" in *.gz) zcat "$SRC" > "$SRCPLAIN" ;; *) cp -f "$SRC" "$SRCPLAIN" ;; esac
    GLEN=$(awk '!/^>/{n+=length($0)} END{print n}' "$SRCPLAIN")
    NREADS=$(( GLEN * DEPTH / (2 * READLEN) ))
    echo "[sim] ${PREFIX}: genome ${GLEN} bp, depth ${DEPTH}x, ${NREADS} pairs of ${READLEN} bp"
    # SIM_ERR, SIM_FRAG and SIM_FRAG_SD exist for T9, which measures how
    # optimistic the error-free default is. The MUTATION rate -r and the indel
    # fraction -R stay at zero in every configuration: those inject differences
    # the source genome does not have, which would corrupt the truth set rather
    # than model the instrument. Only base error and fragment geometry vary.
    "$WGSIM" -e "${SIM_ERR:-0}" -r 0 -R 0 -X 0 -N "$NREADS" \
        -1 "$READLEN" -2 "$READLEN" -d "${SIM_FRAG:-$FRAG}" \
        -s "${SIM_FRAG_SD:-50}" -S 11 \
        "$SRCPLAIN" "$FQ1" "$FQ2" > /dev/null
    echo "[sim] ${PREFIX}: wgsim -e ${SIM_ERR:-0} -d ${SIM_FRAG:-$FRAG} -s ${SIM_FRAG_SD:-50}"
fi
# SIM_BAM_SUFFIX keeps one sample's arms from overwriting each other's BAM.
BAM="${OUTDIR}/${PREFIX}${SIM_BAM_SUFFIX:-}.bam"

# --- map ----------------------------------------------------------------------
# Reference artefacts are built ATOMICALLY: to a temp name, then mv into place.
# Two failure modes have already been hit here. A plain existence check silently
# reused a stale dictionary (228 sequences beside a 129-sequence reference), so
# these now compare mtimes. And six array tasks starting together on a
# freshly-downloaded reference all raced to create the dict; CreateSequenceDictionary
# refuses to overwrite, so five died with FileAlreadyExistsException and
# run_block.sh logged them as "item N failed". mv is atomic on a POSIX filesystem,
# so a loser simply overwrites with an identical file instead of failing.
if [[ ! -f "${REF}.bwt" || "$REF" -nt "${REF}.bwt" ]]; then
    _lock="${REF}.bwaindex.$$"
    cp -f "$REF" "$_lock"
    "$BWA" index "$_lock" >/dev/null 2>&1
    for ext in amb ann pac sa bwt; do   # bwt last: it is the "indexed" signal
        [[ -f "${_lock}.${ext}" ]] && mv -f "${_lock}.${ext}" "${REF}.${ext}"
    done
    rm -f "$_lock"
fi
DICT="${REF%.*}.dict"
if [[ ! -f "$DICT" || "$REF" -nt "$DICT" ]]; then
    _tmpd="${DICT}.$$"
    gatk_run CreateSequenceDictionary -R "$REF" -O "$_tmpd" && mv -f "$_tmpd" "$DICT"
    rm -f "$_tmpd"
fi
if [[ ! -f "${REF}.fai" || "$REF" -nt "${REF}.fai" ]]; then
    "$SAMTOOLS" faidx "$REF"
fi

"$BWA" mem -t "${MTB_THREADS:-8}" -R "@RG\tID:${PREFIX}\tSM:${PREFIX}\tPL:ILLUMINA" \
    "$REF" "$FQ1" "$FQ2" 2> "${OUTDIR}/${PREFIX}.bwa.log" \
  | "$SAMTOOLS" sort -@ 4 -o "$BAM" -
"$SAMTOOLS" index "$BAM"

# mapping summary: the fraction of reads that find a home is itself a bias readout
"$SAMTOOLS" flagstat "$BAM" > "${OUTDIR}/${PREFIX}${SIM_BAM_SUFFIX:-}.flagstat"

# --- call ---------------------------------------------------------------------
# SIM_NO_CALL stops here. Arm D of refbias_t5.sh wants the H37Rv alignment only to
# select reads from; calling on that backbone is T1's question, not its own, and
# HaplotypeCaller over 4.4 Mb is the slowest step in the script.
if [[ -n "${SIM_NO_CALL:-}" ]]; then
    rm -f "$SRCPLAIN"
    echo "[sim] ${PREFIX}: mapped to $(basename "$REF"), no variants called (SIM_NO_CALL)"
    exit 0
fi
# SIM_GVCF additionally emits a GVCF. The reference-composition test needs to
# know, at positions where the sample made NO call, whether that means "reads
# covered this and agreed with the reference" or "nothing was observed here" --
# a plain VCF cannot distinguish them, and conflating the two is what turns a
# missed call into an inherited wrong call.
if [[ -n "${SIM_GVCF:-}" ]]; then
    gatk_run --java-options "-Xmx8g" HaplotypeCaller \
        -R "$REF" -I "$BAM" -O "${OUTDIR}/${PREFIX}.g.vcf.gz" \
        -ploidy 1 -ERC GVCF --native-pair-hmm-threads "${MTB_THREADS:-8}" \
        > "${OUTDIR}/${PREFIX}.gatk.gvcf.log" 2>&1
fi
gatk_run --java-options "-Xmx8g" HaplotypeCaller \
    -R "$REF" -I "$BAM" -O "${OUTDIR}/${PREFIX}.vcf.gz" \
    -ploidy 1 --native-pair-hmm-threads "${MTB_THREADS:-8}" \
    > "${OUTDIR}/${PREFIX}.gatk.log" 2>&1

# Stamp the P0 graph-build provenance into the header. P0 names its asset
# directory after the graph's checksum so assets from different builds cannot
# mix, but a VCF produced against those assets does not otherwise say which build
# made it -- and since the graph is rebuilt when bugs are found, output from a
# superseded build would be indistinguishable from output from the fixed one.
# That distinction is unrecoverable after the fact, so it is stamped at
# production time. The stamper is a no-op when no P0 build directory exists, so
# harnesses predating P0 are unaffected; when a build IS resolved and stamping
# fails it exits non-zero and fails the sample, because a provenance step that
# fails quietly is worse than one that is absent.
_stamp=( "${OUTDIR}/${PREFIX}.vcf.gz" )
[[ -n "${SIM_GVCF:-}" && -s "${OUTDIR}/${PREFIX}.g.vcf.gz" ]] && \
    _stamp+=( "${OUTDIR}/${PREFIX}.g.vcf.gz" )
bash "$(dirname "${BASH_SOURCE[0]}")/stamp_build_id.sh" "${_stamp[@]}"

[[ -n "${SIM_KEEP_FQ:-}" ]] || rm -f "$FQ1" "$FQ2"
rm -f "$SRCPLAIN"
echo "[sim] ${PREFIX}: $(zcat "${OUTDIR}/${PREFIX}.vcf.gz" | grep -vc '^#') raw variant records"
