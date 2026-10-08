# scale200 rerun on the fixed code (build 7713a8d71d8e-fix1)

Approved by the user on 2026-10-08 ("yes, run the scale200 rerun at 380
billing-hours"). Purpose: data to evaluate the audit fixes. scale200's 200
isolates are run again on the CX333 graph, with every fix on main, into new
folders, and then compared with the old run (`refbias/scale200`, finished
2026-10-01/02).

## What is new and what is the same

| | old run | this run |
|---|---|---|
| graph | CX333 (`graphs/CX333.s10k.k23.K15`, sha256 7713a8d7...) | the same graph |
| build | `refbias/build/7713a8d71d8e` | `refbias/build/7713a8d71d8e-fix1` |
| cohort (registry) | `scale200` | `scale200_fix` (same isolate table and CRAM table) |
| outputs | `refbias/scale200`, `refbias/work/scale200_*`, `accessory/scale200`, `assoc/scale200` | `refbias/scale200_fix`, `refbias/work/scale200fix_*`, `accessory/scale200_fix`, `assoc/scale200_fix` |
| code | main as of 2026-10-01/02 | main as of 2026-10-08 (all audit fixes) |
| phenotype | `assoc/scale200/rrdr_carriers.txt` (made 2026-09-25) | the same file, so a change in results comes from the fixes |
| canettii isolate in the cohort tree | kept | kept (ASSOC-1 is still the user's decision) |

The build gets its own id because P0 would otherwise name it after the
graph's checksum, which is the old build's id. Every pass is pointed at it
with `MTB_BUILD_DIR`.

## Build inputs made outside P0

All in `$MTB_WORK/refbias/build_inputs/7713a8d71d8e-fix1/`. The graph folder
and the old build are only read.

| step | script | product |
|---|---|---|
| 1 graph VCF | `01_graph_vcf.sbatch` | `graphvcf/all_variants.collapsed.vcf.gz` (the graph's 2026-09-11 decomposed VCF through the fixed `bin/vcf_collapse.sh`: GRAPHVCF-5 trim-and-union, D22 left-alignment), and `snps/indels/svs/small_variants.vcf.gz` by `bin/vcf_split_classes.sh --outdir` |
| 2 accessory panel | `03_accessory_panel.sbatch` and the blast array | `accessory_panel/` by `docs/PANEL_TREE.md` section 4 |
| 3 panel tree | `04_panel_tree.sbatch` | `tree/cx333.{snps.fasta,sites.tsv,rooted.nwk}` by `docs/PANEL_TREE.md` section 1 |
| 4 IS6110 check | `05_is6110_check.sbatch` | `is6110_check/`: both IS6110 stages rerun into a separate folder and compared file by file with the shared `is6110/assets/` (which the chain reads) |

P0 steps run through `02_p0.sbatch` from `runroot`. Not regenerated:

- the 333 NCBI annotation files: copied from the old build (identical bytes;
  no download);
- the IS6110 crossmaps: the chain reads the shared copies, which step 4
  checks reproduce, and P0's `is6110_intervals` checks against this build's
  references.

The decomposed VCF is reused because the decomposition settings in
`bin/vcf_decompose.sh` (vcfbub `-l 0 -a 100000`, vcfwave `-I 1000`) are the
ones that made it; only the collapse after it changed.

## Then

1. `bash bin/refbias_run.sh scale200_fix --dry-run`, then the full chain
   p1 to p5vcf, from `runroot` with `MTB_BUILD_DIR=refbias/build/7713a8d71d8e-fix1`.
2. The association tail (`assoc/bin/cohort_assoc_tail.sh scale200_fix
   assoc/scale200_fix/rrdr_carriers.txt`).
3. The old vs new comparison, report only.

## Findings so far

- The fixed collapse gives 79,474 SNP/MNP records, 9,027 indels and 4,652
  SVs over 332 genomes, with no duplicate key and every indel left-aligned.
- The accessory candidates are the same 1,643 insertions, but D22 moves most
  of them left: 108 keep their position, most move 1 to 14 bp. Accessory
  locus ids carry the position, so the comparison matches loci by sequence
  and carriers, not by id.
