# scale200: key outputs, durable copy

Copied 2026-10-08 from netscratch (90-day purge) by
`analysis/rerun_scale200/archive_to_mirror.sh` in the repository
`mtbc-graph-variants`. BAMs and per-sample working files are not here; they
remain on netscratch until purged.

**The cohort:** 200 *M. tuberculosis* isolates (100 from scale100 plus 100 more)
and one *M. canettii* isolate, run against the CX333 pangenome graph.
- Each isolate is aligned to its closest panel genome (the matched
  reference), and calls are reported in H37Rv coordinates.
- The phenotype for the association is RRDR carriage: 57 of 202 tips.

| folder | what | pipeline source |
|---|---|---|
| `2026-10-01_original/` | the run on the code of 2026-10-01/02, build `7713a8d71d8e` | `refbias/scale200`, `assoc/scale200`, `accessory/scale200` |
| `2026-10-08_audit_fixes/` | the rerun on all audit fixes, build `7713a8d71d8e-fix1` | `refbias/scale200_fix`, `assoc/scale200_fix`, `accessory/scale200_fix` |
| `COMPARISON.md` | old vs new report | `analysis/rerun_scale200/COMPARISON.md` |

The pipeline's working folders are still named `refbias/` for historical
reasons: the project began as a reference-bias study.

**Inside each run:**

| folder | files | what they are |
|---|---|---|
| `calls/` | `merged.vcf.gz` (+ `.tbi`) | every call in every sample. FORMAT/ST holds the pipeline state: ALT, REF, ABSENT (region missing in the sample) or NOCALL |
| | `sites.tsv`, `keys.tsv` | site and key tables behind the VCF |
| | `validation.tsv`, `sanity.tsv`, `states.meta.tsv` | P5 checks |
| | `build_stamp.txt` | the graph build the run used |
| `structural_variants/` | `sv_intervals.tsv` | the SV interval catalogue |
| | `sv_matrix.tsv` | sample x interval states |
| `reference_choice/` | `refmap.tsv` | the matched reference P1 chose for each isolate, and its SNP distance |
| | `p2_summary.tsv` | per-sample P2 call counts |
| `association/` | `scan.tsv` | the variant-level association scan: three nulls (branch, region, lineage); survivors pass all three at q < 0.05 |
| | `small_gene.tsv`, `sv_gene.tsv`, `is6110_gene.tsv` | gene-level burden tests |
| | `events/` | the per-branch event matrix |
| | `chain_audit.txt` | the end-to-end check |
| | `rrdr_carriers.txt` | the phenotype |
| `accessory_presence/` | one table per isolate | the presence state of each accessory locus |
| `trees/` | `*.rooted.nwk` | the cohort tree |
| | `*.combined.*` | the cohort + CX333 panel tree (IQ-TREE outputs) and the SNP alignments |

In the original run, the cohort + panel tree files are named `scale200_cx333.*`.
