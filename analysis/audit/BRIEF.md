# Pipeline audit brief (2026-10-05)

You are auditing part of the mtbc-graph-variants pipeline for bugs before a new
graph build and a full genotyping rerun. The user wants it gone through "with a
fine tooth comb", because two bugs (below) were found late, in code that had
not been reviewed.

## Locations

- **Repository (reviewed code):** `/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants`.
  Read `CODE_REVIEW.md` (the earlier review; items marked fixed should be
  re-verified, not assumed) and `HANDOFF.md` sections 0b, 0j and 2.
- **Working tree (production data, plus some code never brought into the
  repository):** `/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome`. The
  association chain runs there: `assoc/bin/*`, `bin/vcf_to_alignment.py`,
  `bin/build_snp_tree.sh`, `bin/pggb_build.sh`, `bin/vcf_decompose.sh`,
  `bin/vg_deconstruct.sh`, `bin/build_panel.py`, `bin/snp_nonredundant.py`,
  `sv2frame/bin/sv_twoframe.sh`, `is6110/bin/panisa_score.py`,
  `bin/run_snpEff.sh`, `bin/assemble_from_cram.sh`. Where a script exists in
  both places, find out which copy production actually runs.
- **Production runs:** `bin/refbias_run.sh` from
  `mtbc-graph-variants/runroot`, which links the repository's `bin/`,
  `config/`, `accessory/bin`, `is6110/bin`, `graphframe/bin`, `sv2frame/bin`
  and the working tree's data. The association chain is
  `assoc/bin/cohort_assoc_tail.sh`, run from the working tree root.
- **Production outputs to check claims against:** `refbias/scale200/`,
  `refbias/gwas1000/`, `assoc/scale200/`, `assoc/gwas1000/`, `data/trees/`,
  `refbias/build/7713a8d71d8e/`, and the graph
  `graphs/CX333.s10k.k23.K15/`.

## Rules (strict)

1. **Read-only.** Do not edit any code or data, in either location.
2. **No cluster jobs.** No `sbatch` or `srun`. Small local tests are fine (a few
   minutes of CPU).
3. **Write only:**
   - your report, `analysis/audit/<area>.md`, in the repository;
   - throwaway files in your own scratch directory,
     `/tmp/claude-12043/audit_<area>/`.
4. **Never write into the working tree.**
5. **Report; do not fix.** Each finding gives a suggested fix in words or as a
   short diff in the report.

## Bug patterns already found here (look for more of the same kind)

1. **Decomposed graph VCF: one event, several records.**
   `all_variants.decomposed.vcf.gz` splits an event into one record per allele
   path. The same (pos, ref, alt) can occur several times, with a sample's
   genotype 1 in only one of them. Nested snarls (LV > 0) hold most SNPs.
   - Fault A: `add_outgroup.py` assigned `exact = g` in a loop, so the last
     duplicate won (about 500 wrong outgroup cells per cohort).
   - Earlier misreads: reading top-level (LV=0) records only; pooling by snarl.
   - **Check every reader of a graph VCF** for: dict overwrite; first or last
     record wins; LV filtering; "no record means REF".
2. **Tree node chosen by position, not by content.** Fault B:
   `ancestral_alleles.py` takes "the root's child that is not the outgroup" as
   the MTBC ancestor. In the panel tree that node also contains the second
   M. canettii (GCF_000253375). Look for any code that picks nodes structurally
   (root children, first or largest child) where the meaning depends on which
   taxa are below.
3. **Silent empty outputs or defaults.** Examples:
   - a missing `--lineages` gave an all-blank `p_lineage` column;
   - missing presence tables gave blank `cond_*` columns.

   Look for optional inputs whose absence silently changes results, and for
   `.get(x, default)` or `try/except: pass` that hides missing data.
4. **Steps run by hand and never added to the chain.** The small-variant and SV
   gene burdens were missing from the chain. Look for products that docs or
   HANDOFF reference but no chain script produces.
5. **Skip-if-exists guards on the wrong file.** A partial output was treated as
   finished. Check every "already done" test.
6. **Coordinate and frame errors:**
   - 0-based against 1-based positions (BED, PAF, VCF, numpy indexing);
   - H37Rv frame against reference frame against node frame;
   - the strand of reverse alignments;
   - PanSN names (`GCF_x#1#contig`) against bare contig names.
7. **Allele and genotype handling:**
   - multi-allelic records;
   - `*` alleles;
   - missing genotypes counted as REF;
   - haploid against diploid GT parsing;
   - phased `|` separators.
8. **Ordering and identity:**
   - sample order mismatches between files;
   - joins on keys that are not unique;
   - silent drops on key mismatch.
9. **Docs against code:** documented settings that the code does not use (for
   example, RUNBOOK.md and QC_PIPELINE.md disagreed on the foreign-screen
   background).

## Report format (`analysis/audit/<area>.md`)

**One section per finding:**

- an ID (`<AREA>-<n>`) and a severity:
  - **HIGH:** changes production results materially;
  - **MEDIUM:** changes results in a minority of records, or would at scale;
  - **LOW:** robustness or clarity;
- `file:line`, and which copy runs in production;
- **what is wrong,** and the evidence: the code excerpt, plus whenever possible
  a check against a real output showing the effect (counts);
- **effect on current outputs:** which cohorts and products, and roughly how
  many records;
- **suggested fix.**

**Then:**

- a short "checked and found sound" list, so the coverage is visible;
- the open questions you could not settle.

**Be concrete and verify before reporting.** A false alarm costs the user
time. Mark anything you could not verify as UNVERIFIED.
