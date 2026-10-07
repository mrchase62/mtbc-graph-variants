# Pipeline audit: consolidated findings

2026-10-05. Seven read-only audits, one per area; each report is in this
folder. No code or data was changed, and no jobs were submitted.

**Totals: 97 findings: 10 HIGH, 39 MEDIUM, 48 LOW.**

| area (report) | HIGH | MEDIUM | LOW |
|---|---:|---:|---:|
| graph-VCF readers (`graphvcf.md`) | 1 | 4 | 2 |
| tree and polarity (`tree_polarity.md`) | 0 | 6 | 5 |
| association tests (`assoc_tests.md`) | 1 | 5 | 8 |
| P0-P2 and orchestration (`p0_p2.md`) | 1 | 3 | 10 |
| IS6110 and accessory (`p3_is6110_accessory.md`) | 1 | 4 | 5 |
| P4-P5 and SV genotyping (`p4_p5.md`) | 1 | 7 | 9 |
| panel and graph build (`panel_graph_build.md`) | 5 | 10 | 9 |

**"Verified" below** means I re-checked the finding myself, in code and, where
marked, on production data. The rest rest on the auditor's own counts from
re-running the logic on production files. Every auditor's original-logic run
reproduced the production outputs, which supports those counts.

## A. Findings that change current scale200 and gwas1000 results

These make the current association results provisional, chiefly the IS6110
results, the branch-null survivors and the SV results.

| ID | what | size | verified |
|---|---|---|---|
| **P3IS-1** (HIGH) | IS6110 carriers whose matched reference already has a copy at the site are written REF (GT=0) on H37Rv-frame records. Where H37Rv has no copy they look like non-carriers; where H37Rv has one, polarity is inverted. | scale200: 1,098 wrong cells against 469 correct ALT; gwas1000: 5,103 against 2,424 | yes, on data (1,098 exact) |
| **ASSOC-1** (HIGH) | A read-based *M. canettii* isolate (`canettii`) is in both cohorts with phenotype 0, as a tip of both association trees. It inflates the length-weighted branch null. | its branch is 25.8% (scale200) and 11.4% (gwas1000) of tree length; removing it cuts q_branch passes from 138 to 31 and 68 to 38 | yes, on data |
| **PGB-6** | A third artifact genome: GCF_039770655 (the only lineage 9 genome, the reference for every lineage 9 isolate) has a 377 bp poly-T. It produces false insertion and deletion calls aligned perfectly with lineage 9. | all 16 lineage 9 isolates | yes, on data |
| **P4P5-4** | delly records genotyped 0/0 but FILTER PASS become ALT deletion carriers | 883 rows (scale200), 6,604 (gwas1000) | auditor |
| **P4P5-3** | one clip cluster "confirms" an inherited deletion even when H37Rv-frame depth is full | 995 SV ALT cells (scale200), 7,321 (gwas1000) | auditor |
| **P4P5-1** | any non-exact odgi projection is written ABSENT (GT=2, `*`) | about 12,000 SNP cells (scale200), 109,000 (gwas1000). scale200 SNP records hold 134,308 ABSENT cells in all; how many are real deletions is to be settled. | code yes; size partly |
| **P4P5-2** | core keys stated REF where the matched reference carries the ALT and the sample matches its reference | 1,089 scale200 cells (29% of cells at reference-carried core indels of 50 bp or more) | auditor |
| **P3IS-2** | 674 of 802 accessory-presence loci are sequence H37Rv already has, so the read route cannot see them, but they are still written absent | 481 IS6110-element loci have 0 PRESENT calls; about 575,000 absent cells in gwas1000 | auditor |
| **P3IS-3** | node-frame IS6110 keys on repeated graph nodes merge unrelated insertions (over 1 Mb apart) | gwas1000: 1,861 of 3,585 node-frame rows at risk | auditor |
| **P3IS-4** | DR-array rescue offset error | 1 of 24 gwas1000 rescues | auditor |
| **Fault A / GRAPHVCF-1** | `add_outgroup.py`: the last duplicate record wins | 386 (scale200), 323 (gwas1000) | yes (check 1) |
| **TP-1 / GRAPHVCF-2** | `add_outgroup.py` writes REF where it should write N (all duplicates missing; inside outgroup deletions) or ALT (upstream MNP or indel records ignored) | about 300-500 cells per cohort per cause; events change for 407 more variants | auditor |
| **TP-2** | `vcf_to_alignment.py` drops SNP records written `X,*`: no tree column, no outgroup, no polarity lookup | 11,163 scale200 SNP records (19%); 2,592 left with an unresolved root | yes, on data |
| **Fault B** | the ancestral-allele node includes the second canettii | 132 lineage 1-4 variable sites | yes (check 2) |
| **GRAPHVCF-3 / TP-8** | `panel_polarity.py` (working tree, hand-run) has the last-duplicate fault | 903 keys; 39-42 variants per cohort change polarity; none significant | yes, in code |
| **GRAPHVCF-4** | P1 reference selection keys the panel by position only, ignoring the ALT base; 319 padded SNPs are dropped | 10 isolates on a slightly worse reference (1-12 SNPs) | auditor |
| **GRAPHVCF-5** | the decompose round-trip leaves REF padded on 1,399 records (376 SNPs), which lose their ancestral allele | 112 SNP variants per cohort | auditor |
| **TP-3** | the H37Rv tip is set to REF at node-frame records, where H37Rv has no sequence | 62 / 423 variants gain 2 or more events only because of it | auditor |
| **ASSOC-2 to ASSOC-5** | the region-null leave-one-out never removes anything; deduplicated-record lookup; SV interval tiers not passed (and stale); the SV burden credits only the gene at the anchor base (21% of deletions miss genes) | see the report | auditor |
| **P4P5-5, -6, -8** | catalogue intervals widen to the union; one insertion emitted twice; inherited records overlapping the sample's own call | see the report | auditor |
| **PGB-8** | 79 of 333 panel genomes never had the SNP-outlier screen | panel composition | auditor |

## B. Hazards for the rebuild and rerun (no effect on current outputs)

- **Skip guards that only test whether a file exists:** P1 and P2
  (p0_p2 HIGH), P0 ancestral and the association chain (TP-4), accessory tables
  (P3IS-8). A rerun into the existing folders would keep stale results.
  **Rerun into new output folders, and key the guards on build ID.**
- **Paths hard-coded to the CX333 graph:**
  - panel SNPs (`p0_prepare.sh:46-48`);
  - the node tables and panel allele frequencies (P4P5-7, HIGH);
  - the build ID stamp (PGB-14);
  - the IS6110 tie-break interval files (p0_p2).
- **Repeat-mask provenance:** the mask changed on Sep 28 (181 → 395 intervals)
  under a build whose manifest records the old checksum.
- **Graph build:**
  - `pggb_build.sh` cannot start clean, because it writes provenance before
    the emptiness check (PGB-1, verified);
  - its defaults differ from CX333's (PGB-2);
  - vcfwave settings differ from production (PGB-10).
- **The panel build is not fully in code:**
  - the 484 → 333 reduction (PGB-3);
  - the reference-guided rule (PGB-4);
  - the foreign screen (PGB-5);
  - the RefSeq regex drops 2 M. orygis genomes (PGB-13);
  - P0 references are not the graph's rotated sequences (PGB-15).
- **Duplicate records:** they come from the decompose step (vcfwave). The
  collapse step already removes them almost losslessly (64 cells). **Readers
  of the new graph should use the collapsed VCF,** and the sync should keep
  it.

## C. Code location

- **Production runs code in two places.** The association chain, the tree
  builders, `panel_polarity.py` and all graph and panel build code exist only
  in the working tree. None of it was in the earlier review.
- **Stale duplicates:** the working tree holds old copies of
  `cohort_assoc_tail.sh` and `audit_chain.py` (ASSOC-6), and an old
  `sv_twoframe.py` (P4P5-16). Running the old chain script brings back the
  three omissions fixed in HANDOFF 0b.
- **Fix:** bring every production script into the repository, run only from
  the repository (through runroot), and remove the ambiguity.

## D. LOW findings (48)

These are robustness items, documentation mismatches, missing tests and minor
counts. They are listed per report.

## Proposed order of work

1. **Bring all production code into the repository** (section C), one commit
   with no behaviour change.
2. **Fix section A, with a test for each finding.**
   - Each test reproduces the bug on a small case and on real data copies, in
     scratch or `analysis/`.
   - The cohort-design item ASSOC-1 needs your decision first.
3. **Fix section B, and pass every path through the build** rather than
   hard-coding it.
4. **A second review of the fixes** (code review of the diff), then a check of
   the LOW items.
5. **Choose and build the new panel and graph:** the artifact genomes go,
   section B's build fixes are in, and the test-graph results set the
   composition.
6. **Rerun into new folders:** scale200 first, checked in full against these
   findings (each fixed count should move as predicted), then gwas1000. Costs
   go to you before each submission.
