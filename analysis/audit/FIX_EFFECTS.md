# How the audit fixes change the results

Last updated 2026-10-07. Fixes are on branch `audit-fixes`.

**What this document is.** Each fix is compared with the current production
outputs *at the stage the fix changes*: P1 reference choice, P4 variant keys,
and association burden units. These are measurements, not projections.

**What it is not.** It does not yet show the end results: the P5 matrix,
trees, association survivors and drug-resistance positive controls. Those
come from one scale200 rerun after all 8 queued fixes are in (cost estimate
first). That comparison will be added here as Part 2. In that rerun the fixes
act together, so a change in an end result can be pinned on one fix only by
rerunning with that fix off.

## Summary

| Fix | Stage | What changes | scale200 | gwas1000 |
|---|---|---|---|---|
| D41: node alleles on the forward strand, plus node-indel left-alignment | P4 keys | one event read in opposite directions got two keys | 11 bases with two keys → 0; homopolymer events with two keys 3 → 0; 73 records now keyed in H37Rv | not run (P4 not rerun for gwas1000) |
| Fix 1 (D18, D39): genes credited by the bases a variant changes | burden units, genic stratum | SNPs in overlapping genes, MNPs and insertions at gene ends | 271 of 63,988 records change unit (0.4%) | 584 of 150,118 (0.4%), 9 of them IS6110 |
| Fix 2 (D20, D24): reference picked by mismatches per site both sides called | P1 reference | uncovered sites counted as matching; raw counts not normalised | 24 of 200 isolates get a different reference (12%) | 112 of 997 (11%) |
| Fix 3 (D21): H37Rv as a candidate reference, and a path projected onto itself is the identity | P1 reference; P3-P5 and IS6110 projections | H37Rv could never be chosen; odgi mapped a path onto another copy of itself in repeats | 1 of 200 isolates now maps to H37Rv | 0 of 997; across the build's IS6110 projection store, 22 of 2,707 same-reference lookups were on the wrong copy |
| Fix 4 (D22): graph VCF left-aligned, then collapsed | graph VCF (panel_af, inherited half) | indels in repeats at vcfwave's position, not the cohort's | build asset: 8,304 of 93,214 keys move; 49 records fold into 41 keys | same asset |
| Fix 5 (D29): UniVec vector check in the foreign screen | panel QC (which genomes are reviewed) | a construct shared by 2+ genomes could pass as native | — | — (panel: of 2,703 inserts, exactly the 2 known constructs flagged VECTOR) |
| Fix 6 (D32): foreign-screen background chosen by quality | panel QC (the screen's background) | each sublineage's background genome was the first accession, not the best | — | — (CX333: 44 of 65 slots have a choice; the effect depends on the new panel's score file) |

## D41: one event, one key (scale200, P4 rerun on all 199 samples)

**Problem.** Variants on sequence that is off the H37Rv path are keyed as
`node:<id>:<offset>:REF>ALT`. Each sample's reference wrote them in the
direction that reference walks the node. Two references walking a node in
opposite directions gave one event two keys, e.g. `C>T` and `G>A`. The
cohort matrix then split one variant into two rows, each carried by part of
the cohort.

**Measured effect:**

- Same-base events with two keys: 11 → 0.
- Records restated or re-keyed: 22 on the reverse-strand group and about
  116 on the along-strand group.
  - 73 are now keyed in H37Rv coordinates, so they join the H37Rv rows other
    samples already wrote. Where another sample carries the event, the keys
    agree.
  - 30 core records are left to the direct arm.
  - 45 stay as the reference reads them: inversion junctions, nodes visited
    twice, and anchors off the path.
- Node indels are now left-aligned on the node's own sequence. 59 records
  moved. Keys that were not left-aligned: 42 → 0. Homopolymer events with two
  keys: 3 → 0 (example: node 107910, `G>GAG` at 3 and `T>TGA` at 2 were one
  insertion).

**Expected downstream effect: small.** A handful of matrix rows merge, so
those variants have higher carrier counts. Node-frame rows are a small share
of the matrix.

## Fix 1: the bases a variant changes decide its genes

**Problem.**

- A SNP in two overlapping genes counted only for the gene that starts
  first.
- An insertion was assigned by its anchor base, so an insertion just past a
  gene's last base counted inside that gene.

The gene burdens and the scan's genic/intergenic split both used the anchor
rule.

**Measured effect**, on production event tables (records the burden admits;
no callability floor applied):

| | scale200 | gwas1000 |
|---|---|---|
| records checked (small + SV + IS6110) | 63,988 | 150,118 |
| gain an overlapping gene | 260 | 553 small + 7 IS6110 |
| move to a different single unit | 10 | 10 small + 2 IS6110 |
| other change | 1 | 3 |
| SVs changed | 0 | 0 |

Examples:

- `h37rv:95415 T>C`: hycQ → hycQ **and** hycE (overlapping genes).
- `h37rv:60400 C>T`: Rv0057 → Rv0057 **and** dnaB.
- `h37rv:333136 T>TTGCCGTTCTG`: was in vapC25; the insertion follows the
  gene's last base, so it now goes to vapC25's promoter (`up:vapC25`).
- IS6110 at 850040: was in PPE12; it lies past the gene, so now
  `up:PPE12`.
- IS6110 at 3494459: was intergenic (Rv3128c-Rv3129); now Rv3129's
  promoter.

**Expected downstream effect:**

- Burden tests for overlapping gene pairs gain carriers, e.g. toxin-antitoxin
  pairs (vapB1/vapC1) and operons such as dnaB/Rv0057.
- A few promoter units gain or lose IS6110 insertions.
- About 0.4% of burden records change, so survivor lists should be mostly
  stable. Any change will be in overlapping-gene units.

## Fix 2: the matched reference, by mismatches per called site

**Problem.**

- The reference with the fewest SNPs to the isolate was chosen, but sites
  the isolate had no reads at counted as matching H37Rv.
- Counts were not normalised by how many sites were compared, so a reference
  with fewer comparable sites looked closer.

**New rule:**

- Compare only panel sites both sides called. The isolate side comes from
  its P1 H37Rv BAM (≥5 reads, MAPQ ≥20).
- Mask filtered and indel/complex isolate records.
- Rank by mismatches per compared site.
- References compared over fewer than half the best site count are not
  eligible.

**Measured effect** (new selector run on every isolate):

| | scale200 | gwas1000 |
|---|---|---|
| isolates | 200 | 997 |
| different reference chosen | **24 (12%)** | **112 (11%)** |
| new pick was the old 2nd choice | 18 | 87 |
| old pick now worse by ≤2 SNPs | 13 | 59 |
| … by ≤5 SNPs | 21 | 91 |
| … by >10 SNPs | 1 | 5 |
| panel sites compared per isolate (median, of 75,587) | 74,145 | 74,140 |

Examples (scale200):

- SAMEA114663528: GCF_932527465 → GCF_977011305; old pick 84 vs new 81
  mismatches (a near-tie).
- SAMEA1100846: GCF_014901095 → GCF_977011305; 209 vs 208.
- SAMEA1403838: GCF_014899745 → GCF_000193185; 74 vs 74 per about 74,000
  sites, decided by the rate.

**Expected downstream effect:**

- Most changes swap one close relative for another. That changes which
  reference-frame calls an isolate contributes but should change few
  H37Rv-frame matrix cells.
- The few isolates whose old reference was more than 10 SNPs worse should
  see the clearest gain: fewer reference-bias calls near their own
  diversity.
- This is the fix most likely to move end results, because it changes the
  whole P2-P4 input for about 11% of isolates.

## Fix 3: H37Rv as a candidate reference (D21)

**Problem.** H37Rv is the reference the graph was decomposed against, so it
had no column in the panel VCF and could never be chosen as a matched
reference, even for an isolate closer to H37Rv than to any other panel
genome.

**Measured effect on reference choice** (same selector run as fix 2, with
H37Rv added):

- scale200: 1 isolate changes. SAMEA7526648 moves from GCF_013267635 to
  H37Rv: 29 vs 36 SNPs over about 74,500 sites, a clear win rather than a
  near-tie.
- gwas1000: none. H37Rv is not in any isolate's top two.
- No other isolate's choice changes.

**Found while testing R = H37Rv through P2-P5.** Projecting a position onto
its *own* path through the graph is not always the identity. Where a path
passes the same graph nodes more than once (tandem repeats), odgi answers
with one of the copies.

- Today's P4 on the earlier H37Rv-pinned scale200 arm (100 isolates): 729 of
  51,139 composed records (1.4%) were keyed away from their own position,
  by up to 1 kb, all in PE/PPE repeats (552 at 3.93-3.95 Mb). After the
  fix: 0.
- The P5 direction (H37Rv keys back onto R = H37Rv): 16 of 52,007 positions,
  by a median of 4 bp.
- **This also affects current results, for any reference.** In the IS6110
  arm's stage 2, a key whose carrier shares the sample's reference is looked
  up on that reference itself. 22 of 2,707 such lookups in the current
  projection store landed on another copy, by up to 69 bp, so those IS6110
  cells were read at the wrong place.

**Expected downstream effect:** one isolate in scale200 is called against
H37Rv instead of its old reference. IS6110 genotypes change in a handful of
repeat-region cells.

## Fix 4: the graph VCF, left-aligned (D22)

**Problem.** The graph VCF kept vcfwave's indel positions, trimmed but not
left-aligned, while every cohort key is left-aligned. An indel in a repeat
was therefore one event under two positions:

- P5's `panel_af` lookup came back blank for it;
- vcfwave sometimes wrote one event at two positions, as two records with
  the carriers split between them.

**Measured effect** (today's collapse on the production decomposed VCF,
with and without left-alignment):

- 8,304 of 93,214 keys move to their left-aligned position;
- 49 records fold into 41 keys, with their carriers unioned;
- no carrier is lost, and the result is identical to left-aligning
  afterwards.

**Also found:** the production graph VCF predates the earlier collapse fix
(GRAPHVCF-5, already committed). Against today's code it has 1,461 padded
keys (e.g. `CT>GT` for the SNP `C>G`) and 61 missing carrier cells. The new
build regenerates it.

**Expected downstream effect:**

- `panel_af` filled for indels in repeats;
- the panel and cohort agree on indel positions;
- the 41 merged events show their full carrier counts.

## Fix 5: a vector check in the foreign screen (D29)

**Problem.** The foreign screen calls an insert native when other genomes
carry it. A construct carried by two genomes passes as REVIEW (each vouches
for the other), and one carried by three would pass as native.

**New rule:** every insert, native or not, is searched against NCBI's
UniVec_Core with VecScreen's blastn settings. A strong match makes it
`VECTOR`.

**Measured effect:**

- Of CX333's 449 inserts, exactly two get strong matches: the two known
  constructs, GCF_044324775 (pJEB, 2,295 bp of vector) and GCF_021535155
  (attB vector, 2,107 bp).
- None of the external assemblies' 2,254 inserts matches, not even
  moderately.
- End to end on 7 genomes, those two move from `REVIEW_one_homologue` to
  `VECTOR`; every other verdict is unchanged.

**Expected downstream effect:** none on the current panel, where both
genomes were already excluded by hand. For the new panel build, a construct
can no longer pass by being common.

## Fix 6: the foreign screen's background, chosen by quality (D32)

**Problem.** The foreign screen compares each insert with a small
background: one genome per sublineage. The genome taken for each sublineage
was simply the first accession in sorted order, so the background's quality
was an accident of numbering.

**New rule:** take the best genome per sublineage by `--rank-by`, the same
score file and rule that pick a clone cluster's representative (D27). A
competing genome with no score stops the run.

**Measured effect:** the rule matters in 44 of CX333's 65 sublineage slots,
where 305 genomes compete. Which genomes it picks depends on the score file,
which will be defined with the new panel's clone-collapse step. Nothing
changes in current results: the screen is panel QC and runs again with the
new panel.

## Still to come

**Fixes 7-8**, each to be added here as it is done:

- D38: refuse outputs from before the guards.
- D43: a missing catalogue path stops with an error.

**Part 2 (after the scale200 rerun):** old vs new for:

- reference choice;
- P4 keys;
- P5 matrix cells by state;
- tree topology (Robinson-Foulds);
- association survivors and q-values;
- the drug-resistance positive controls (katG, rpoB, embB, gyrA, the fabG1
  promoter and others).
