# RD inheritance and SV ancestral reconstruction, new 1K cohort (2026-10-16)

Analysis only. VCF `runroot/refbias/new1k_p200/p5/merged.vcf.gz` (1,000 L1-4 isolates, P200 graph
3893a3762ed1), read through `analysis/sv_eval_1k/out/new.pkl` (that directory's `extract.py`: every
non-small record plus CLASS=small indels >= 50 bp, per-cell ST and EV). Nothing in the pipeline was
changed.

## Methods, in standard terms

- **Tree and root.** The cohort tree is `analysis/sv_eval_1k/out/new1k.fasttree.nwk` (FastTree GTR+gamma
  on the VCF's biallelic SNPs outside the repeat mask). L1, L2, L3, L4 and L2+L3 are each monophyletic
  in it (0 misplaced isolates). **Root: on the branch that separates L1 (250) from L2-4 (750)**; the
  cohort has no L5/6 outgroup, and L1 is the first split within L1-4. One consequence: a character that
  splits exactly L1 vs L2-4 ties at the root (`TIED`). Only the panel's L5/L6 outgroups can resolve
  that, so the panel AA_STATE is the reference for those.
- **Fitch parsimony** (`common.py`, vectorised): leaves ALT=1, REF=0, and ABSENT/NOCALL/DUP unknown.
  The up pass gives the minimum change count, using the k - max rule at multifurcations (the same rule as
  `analysis/sv_homoplasy/count.py`). The down pass starts from the root state (REF when tied), keeps the
  parent's state wherever the child's set allows it, and counts changes away from the root state
  ("derived") and back to it ("reversals"). A record is scored when it has >= 2 ALT and >= 2 REF isolates.
- **RD matching** (`rd_lineage.py`): `known_RDs.tiered.tsv` (135 RDs) against H37Rv-frame DEL records
  (CLASS sv/is6110, any COMPONENT, including small_allele). A match needs reciprocal overlap >= 0.5 of the
  deleted bases. The primary record is the match with the highest overlap; with no match, the best
  overlapping record is shown, flagged `n_match=0`. Expected carriers are isolates whose TB-Profiler
  deepest sublineage is the RD label or below it (`(except x)` labels honoured). Isolates typed only
  above the RD's level are "unresolved". Labels that name no L1-4 group (L5/6/8/9, animal, canettii)
  expect no carriers.
- **Near-identical pairs** (`pairs.py`, one Slurm job of ~2 min on 4 CPUs): the SNP distance comes from
  `snps.fasta` (89,791 columns, about 89,200 compared per pair, N ignored). For every pair at <= 5 or
  6-12 SNPs, and every SV-like record: a cell pair is "both called" when both cells are REF/ALT, and
  "discordant" when one is ALT and the other REF. Cells are split by class, COMPONENT, repeat mask, EV
  and whether the two isolates share the same P1 matched reference. `odd_side.py` finds which isolate
  is the odd one out, using the record's ALT frequency in the rest of the pair's sublineage (>= 0.8
  means the REF side is odd; <= 0.2 means the ALT side is odd). `private_alt.py` counts, per isolate,
  the DEL-interval ALT calls its own sublineage does not share.
- **Ancestral reconstruction** (`ancestry.py`): DEL and INS (CLASS sv, any COMPONENT), plus ACC and
  IS6110 for context, with 3,566 records scored. Gene annotation comes from the H37Rv snpEff dump; IS6110
  distance is measured to the H37Rv copies (`is6110_elements.gff`) and to the cohort's IS6110 insertion
  sites. `locus.py start end` prints per-lineage states of the records in a window.

## 1. Known RDs by lineage (`out/rd_records.tsv`, `rd_by_lineage.tsv`, `rd_by_sublineage.tsv`, `rd_exceptions.tsv`)

**Lineage-defining RDs behave as expected.** In the table, ST counts are ALT/REF/ABSENT/NOCALL, and
"Fitch" is the change count on the cohort tree.

| RD | expected | L1 | L2 | L3 | L4 | Fitch | AA_STATE / flag |
|---|---|---|---|---|---|---|---|
| RD239 | L1 | 250/0/0/0 | 0/250/0/0 | 0/249/0/1 | 0/250/0/0 | 1 (root TIED) | . |
| RD147c | L1 | 250/0/0/0 | 0/250/0/0 | 0/250/0/0 | 0/249/0/1 | 1 (TIED) | . |
| TbD1 (ACC record, 2,153 bp at 1,761,789, mmpL6) | present in L1 only | 250 present | 0/242 | 0/235 | 0/247 | 1 (TIED) | AA_SOURCE assumed |
| RD105 | L2.2 Beijing | 0/250 | 215 ALT, 35 NOCALL (L2.1) | 0/250 | 0/250 | 1 | REF / HOMOPLASTIC (panel 2) |
| RD105ext | L2.1 | 0/250 | 35 ALT / 215 REF | 0/249 | 0/250 | 1 | . |
| RD181 | L2.2.1 | 0 | 182/68 | 0 | 0 | 1 | . |
| RD311 | L2.2.1 "modern" | 0 | 161/89 | 0 | 0 | 1 | . |
| RD750 | L3 | 0/250 | 0/250 | 250/0 | 0/249/0/1 | 1 | . |
| RD316 | L3 | 0 | 0 | 250/0 | 0 | 1 | . |
| RD9 | L5/6/animal (negative control) | 0/250 | 0/250 | 0/250 | 0/250 | n/a | exact graph record, REF in all 1,000 |

Over all 61 RDs with expected carriers in the cohort: 59 have a matching record, 30 match perfectly
(every expected isolate ALT and no ALT outside), 44 match within 10% / 2%, and 37 need exactly one change
on the tree. RD207 is not in the RD file. Its CRISPR region is masked and fragmented into dozens of
caller intervals; the best candidate (3,119,536-3,122,522) is ALT in 208 L2 isolates but also in 48 L1
and 52 L3, so it is not usable as a check.

**Exceptions, by cause:**

- **Label vs TB-Profiler nomenclature, with the record still a clean clade (1 change):**
  - RD505 "1.2.2" is ALT in all 37 L1.2.1.2(.1) isolates and in none of the 71 L1.2.2.
  - RD524 "3.2" sits in 13 of 15 L3.1.3.
  - RD311 "2.2.1 modern" covers 161 of 182 L2.2.1; the 21 REF are one clade.
  - RD236a "1.1.1 except 1.1.1.2" covers 27 of 40.
  - RD713/RD317 (L5 labels) coincide with L3.1.2 / L4.4.1 deletions at the same coordinates.
- **Independent deletions at one locus folded into one interval record** (a real homoplasy hotspot,
  but the record looks multi-origin):
  - RD516 (separate L1, L2.1, L3 and L4.1 alleles): the 135 bp interval is ALT in L1, L2.1, L3 and L4.1
    together (Fitch 2-3). The allele-specific small_allele records do separate the alleles.
  - N-RD25 (L2.1 and L3 alleles share the 950 bp interval record).
  - RD5 family (phiRv2 prophage, IS6110-flanked, 80 ALT, 11 changes).
  - RD11 (EXCLUDE_REPEAT tier, 13 changes).
- **RD182 (L4.1.2.1, 6.5 kb) has no DEL record at all.** Its 22 L4.1.2.1 isolates are ABSENT on the
  IS6110 records inside it (2,547,053 and 2,550,014). The event shows up only as 67-68 bp
  small_allele INS records at its left breakpoint (2,545,196; 21/20 carriers, all L4.1.2.1). A
  lineage-defining deletion is therefore not written as a deletion. **This looks wrong.**
- **Partial within-lineage carriage:**
  - RD511 (L2.2.2): 17/33 ALT, 7 changes.
  - RD724 (L4.6.1) shares that region's interval and is also ALT in 18 L2.2.2 (9 changes).
  - RD167 (L1.2.1): only the 5 L1.2.1.1 isolates.
  - RD2_ER: no matching record.
  - RD513: 4 of 5 NOCALL (EV ambiguous_depth).
- **Out-of-clade singletons on otherwise clean RDs**, 1-5 isolates each: RD115, RD150, RD163, RD174,
  RD252, RD310, RD313, RD509, RD510, RD517, RD525. Most carry EV `.`; RD150 and RD174 carry
  `inherited_unrecorded_matched_only`.
- **Inconsistent nesting semantics.** On RD105, L2.1 isolates (which carry the larger RD105ext) are
  NOCALL. On RD516, isolates with an overlapping, differently sized deletion are ALT. Neither uses
  ABSENT (`*`), which the header defines for "region deleted in this sample".

## 2. Near-identical isolate pairs (`out/pair_*.tsv`, `odd_side*.tsv`, `private_alt_per_sample.tsv`)

There are 20 pairs at <= 5 SNPs (35 isolates; 8 L1, 5 L2, 1 L3, 6 L4) and 39 pairs at 6-12 SNPs
(65 isolates). All <= 5 pairs share the same matched reference, as do 36 of the 39 at 6-12.

Discordant cells (ALT vs REF) per 10,000 both-called cells:

| class | in repeat mask | <= 5 SNPs: called / discordant / per 1e4 | 6-12 SNPs: per 1e4 |
|---|---|---|---|
| DEL | outside | 10,441 / 2 / **1.9** | 2.9 |
| DEL | inside | 23,707 / 471 / **198.7** | 308.2 |
| INS | outside | 8,696 / 9 / 10.3 | 15.2 |
| INS | inside | 10,755 / 0 / 0 | 11.4 |
| IS6110 (H37Rv frame) | outside | 24,223 / 2 / 0.8 | 1.7 |
| IS6110 node frame | - | 6,857 / 0 / 0 | 1.6 |
| ACC | inside | 454 / 2 / 44 | 178 |
| ACC | outside | 645 / 0 / 0 | 0 |

- **Known RDs never differ within a <= 5 pair** (0 of 2,106 called RD cells). At 6-12 SNPs there is
  1 of 4,124: RD509's in-mask record, SAMN03647301 ALT vs SAMN03647464 REF.
- **Almost all DEL discordance is in-mask `interval` records.** These are 464 of the 473 discordant
  DEL cells at <= 5 SNPs, from 396 records (median 105 bp), and come from both graph (A/B tier) and
  caller intervals. Outside the mask, DEL discordance is 2 cells in 20 pairs.
- **The discordance is one-sided per pair.** In 998 of 1,021 resolvable in-mask DEL cells, the odd
  isolate is the ALT side: it carries a deletion its sublineage does not. Specific isolates drive this,
  for example SAMEA1117745 (138 cells), SAMEA5569818 (137), SAMN06010466 (89), SAMEA3558534 (89) and
  SAMEA95434918 (85). Read length does not explain it: some of these pairs are both 150 bp.
- **Per-isolate private in-mask DEL-interval ALT calls have a long tail:** median 5, 90th percentile 53,
  maximum 154. Outside the mask the median is 0 and the maximum 6. This is a candidate per-sample QC
  filter, in line with excluding noisy samples rather than patching the caller.
- **The remaining out-of-mask differences are mostly 50-75 bp tandem-repeat copy-number INS.** Examples
  are 2,990,584; 802,463; 4,348,7xx; 1,955,672; 2,613,924; and 39,030. The odd side is the ALT side in
  29 of 34. These are plausibly real VNTR changes (see section 3). There are also about 10 IS6110 sites,
  plausibly real transpositions, and one 457 bp DEL at 4,338,075 (SAMN21855968 vs SAMN33756512).
- **A different matched reference raises INS discordance:** 112 vs 5.6 per 1e4 at 6-12 SNPs. Most of it
  comes from one isolate, SAMN13568637, mapped to a different reference. It includes IS6110-length
  (1,358 bp) small_allele INS at 2,163,646 and 2,367,671.
- EV is `.` on nearly all discordant cells; the inherited_* EVs account for about 25 cells.

## 3. Ancestral reconstruction of SV records (`out/ancestry_records.tsv`, `ancestry_agreement.tsv`, `homoplastic_top.tsv`)

| class | mask | records | 1 change | >= 2 changes | median changes | derived | reversals |
|---|---|---|---|---|---|---|---|
| DEL | outside | 252 | 149 | 41% | 1 | 571 | 18 |
| DEL | inside | 1,304 | 373 | 71% | 4 | 25,366 | 5,058 |
| INS | outside | 244 | 124 | 49% | 1 | 670 | 196 |
| INS | inside | 1,113 | 539 | 52% | 2 | 2,540 | 894 |
| IS6110 | outside | 477 | 271 | 43% | 1 | 868 | 19 |
| ACC | all | 46 | 1 | 98% | - | 599 | 236 |

Outside the mask, deletions behave as irreversible events: 571 derived changes against 18 reversals.
Inside the mask, the 5,058 reversals and the median of 4 changes match the per-isolate noise found in
section 2.

**Cohort root vs panel AA_STATE** (records that carry an AA_STATE; tied roots excluded from the rates):

| class | mask / NOCALL share | records | tied | agree | disagree |
|---|---|---|---|---|---|
| DEL | outside | 51 | 0 | 51 | 0 |
| DEL | inside, NOCALL < 50% | 346 | 14 | 329 | 3 |
| DEL | inside, NOCALL >= 50% | 156 | 42 | 108 | 6 |
| INS | outside, NOCALL < 50% | 42 | 15 | 18 | 9 |
| INS | outside, NOCALL >= 50% | 6 | 0 | 1 | 5 |
| INS | inside, NOCALL < 50% | 283 | 39 | 220 | 24 |
| INS | inside, NOCALL >= 50% | 656 | 288 | 186 | 182 |

By AA_FLAG: DEL clean 335/336 agree, DEL HOMOPLASTIC 153/160, INS clean 292/438, INS HOMOPLASTIC
130/201. DEL polarity agrees essentially everywhere. INS disagreements are concentrated in two groups:

- Allele-split small_allele INS where about 90% of isolates are NOCALL (the carriers of the other
  alleles), so the root rests on few leaves.
- Tandem-repeat copy-number loci, where the cohort makes the longer allele ancestral (ALT root) and the
  panel says REF.

Most of the 220 ALT-root / REF-panel cases are in-mask PE_PGRS/PPE allele records. 3,566 records were
scored in all; 2,026 have no AA_STATE (AA_SOURCE assumed): 1,003 of 1,556 DEL (mostly
interval records), 370 of 1,357 INS, all 46 ACC and all 607 IS6110.

**The most homoplastic records outside the mask** (full list of 40 in `out/homoplastic_top.tsv`):

- **Tandem-repeat (VNTR-like) INS of about 53-58 bp units and their multiples:**

  | position | gene | changes |
  |---|---|---|
  | 2,613,926 | moeW | 44 |
  | 3,351,090 | | 23 |
  | 3,690,950 | | 19 |
  | 2,352,065 | | 19 |
  | 2,990,584 | aftC | 17 |
  | 1,779,278 | | 17 |
  | 39,017 | | 17 |
  | 4,348,724 | eccCa1 | 15 |
  | 2,372,436 | dop | 14 |
  | 802,427 | rplW | 13 |
  | 1,955,581 | | 12 |
  | 1,644,248 | | 9 |
  | 4,155,459 | leuA | 9 |
  | 3,232,650 | | 13 |

  Several of these positions coincide with known MIRU-VNTR loci (e.g. MIRU40 about 802 kb, MIRU16 about
  1,644 kb, Mtub21 about 1,955 kb, Mtub39 about 3,690 kb, MIRU39 about 4,348 kb, QUB-4156 about
  4,156 kb, QUB-3232 about 3,232 kb). Repeated copy-number change is expected there. None of them is
  IS6110-proximal.
- **plcD deletions** at 1,987,018-1,987,481 (175 / 463 / 95 bp; 25, 20 and 16 changes). They are
  IS6110PROX, 220-370 bp from H37Rv IS6110-3, at a cohort IS6110 insertion site: the known
  IS6110-mediated plcD deletion hotspot, so the homoplasy is real biology.
- **Small, about 70 bp deletions, not IS6110-proximal:** nrdB 279,515 (71 bp, 23 changes, 98 ALT),
  Rv1157c 1,283,485 (69 bp, 20 changes) and ftsK 3,059,787 (74 bp, 14 changes). They are most likely
  tandem-repeat contractions as well. They are worth a read check before use in convergence tests.

## Conclusions

1. RD inheritance is confirmed on the main lineage markers: RD239, RD147c, TbD1, RD105, RD181, RD750
   and RD316 each sit on one branch with no ALT outside their clade, and RD9 is REF everywhere. In
   near-identical pairs, known RDs never disagree.
2. Outside the repeat mask, SV calls are very reproducible between near-identical isolates: DEL about
   2 and IS6110 about 1 discordant cell per 10,000 called cells. The INS rate (about 10) is mostly VNTR
   copy number, which is likely real.
3. Inside the mask, DEL `interval` records disagree about 100 times more often (about 2-3%). Specific
   isolates adding private deletions explain this, and a per-isolate private-call count would flag them.
4. Things that look wrong or need a decision:
   - RD182 is not emitted as a DEL.
   - Nested and overlapping deletions are coded inconsistently (ALT on RD516's interval, NOCALL for
     RD105 in RD105ext carriers, never ABSENT).
   - Independent alleles at RD516 / N-RD25 / RD511-RD724 share one interval record.
   - 829 of 948 scored DEL `interval` records (graph and caller) carry no AA_STATE. That includes most
     lineage RDs (RD239, RD750, RD181 ...), so their polarity rests on the 'assumed' rule.
   - The RD file's sublineage labels do not match TB-Profiler names for several L1/L3 RDs (RD505, RD524,
     RD167).
