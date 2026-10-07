# Decisions to confirm: audit fixes, section A (2026-10-05)

All five fix groups are merged on branch `audit-fixes` (pushed, not yet merged
to main). 117 tests pass together.

Each group implemented the conservative option, meaning "unknown" rather than
a guess, wherever correct behaviour needed a scientific choice. **Each item
below is implemented as described in the "Implemented" column. Confirm it, or
pick the alternative.** The figures come from the groups' validation runs on
copies of production data.

## The user's answers (review started 2026-10-07)

| Group | Items | Answer |
|---|---|---|
| Calling and genotyping | D1-D7 | all confirmed as implemented (2026-10-07) |
| IS6110 | D8-D10 | all confirmed as implemented; D8's alternative (key by the nearest single-copy node) tracked for the new build in HANDOFF 0q (2026-10-07) |
| Node-frame alleles | D41 | option b, done on `audit-fixes` (2026-10-06) |

## Calling and genotyping (P4-P5)

| # | Choice | Implemented | Alternative | Size |
|---|---|---|---|---|
| D1 | A projected H37Rv position with no clean junction in the reference | NOCALL | ABSENT (old behaviour), or ABSENT whenever the odgi distance is over 5 | Of the old ABSENT SNP cells, about half are confirmed deletions and stay ABSENT; the rest go NOCALL or are genotyped |
| D2 | No projection at all | NOCALL | ABSENT | |
| D3 | The reference carries a third allele, or matches neither haplotype | NOCALL | REF (old) | **The largest volume:** about 600 cells per sample, mostly masked or PE indels in homopolymers and repeats |
| D4 | A deletion of 50 bp or more that the reference carries, evidenced only by a reference block | NOCALL | ALT | |
| D5 | Full H37Rv-frame depth plus a clip cluster | NOCALL. Also, H37Rv-frame NOCALL plus clips gives NOCALL. | REF; ALT when both ends are clipped | scale200 1,918, gwas1000 14,087 cells leave ALT |
| D6 | Accessory loci H37Rv already carries | unmeasurable, no record (including PRESENT calls: 763 / 4,760 ALT cells) | keep PRESENT as a presence-only ALT | level-1 records 802 → 128 |
| D7 | A caller insertion duplicating an IS6110 or accessory record | dropped | kept and flagged INFO/DUPOF | 157 / 209 rows |

## IS6110

| # | Choice | Implemented | Alternative | Size |
|---|---|---|---|---|
| D8 | Insertions on graph nodes their reference visits more than once | not keyed: they leave the cohort records, counted and reported | key them by the nearest single-copy node (a new method) | **scale200: 456 of 772 node-frame rows; gwas1000: 1,945 of 3,585** |
| D9 | Sites where H37Rv already holds a copy | carriers REF; non-carriers NOCALL, so the record is effectively uninformative | emit as a loss of H37Rv's copy (`<DEL:ME:IS6110>`); or drop these keys | 2,028 / 11,725 non-carrier cells REF → NOCALL |
| D10 | Conflicting duplicate (sample, key) rows | stop with an error | ALT over REF | 0 duplicates in current data |

## Ancestral states and the outgroup (polarity)

| # | Choice | Implemented | Alternative | Size |
|---|---|---|---|---|
| D11 | Outgroup uncalled over a covering record | N | REF (old "no record" rule) | 86 / 87 cells |
| D12 | `*` alleles and sibling-allele carriers count toward the alignment's 10% missing limit | yes | no: this would admit about 1,400 more X,* SNPs (scale200) and keep 37 panel sites | |
| D13 | Panel allele frequency in `panel_polarity` | carriers / genomes called, from GT | sum(AC) / max(AN) from INFO | fourth-decimal changes |
| D14 | The H37Rv tip at node-frame IS6110 records | stays REF | unknown, as for other node-frame classes | 207 records |

## Association tests

| # | Choice | Implemented | Alternative | Size |
|---|---|---|---|---|
| D15 | **Survivor rule** | requires q_branch < 0.05 for every variant (HANDOFF's "all three nulls"); a missing lineage null counts as a failure | the old rule (region + lineage only) | 0 survivors change today, with or without the canettii isolate |
| D16 | Deletion overlap needed to credit a gene | any overlap, 1 bp or more | at least 50% of the gene; or whole-gene removal only | SV-burden units: scale200 179 → 263 |
| D17 | Burdens get the scan's 80% callability floor | yes | report floor-passing origins alongside instead | scale200 small units 5,189 → 5,140 |
| D18 | Point variants in overlapping genes | the earlier-starting gene only | credit every overlapping gene (17,874 bp of H37Rv counted twice) | |
| D19 | SV evidence tiers | each cohort tiered from its own genotypes | tier every cohort from the largest cohort | |

## Reference selection and the graph VCF

| # | Choice | Implemented | Alternative | Size |
|---|---|---|---|---|
| D20 | Normalise the selection distance by compared sites | no (documented) | report d/n_compared; or score only sites both sides called | 17 of 997 gwas1000 rank-1 choices would change |
| D21 | H37Rv as a possible matched reference | not a candidate | add it as an all-REF column | 1 scale200 isolate |
| D22 | Left-align indels in the collapsed graph VCF | no, trim only (positions stay vcfwave's) | `norm -f`: 8,934 records move, 41 events re-duplicate | some indels in repeats get a blank panel_af |
| D23 | Duplicate records disagreeing 0 against `.` | 0 wins | `.` wins | |
| D24 | Uncovered isolate sites in P1 selection | still counted REF (CODE_REVIEW 6.11) | needs isolate coverage data | |

## Deferred by the user

- **ASSOC-1, the canettii isolate:** see `analysis/canettii_effect/README.md`.
  Survivors are unchanged without it; it inflates only the branch null.

## Found during the fixes, still open

- **Possible odgi strand off-by-one in P4:**
  - At inverted path steps, odgi's position output is one base off and
    complemented. `p5_states` now corrects for this.
  - `p4_place`'s projection ignores the strand column. Whether it has the
    same off-by-one is unverified. It affects 13,508 of GCF_000193185's
    67,901 keys.
  - **To check in step 3.**
- **Indels in the small and IS6110 burdens** are still credited by their
  anchor base.
- **`retier_intervals.py`:** a possible off-by-one in its repeat and RD
  overlap.
- **The merged-VCF header descriptions** need `mixed` added for IS6110
  EVIDENCE and SITECLASS.
- **`p5_sanity.py`** should write NA, not 0, when the H37Rv counts are
  missing.
- **`is6110_project_sites.py`** defaults its ISMapper directory to the pilot's
  (`refbias/p1f`).
- **Rerun prerequisites:**
  - regenerate `panel_polarity.tsv`, which is now required;
  - P0 ancestral (move aside the old guard file);
  - p1iv in full (the new `node_occ` column);
  - interval IDs change (twoframe and svgt rerun);
  - P1 to P5 into new output folders.

---

# Section B fixes (step 3): additional decisions

All merged on `audit-fixes` (200 tests pass, f66dc57).

## Graph build

| # | Choice | Implemented | Alternative |
|---|---|---|---|
| D25 | Sparse mapping in `pggb_build.sh` | off by default; explicit `-x` only | default `-x auto`, after a record-by-record comparison of test arms B and D |
| D26 | The durable mirror keeps the decomposed graph VCF (175 MB) | yes | drop it; rebuilding the collapsed file then costs a vcfwave rerun of up to 12 h |
| D27 | `snp_nonredundant.py --rank-by` | required | keep the accession tie-break, documented |

## Panel construction

| # | Choice | Implemented | Alternative |
|---|---|---|---|
| D28 | Unknown sequencing technology | passes; noted for review (18 genomes) | SHORT_READ_SUSPECT, or a required manual call |
| D29 | An insert counts as "native" only if at least 2 other genomes carry it | yes; otherwise REVIEW | a vector database (UniVec, not available locally); position rules |
| D30 | IS6110 profile site tolerance | 50 bp (10 bp broke 46 production calls) | |
| D31 | The new screens never exclude a genome on their own | yes; every flag needs a recorded decision | automatic exclusion |
| D32 | Foreign-screen background | first accession per sublineage, plus genomes with no call | chosen by quality |
| D33 | The canettii outgroup in the foreign screen | not exempted in code; reviewed by hand | an exemption list |

## Build safety

| # | Choice | Implemented | Alternative |
|---|---|---|---|
| D34 | Panel allele-frequency source | the build's collapsed graph VCF (identical values where both exist; 180 / 268 more keys) | also copy the nolab file when present |
| D35 | `add_outgroup` panel VCF | the build's collapsed VCF (scale200: 51 of 46,379 cells move between called and N; no base changes) | copy the decomposed VCF into the build |
| D36 | Node table scope | nodes off the H37Rv path (complete for every current key) | all nodes (about 60M rows) |
| D37 | References | the deposited sequences, not the graph's rotated paths. All 333 are identical up to rotation and strand; GFFs, crossmaps and the projection store rely on the deposited frame. | stage the references from the graph FASTA and rotate the GFFs (a new method) |
| D38 | Outputs from before the guards existed | adopted only on evidence (VCF stamp, GATK reference, P2 stamps, timestamps); everything else refused | refuse all legacy outputs |

## Leftovers

| # | Choice | Implemented | Alternative |
|---|---|---|---|
| D39 | Complex or MNP records in the small burden | stay at the anchor (14 / 15 records) | credit POS to POS+len(REF)−1 |
| D40 | A small deletion inside overlapping genes | credited to both genes, as SV deletions are; this departs from D18 for small deletions | D18's one-gene rule |

## Found and fixed in step 3 (not in the audit)

- **The odgi strand off-by-one in P4 and P4b was confirmed and fixed.**
  - At inverted path steps, P4 placed called records one base off with the
    uncomplemented base.
  - The placed REF matched H37Rv at only 22 of 59 SNPs (chance level) before
    the fix, and 59 of 59 after.
  - Affected: scale200 234 records in 52 samples; gwas1000 about 1,660
    records in 276 samples.
  - In P4b, 1,285 scale200 SV keys move by 1-2 bp.
- **Node-offset orientation** (in the clean-up group, now running): the same
  base can get two `node:<id>:<offset>` keys when references walk a node in
  opposite directions. 13 / 71 nodes are affected.

## Clean-up pass (merged, fb199f4; 222 tests)

**Fixed: node keys now use the node's forward offset,** so one base has one
key whichever way a reference walks the node.

- **Wrong offsets before:** the old walking offset named the wrong base at
  345 scale200 and 2,097 gwas1000 P4 node-frame records.
- **P5 cells read at the wrong base before:** 536 of 119,480 in scale200 and
  15,972 of 1,594,916 in gwas1000.

| # | Choice | Implemented | Alternative |
|---|---|---|---|
| D41 | **Node-frame alleles** | offsets fixed only; alleles stay in each reference's own orientation, so the same event read from opposite walks still has complementary alleles (C>G vs G>C) and two full keys (11 of 11 shared bases in scale200) | restate node-frame alleles on the node's forward strand (complement when walked '-'; re-anchor indels; matching flip in p5_states) — **recommended as a separate change; needs your decision** |
| D42 | A record on a reverse-walked node with no length available | no key: dropped and counted (IS6110: fatal) | a P0 node table covering every node |
| D43 | Accessory catalogue inputs | no cohort census; insgt routing only when given | |
| D44 | Ancestral outgroups (`ANC_OUTGROUPS`) | optional; the default is both canettii (fatal if they are not tree leaves) | required |
