# Step 4 second review: integration (2026-10-06)

Area: integration of the eleven fix groups on `audit-fixes` (ffd7bb9).

What this covers:
- (a) coverage of every HIGH and MEDIUM audit finding, and a spot-check of D1-D44;
- (b) the interface contracts between the files several groups edited;
- (c) test-suite quality, including a skip census and mutation tests;
- (d) costs at 10,000 samples.

How it was done:
- Code was read from a `git archive audit-fixes` export at `/tmp/claude-12043/review2_integration/src`. File:line references below are to that export, which is the same as the branch.
- Scratch is in `/tmp/claude-12043/review2_integration/`. Nothing was written to the repository, the worktrees, runroot or netscratch, except this file.
- Parts of (a) and (b) were done by four read-only helper reviewers. I re-checked every finding below that is rated MEDIUM or above myself, in code and, where marked, on data.

## Summary

- **Correctness of individual fixes:** good.
  - 9 of 10 mutation tests are caught by the suite.
  - Every one of the 44 decisions is implemented as described.
  - No fix implements the wrong behaviour.
- **The defects are at the joins between groups:**
  - **The panel polarity table's indel keys are normalised differently from the merged VCF's (HIGH).** About 3,000 scale200 indels get no polarity lookup. 103 of them are left with the wrong ancestral allele.
  - **D6's "unmeasurable" presence state makes the level-2 conditioning blank every variant in a copy-number locus (MEDIUM).** 249 of 308 conditional variants lose all their events, silently.
  - **The P0 manifest omits `ancestral.tsv`, and a stale `manifest.done` survives an incomplete re-manifest (MEDIUM).**
  - **The association tail's provenance records are too narrow (MEDIUM).** Stale event matrices are reused, and the combined CX333 + cohort tree (TP-5, not addressed) cannot be fed in.
- **The suite runs 254 tests:**
  - 0 skipped with `config/project_env.sh` sourced;
  - 18 skipped silently with the MTB variables unset, including every `add_outgroup` and event-writer polarity test.
- **Scale:** no fix adds a prohibitive per-sample or per-record cost.

---

## Findings

### R2-INT-1 (HIGH): panel_polarity keys are not normalised like the merged VCF's, so indels lose or invert their polarity

**Where:**
- `assoc/bin/write_event_matrix.py:508-521`: exact lookup on `(chrom, pos, ref, alt)`.
- `bin/panel_polarity.py`: keys taken verbatim from `graph_collapsed.vcf.gz`.
- `bin/vcf_collapse.sh:89-95`: trims only, does not left-align (D22).
- `bin/p5_keys.py` / `bin/mtb_norm.py:25-58`: the merged VCF's keys are left-aligned against H37Rv.

**What is wrong:**
- D22 kept the collapsed graph VCF trim-only, so its indels sit where vcfwave put them, often right-shifted (e.g. `1817 G>GG`).
- The cohort keys are left-aligned (`1815 T>TG`).
- The GRAPHVCF-3/TP-8 fix made the event writer look polarity up per allele by exact key. A shifted indel therefore finds no row and stays "unpolarised", which means REF is taken as ancestral.

**Evidence (my reproduction, real data):**
- Inputs:
  - the scale200 `merged.vcf.gz` (its keys come from unchanged `mtb_norm`);
  - `panel_polarity.py` from audit-fixes, run on CX333 `all_variants.collapsed.vcf.gz` (`iface/pol.tsv`).
- Every panel indel key was left-aligned with `mtb_norm.normalise` against the build's H37Rv.
- Results for the 8,711 H37Rv-frame small indels:
  - 1,492 match exactly;
  - **2,951 match only after left-alignment;**
  - **103** of those 2,951 have the outgroup carrying the ALT at `panel_af >= 0.05`. They are written ALT-derived when ALT is ancestral, which is backwards.
- The helper reviewer measured 3,022 and 101 independently, by a different route.

**Effect:**
- About 100 indels per cohort are reconstructed with the wrong ancestral state, so their events sit on the wrong branches.
- About a third of the table's indel coverage is lost.
- `ancestral.tsv` and `add_outgroup` are SNP-only and are not affected.
- D22 recorded "some indels in repeats get a blank panel_af"; this consequence for polarity was not recorded.

**Suggested fix:**
- In `panel_polarity.py`, normalise each key with `mtb_norm.normalise(h37, pos, ref, alt)` (the function P5 uses), then pool any keys that collide.
- Alternatively, normalise both sides at lookup.
- Add a test with a right-shifted panel indel.
- `p5_matrix.load_panel_af` has the same mismatch, but only for the `panel_af` flag.

### R2-INT-2 (MEDIUM): with D6, level-2 conditioning blanks every variant inside a copy-number accessory locus

**Where:**
- `assoc/bin/write_event_matrix.py:779-783, 794-800, 826-834`;
- `assoc/bin/assoc_scan.py:303-320`;
- the same pattern in `bin/audit_chain.py`.

**What is wrong:**
- After the P3IS-2 fix (D6), `locus_presence.py` writes UNMEASURABLE for every read-route-blind locus, and the merge drops those loci.
- The level-2 readers still condition every node variant that `node_locus` places in such a locus. They see no PRESENT carriers, so `car_leaf.get(lid) is None`, which leads to `b[:] = UNK` and `n_applicable = 0`.
- The variant can never have an event and drops out of the scan and the audit silently.
- The merge applies `read_route_blind`; the level-2 readers do not.

**Evidence:**
- Excerpt:
  ```python
  keep_l = car_leaf.get(lid)
  if keep_l is None:
      b[:] = UNK            # locus never level-1 genotyped
      v["n_applicable"] = 0
  ```
- Helper reviewer, real data (scale200 node_locus and presence tables rewritten to UNMEASURABLE, which is what a fresh run produces):
  - 51 of the 61 placed loci are blind;
  - 249 of 308 conditional variants get `n_applicable = 0`, against 171 of 254 with the old tables.
- I re-read the code path. The counts are UNVERIFIED by me.

**Effect:** every node-frame variant inside an IS6110 copy or another copy-number locus is lost from the association test, with no count anywhere.

**Suggested fix:**
- In all three readers, load the build's catalogue and skip conditioning for a `read_route_blind` locus: leave `acc_locus` empty and the leaves as called.
- Alternatively, treat a locus with no measured cell as unconditioned.
- Count both cases.

### R2-INT-3 (MEDIUM): P0 manifest completeness and the build guards

**Where:** `bin/p0_prepare.sh:697-699` (step_manifest), `bin/p0_check.py` verify.

**What is wrong:**
- The completeness list omits `ancestral`, although `p5_finish.sh` requires `ancestral.tsv`. `manifest.done` can therefore be written before the ancestral table exists, and that table is then never re-hashed by `verify`.
- `verify` also does not cover `build_info.tsv`. The tail takes `build_id` and `outgroup` from it.
- On the incomplete path, step_manifest rewrites `manifest.tsv` but leaves an older `logs/manifest.done`, so `mtb_resolve_build` still treats the build as complete.

Found by the coverage and interface helpers; I confirmed the list in code.

**Effect:** the P0P2-2 guarantee ("the build is verified at use") does not hold for the ancestral table or the build stamp.

**Suggested fix:** add `ancestral` to the list, hash `build_info.tsv`, and `rm -f logs/manifest.done` on the incomplete path.

### R2-INT-4 (MEDIUM): the association tail's provenance is too narrow, and TP-5 is not addressed

**Where:**
- `assoc/bin/cohort_assoc_tail.sh:109-128` (`_prov` / `_current`);
- `:152-162` (the tree is always rebuilt from the cohort's `og.fasta`).

**What is wrong:**
- `.prov` records only build_id and vcf_sha, plus polarity_sha for the events.
- Running once without presence tables and then adding them reuses the event matrix with no level 2, while step 5 passes `--accessory-presence` to the scan.
- Replacing the tree also keeps the old events.
- An externally made `rooted.nwk` with no `.prov` is refused (FATAL), so the combined CX333 + cohort tree route (HANDOFF 0b, TP-5) cannot be used.
- TP-5 itself (`analysis/combined_tree/build_alignment.py`: last duplicate wins, missing reads as REF, reads the decomposed VCF) is untracked and unchanged, and no decision covers it.

**Effect:** results become silently stale on rerun, and the production tree route is blocked.

**Suggested fix:**
- Record the sha256 of the tree, `og.fasta`, sites, node_locus and the presence-table list (with their sidecars).
- Add a `TREE=` input that records its own prov.
- Bring the combined-tree route into the repository, or record a decision on TP-5.

### R2-INT-5 (MEDIUM, pre-existing, not a regression): p3acc can start before P1 writes the refmap

**Where:**
- `bin/refbias_run.sh:452-454`: submitted with an empty dependency.
- `accessory/bin/locus_presence_array.sh:41` and `locus_presence_one.sh:17-22`: read the sample and reference from `$REFMAP`.

**What is wrong:** on a fresh full-chain submission the presence tasks fail on a missing refmap, and `p5vcf`, which is `afterok:p3acc`, then never runs. The array is also sized from the cohort table but indexed by refmap row.

**Suggested fix:** depend on `${JOB[p1]}` and size the array from the refmap. Worth fixing before the 10K run.

### R2-INT-6 (LOW-MEDIUM): build identity is not checked on the IS6110 inputs to the merge

**Where:**
- `bin/merge_cohort_vcf.py:545-597, 605-613`;
- `bin/p1i_vcf.sh:48` and `bin/p1i_p5states.sh:49` (`RES=is6110/results`).

**What is wrong:**
- The IS6110 key table and the stage-2 table are found by cohort name under `is6110/results/`. That folder is outside the `_stamp_root`-guarded output root.
- The key table carries a `build_id` column that nothing reads.
- If a cohort name is reused on a new build and p1iv is not rerun, the merge joins the old build's node keys.

**Suggested fix:** in the merge, refuse key-table rows whose `build_id` is not `--build-id`. Also write the stage-2 table's build into it.

### R2-INT-7 (LOW): sharded and unsharded merges can differ at shard edges

**Where:** `bin/merge_cohort_vcf.py:950-983`.

**What is wrong:** `dup_of` compares a caller insertion only with this shard's IS6110 keys and level-1 records. A duplicate within 10-20 bp of a shard boundary therefore survives a sharded run.

**Effect:** at most a few records at 20 shards.

**Suggested fix:** build `is_starts` from the whole key table, and `acc_starts` from the whole catalogue (minus blind loci), before the shard filter.

### R2-INT-8 (LOW): remaining interface loose ends

- **ASSOC-5 is only partly fixed.**
  - `is6110_gene_burden.py:49-63` (`deleted_span`) parses `svi:` IDs and a `<DEL:-len>` ALT. Caller deletions are written with ID `sv:DEL:<pos>:<len>` and ALT `<DEL>`, so the span fallback never fires and they are credited at the anchor.
  - The old scale200 merged VCF holds 2,947 UNCATALOGUED caller deletions (my count).
  - Today they are presence-only and have no gains, so there is no current effect. Parse `sv:DEL:` IDs.
- **Mask not taken from the build.** `cohort_assoc_tail.sh:236-240` calls `retier_intervals.py` without `--mask`/`--rds`, so it reads the working tree's `data/annotation/H37Rv_repeat_mask.bed` (`retier_intervals.py:103-104`), not `${BUILD}/assets/repeat_mask.bed`.
- **Node-contig sibling groups (TP-7, node frame).**
  - Merged-VCF small records on `node_<id>` contigs are all at POS=1, while IS6110 node records are at offset+1.
  - `vcf_to_alignment.with_siblings` groups by the ID-derived locus with `itertools.groupby`, so a node:offset group that is not contiguous is split.
  - The old scale200 VCF has 8 such cases in 11,740 node records.
- **Reference tip at node-frame sites.** `vcf_to_alignment.py:134-135` writes REF for the `--ref-sample` tip at node-frame sites. The event writer, after TP-3, treats those sites as unknown (`.`). The tree alignment should write N there too; 20 such columns in scale200.
- **D41 allele orientation.** `p4_place.on_strand` complements, and for indels re-anchors, node-frame alleles whenever the projection strand toggles, which now includes the odgi inverted-step toggle. The node offset stays at R's original anchor. For a reverse-strand indel the key's offset and its alleles therefore describe different anchor bases. Keys stay self-consistent between P4 and P5 (both use the same transform), so this is a semantics issue that should be settled together with D41. UNVERIFIED on data.
- **IS6110 assets outside the build.** The IS6110 arm still reads crossmaps and GFFs from `is6110/assets/isclean_matched` and `matched_gff`, which are outside the build and only partly covered by `p0_check` (`is6110_intervals`). A new panel's references must be added there by hand.
- **Inconsistent build pick.** `p5_finish.sh` picks its build with `find -type d` when `MTB_BUILD_DIR` is unset, not with `mtb_resolve_build`, so a half-built directory can be chosen when it is the only one.

### R2-INT-9 (LOW-MEDIUM): test-suite gaps

- **Silent skips without the environment.** With the MTB variables unset, the suite reports OK with 18 skips:
  - all 4 `AddOutgroup` tests;
  - 5 `EventWriter` tests (TP-2, TP-3, TP-6 and the polarity table);
  - `PanelPolarityTable`;
  - 2 IS6110 merged-VCF tests (P3IS-5, -6);
  - 3 sharded-merge tests;
  - 2 alignment-archive tests;
  - 1 SV merge test.

  With the environment sourced, 0 are skipped.

  Tool discovery is inconsistent:
  - `test_audit_cleanup2.bcftools()` and `test_audit_rerun_safety.bcftools()` fall back to the lab path and so run;
  - `test_audit_polarity` reads only `MTB_BCFTOOLS` and skips.

  The PGB test files error under the system python3 (3.6, no `capture_output`).

  Suggest that `run_tests.py` fails, rather than skips, when `MTB_PY`/`MTB_BCFTOOLS`/`MTB_BGZIP`/`MTB_PY_VT` are missing, unless `--allow-skips` is given.
- **D15 is untested.** Mutation M5 removed `o["q_branch"] < 0.05` from `assoc_scan.survives()` and the full suite still passed. `SurvivorRule` tests only the missing-lineage half.
- **No test covers D2 or D4** (`p5_states.py:749-752` and `:613-615`).
- **Wiring tests check source text.** About 25 of the 221 new tests assert that a string appears in a script's source, for example that `p5_merge.sh` contains `--node-paths "$NODE_PATHS"`. They catch a deleted line but not wrong semantics, and no test runs a `.sh` wrapper chain end to end.

---

## (a) Coverage of HIGH and MEDIUM audit findings

"Fixed" was checked in code on every path the finding names. Tests were read to confirm they assert the fixed behaviour.

| ID | sev | status | where (audit-fixes) | test | note |
|---|---|---|---|---|---|
| GRAPHVCF-1 / Fault A | HIGH | fixed | add_outgroup.py:106-111 | test_audit_polarity.AddOutgroup.test_any_exact_duplicate_with_gt1_is_alt | D35 |
| GRAPHVCF-2 | MED | fixed | add_outgroup.py:112-142, 201-204 | AddOutgroup (3 more) | D11 |
| GRAPHVCF-3 / TP-8 | MED | fixed, but see R2-INT-1 | panel_polarity.py:87-111; p0_prepare.sh:587-610; write_event_matrix.py:493-521 | PanelPolarityTable, EventWriter | D13 |
| GRAPHVCF-4 | MED | fixed | t8_select_reference.py:43-104 | test_audit_selection.SelectionByAllele | D20, D21, D24 documented |
| GRAPHVCF-5 | MED | fixed | vcf_collapse.sh:90-156 | CollapseKeepsAllelesTrimmed, DecomposeMatchesProduction | D22, D23 |
| Fault B | — | fixed | ancestral_alleles.py:110-150 | AncestralFaultB (2) | D44 |
| TP-1 | MED | fixed | add_outgroup.py:90-142 | test_upstream_mnp_and_deletion | |
| TP-2 | MED | fixed | vcf_to_alignment.py:83-117; write_event_matrix.py:222-227 | test_star_record_is_a_column, _reads_the_outgroup | D12; mutation M8 caught |
| TP-3 | MED | fixed | write_event_matrix.py:230-240, 820, 924 | test_reference_tip_unknown_at_node_frame | D14; tree alignment still REF (R2-INT-8) |
| TP-4 | MED | fixed (guards), residual | p0_prepare.sh:144-186, 631-664; cohort_assoc_tail.sh:109-128 | P0 / AssocTailGuard tests | prov too narrow (R2-INT-4) |
| TP-5 | MED | **not addressed** | — | — | R2-INT-4 |
| TP-6 | MED | fixed | write_event_matrix.py:741-743 | test_dedupe_first_keeps_outgroup_aligned | |
| ASSOC-1 | HIGH | decision (deferred by user) | — | — | |
| ASSOC-2 | MED | fixed | assoc_scan.py:72-87, 640; burden :307 | LeaveOneOut | |
| ASSOC-3 | MED | fixed | assoc_scan.py:59-69, 447; burden :255 | DedupeSuffixKey (3) | |
| ASSOC-4 | MED | fixed | cohort_assoc_tail.sh:231-251; assoc_scan.py:495-517 | SvTiers (3), ChainScript | D19; mask (R2-INT-8) |
| ASSOC-5 | MED | **partly** | is6110_gene_burden.py:49-63 | SvSpanCredit (svi only) | caller `sv:DEL` (R2-INT-8) |
| ASSOC-6 | MED | fixed (repo side) | cohort_assoc_tail.sh:4-56 | ChainScript | working-tree copies remain |
| P0P2-1 | HIGH | fixed | p1_select_reference.sh:113-185; p1_summary.py; p2_call.sh:84-196; refbias_run.sh:300-334 | P1Guard, P2Guard, RunnerGuard | D38 |
| P0P2-2 | MED | fixed, residual | p0_prepare.sh:176, 423-482; refbias_run.sh:155-165 | BuildVerify, P0 | R2-INT-3, retier mask |
| P0P2-3 | MED | fixed | p0_prepare.sh:59-77, 440-443 | P0 tests | |
| P0P2-4 | MED | fixed | p0_prepare.sh:568-579; p1_summary.py | P1Summary (3) | |
| P3IS-1 | HIGH | fixed | is6110_write_vcf.py:300-340; is6110_p5_stage2.py:111-128, 395 | Is6110CohortFrameState, Is6110Stage2Occupied | D9, D14; mutations M2 and M10 caught |
| P3IS-2 | MED | fixed in VCF; **breaks level 2** | locus_presence.py:49-62, 192-204; merge_cohort_vcf.py:911 | MergeInsertionsAndPresence | D6; R2-INT-2; M1 caught |
| P3IS-3 | MED | fixed | is6110_project_sites.py:139-165; is6110_write_vcf.py:213-219 | Is6110RepeatedNode (3) | D8 |
| P3IS-4 | MED | fixed | is6110_repeat_rescue.py:175-194 | Is6110RescueOffset | |
| P3IS-8 | MED latent | partly | locus_presence_array.sh:46-87; refbias_run.sh _stamp_root | PresenceGuard, RunnerGuard | elementdepth not in P3 guard; consumers do not check sidecars |
| P4P5-1 | MED | fixed | p5_states.py:748-779 | ProjectionAndReferenceAllele (5) | D1, D2; M9 caught; D2 untested |
| P4P5-2 | MED | fixed | p5_states.py:586-625 | test_insertion_r_carries_is_alt_not_ref | D3, D4; D4 untested |
| P4P5-3 | MED | fixed | p5_sv_genotype.py:290-296 | test_clip_cluster_does_not_override_full_depth | D5; M7 caught |
| P4P5-4 | MED | fixed | p4b_place_sv.py:53-59 | test_hom_ref_records_are_not_read | rerun P4b into new folders |
| P4P5-5 | MED | fixed | sv_intervals.py:131-155 | IntervalCatalogue.test_no_widening | |
| P4P5-6 | MED | fixed | merge_cohort_vcf.py:942-1015 | MergeInsertionsAndPresence | D7; shard edge (R2-INT-7) |
| P4P5-7 | HIGH (rebuild) | fixed | p0_prepare.sh:546-560; p5_states.py:633-651; p5_merge.sh | NodeTables, P4P5ReadTheBuildsCollapsedVcf | D34, D36; a missing node gives a silent NOCALL |
| P4P5-8 | MED | fixed | p4_place.py:409-462 | InheritedOverlap | M6 caught |
| PGB-1 | HIGH | fixed | pggb_build.sh:179-244 | PggbCleanBuild (3) | |
| PGB-2 | HIGH | fixed in code | pggb_build.sh:4, 109-163, 226 | PggbSettings (4) | INPUTS.md stale |
| PGB-3 | HIGH | fixed | build_panel.py; make_pansn_fasta.sh:42-51 | BuildPanelFromRecordedDecisions (7) | reproduction of the 333 UNVERIFIED |
| PGB-4 | HIGH | fixed | assembly_provenance_screen.py:116-242 | ProvenanceRule (3) | D28, D30 |
| PGB-5 | HIGH | fixed except UniVec | foreign_insertion_screen.py | ForeignScreen (5) | D29, D32, D33; D31 only partly enforced in code |
| PGB-6 | MED (current data) | **partly: detection only** | assembly_qc_stats.py:69-97 | SingleBaseRuns | no blacklist or masking for current outputs; no D# |
| PGB-7 | MED | partly | assembly_qc_stats.py; assembly_provenance_screen.py; variant_counts.py | several | repeat_checks / hybrid_vs_sr not ported |
| PGB-8 | MED | fixed in code | variant_counts.py; snp_outlier_screen.py | OutlierScreen | `--accessions` optional |
| PGB-9 | MED | fixed for repository readers | add_outgroup, panel_polarity, sv_intervals, merge_catalogues, p0_prepare, sync_back | ReadersTakeCollapsed etc. | 2 working-tree scripts still default to decomposed |
| PGB-10 | MED | fixed | vcf_decompose.sh:88-186 | DecomposeMatchesProduction | |
| PGB-11 | MED | **not addressed** (comment only) | vcf_decompose.sh:74-76 | — | p4b_place_sv.py:294, 345 still branch on LV; no D# |
| PGB-12 | MED | fixed | snp_nonredundant.py:44-88 | NonRedundantInput (3) | D27 |
| PGB-13 | MED | partly | refresh_assemblies.sh:83-98 | RefseqSelection | GenBank not searched; no D# |
| PGB-14 | MED | fixed (main hazard) | stamp_build_id.sh:70-98 | Stamp | exits 0 with no build |
| PGB-15 | MED | decision D37 | p0_prepare.sh:232-236, 377-400 | RunnerGuard | |
| repeat-mask provenance | B | partly | p0_prepare.sh:412-459 | BuildVerify | retier reads the working-tree mask |

**D1-D44:** all 44 are implemented as DECISIONS.md describes.
- File:line evidence for each is in the helper's table, in the scratch transcript.
- Synthetic runs were done for D5, D6, D11, D23, D39, D40 and D44.
- Caveats:
  - D41's alleles follow the projection-strand transform rather than "each reference's own orientation" (R2-INT-8).
  - D43: a mistyped `INSGT_*` path is treated as absent, not refused (`merge_catalogues.py:101-105`).
  - D38: the accessory folder is adopted without evidence when it holds no stamped VCF, and the stamped-VCF check samples only 20 files.
- Of the "found during the fixes, still open" bullets:
  - the odgi strand off-by-one, the retier off-by-one, the `mixed` header, p5_sanity NA and the ISMapper default are fixed;
  - burden indels at the anchor are partly fixed, by D39 and D40;
  - the rerun prerequisites are still operational, and are now enforced by guards.

## (b) Interface contracts checked and found sound

| hand-off | producer | consumer | result |
|---|---|---|---|
| P4 → P5 keys | p4_place.py writes `node_offset` = forward offset (`mtb_norm.forward_offset`) and key `node:<n>:<fwd>` | p5_keys.py `node_key(node, node_offset, …)`; p5_states.py `by_key` uses the same; non-carrier read at `start + forward_offset(off, strand, L)` | consistent. P4 drops a reverse-walked record with no length, and p4_place.sh requires the node table |
| IS6110 producers → key table | project_sites (forward offset, node_occ via `odgi paths -H`) → place_by_flank (carries node_occ) → write_vcf (`frame` h37rv/node/repeat_node; state ALT/REF per H37Rv occupancy; duplicate fatal) | p5_merge, stage2 (`carriers`, `occupied`, `short`), merge_cohort_vcf (`repeat_node` skipped; duplicate conflict fatal; EVIDENCE/SITECLASS `mixed`) | consistent key strings (`h37rv:<p>`, `node:<n>:<off>`, short form `…:`). Stage 2 projects from carrier R positions, not node offsets |
| P5 states → merge | sparse-v1 / `states.u8.npy` with keys_sha1 | read_array / read_sparse; place() uses POS=1 for node small variants (node_offset dropped by design) | consistent |
| P4b / svgt / intervals → merge | `sv:` keys, `svi:DEL:<first>:<len>` | sv_state keyed on `key`; iv_states on `interval`; retier `interval` | consistent (svi ID = retier interval, checked on scale200 by the helper) |
| accessory catalogue → presence → merge | `ACC_%07d`, UNMEASURABLE | merge `acc:<id>`, blind loci dropped | consistent, **except the level-2 readers (R2-INT-2)** |
| merged VCF → vcf_to_alignment → add_outgroup | sites `column chrom pos ref alt` | add_outgroup row per site | real-data check (helper): 6,073 columns, every row the same length, 30,365 random cells match the VCF GT, outgroup differs from bcftools only where N is written by design |
| → write_event_matrix → assoc_scan / burden | row_key `locus\|REF\|altkey[\|#k]`, `primary_alt` for `X,*` | `event_key` round-trips | consistent; the event writer ran on full scale200 with the tail's arguments |
| build-id propagation | build_info.tsv → `--build-id` (p4, p4b, p1i_vcf, p5_finish → merge, required); stamp_build_id; `_stamp_root` | cohort_assoc_tail compares the VCF stamp with the build | consistent; gaps in R2-INT-3 and R2-INT-6 |
| required arguments | every new `required=True` (merge `--build-id`, is6110_write_vcf `--refs/--h37rv/--build-id`, is6110_p5_merge `--h37rv`, add_outgroup `--panel-vcf/--outgroup`, write_event_matrix `--panel-polarity`, p5_states `--node-*`, panel_polarity, graph_frame_offsets) | callers: p5_finish.sh, p1i_vcf.sh, p1i_p5states.sh, p5_merge.sh, cohort_assoc_tail.sh, p0_prepare.sh, tests | all passed. `snp_nonredundant.py --rank-by` has no repository caller (hand-run) |
| guards on a fresh run | refbias_run `_stamp_root` on a new outroot, p2/p1 markers, p5_finish freshness, locus_presence sidecar | — | a clean new outroot and accessory folder are created and recorded (cohort tables live outside the outroot). The legacy CX333 build cannot be brought up to date without a manual reset (helper: date-only `.done` markers, missing panel_manifest.genomes.txt); a new build runs clean |

## (c) Tests

**Full suite:**
- With the environment sourced (`$MTB_PY` 3.11): 254 tests, OK, 0 skipped. 68 s on an idle node; 2-6 min under login-node load.
- With MTB variables unset (`env -i`, same interpreter): OK, 18 skipped silently (list in R2-INT-9).
- The skip set is deterministic. I found no test that passes vacuously: no early `return`, no swallowed exception.

**Mutation tests:** one line of a fix reverted in a scratch copy, then the full suite run with the environment sourced.

| # | mutation | caught by |
|---|---|---|
| M1 | merge ignores `read_route_blind` (D6) | test_audit_p4p5.MergeInsertionsAndPresence.test_merge |
| M2 | IS6110 carrier state back to `REF if shared` (P3IS-1) | test_audit_is6110.Is6110CohortFrameState |
| M3 | P4 inverted-step `t -= 1` removed | test_audit_leftovers.OdgiInvertedStepP4 (2) |
| M4 | `forward_offset` returns the walking offset (D41) | test_audit_cleanup ForwardOffset, P4NodeKey (2), IS6110NodeKey, P5ReadsForwardOffset |
| **M5** | **`q_branch < 0.05` removed from the survivor rule (D15)** | **nothing; the suite passes** |
| M6 | inherited-overlap check off (P4P5-8) | InheritedOverlap |
| M7 | full depth + clips → ALT (P4P5-3, D5) | TwoFrameConfirmation |
| M8 | `X,*` records dropped again (TP-2) | AlignmentStarAndSiblings.test_star_record_is_a_column |
| M9 | always ABSENT off R's path (P4P5-1, D1) | ProjectionAndReferenceAllele.test_unrelated_target_is_nocall_not_absent |
| M10 | occupied-locus non-carrier back to REF (D9) | Is6110Stage2Occupied |

Script and results: `/tmp/claude-12043/review2_integration/mut/{run.py,results.txt}`.

## (d) Scale (10,000 samples)

No fix introduced a per-sample or per-record cost that is prohibitive at 10,000 samples.
- **Node table:** `node_positions.tsv` grows from 438K rows (old accessory table) to 1.38M rows under D36.
  - Each P5 task now loads it whole: about 2 s CPU and about 470 MB RSS, measured.
  - That is about 6 CPU-h across 10K tasks, within the current 2 GB task request.
  - P4 and the IS6110 projection read it filtered.
- **`odgi paths -H`** (P3IS-3) runs once per cohort over 333 paths × about 276K nodes (CX333). Acceptable.
- **`leave_one_out` (ASSOC-2)** is O(|pool|) Python per variant per rung, the same order as the buggy one-liner it replaced.
  - That is about 10^9-10^10 steps at 10K (helper estimate, 1-3 h).
  - It is not new. A bincount precompute would remove it.
- **The new merge blocks** use bisect only, and the level-1 block reads every presence table once per shard.
- **Pre-existing costs worth planning for** (not introduced by the fixes; helper estimates):
  - `assoc_scan` `_dep_cache`: n_leaves × 20,000 float32 per conditional locus, about 48 GB at 10K.
  - `vcf_to_alignment`: about 3×10^9 pysam calls and Python character lists of about 12-16 GB, run in the foreground of the tail.
  - IQ-TREE at 10K tips within `build_snp_tree.sh`'s 24 h allocation.

## Audit findings not addressed

- **TP-5:** the combined tree route is outside the repository and keeps the add_outgroup faults. The new tail guard now refuses its tree.
- **PGB-11:** only the comment changed; P4b still branches on LV.
- **PGB-6:** detection only, with no blacklist or masking for current outputs and no decision.
- **Partly addressed, with no decision recorded:**
  - ASSOC-5 (caller deletions);
  - PGB-7;
  - PGB-13 (GenBank);
  - P3IS-8 (elementdepth; consumers do not check presence sidecars).

## Open questions

1. **R2-INT-1:** should panel polarity be normalised in `panel_polarity.py` (left-align against H37Rv, then pool), or should D22 be revisited (`norm -f` in the collapse)?
2. **R2-INT-2:** for a copy-number locus the read route cannot see, should within-insert variants be unconditioned (the conservative reading), or excluded with a count?
3. **D41 allele orientation for node-frame indels** (R2-INT-8): decide together with the pending D41 alternative.
4. **UNVERIFIED:**
   - an end-to-end real-data trace of an IS6110 insertion, an SV deletion and an accessory locus through P4 → merge → events → scan on new code. It needs P4/P5 reruns, which are out of scope without jobs; only the merged VCF → alignment → outgroup → event-matrix leg was run on real data;
   - the helper's R2-INT-2 counts (249/308);
   - whether IS6110 GT=2 at `ref_shared` sites means element loss, which the writer's default `--absent unknown` would hide (scale200: 6,832 such cells).
