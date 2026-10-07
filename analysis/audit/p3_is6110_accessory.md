# Audit: P3, IS6110 and accessory presence (`p3_is6110_accessory`)

Audited 2026-10-05, read-only. Scratch: `/tmp/claude-12043/audit_p3/`.

**Which copy runs.** `runroot/is6110/bin`, `runroot/accessory/bin` and
`runroot/graphframe/bin` are symlinks to the repository. Every repository
script in this area differs from its working-tree namesake, so production runs
the repository copies. The current key tables carry the `h37rv_pos_placed`
column, which only the repository writer emits, so they came from repository
code.

**`is6110/bin/panisa_score.py`.** Production does not call it. In the
repository it appears only in the docstring of `is6110_junctions.py:17`. In the
working tree, `p1c_panisa.sh` and `p1c_summary.py` use it, but pass p1c is not
in `refbias_run.sh`'s `ORDER`.

**Counts.**

| severity | findings |
|---|---:|
| HIGH | 1 |
| MEDIUM | 4 |
| LOW | 5 |

---

## P3IS-1 (HIGH): IS6110 carriers that share their matched reference's copy are written GT=0 in the H37Rv-frame VCF

**Where** (repository copies run):

- `is6110/bin/is6110_write_vcf.py:283-286`: `state="REF" if shared else "ALT"`;
- `is6110/bin/is6110_p5_merge.py:91-93,120`: the carrier's state is copied as is;
- `bin/merge_cohort_vcf.py:560` plus `STATE_GT` at line 41: REF becomes `0`.

**What is wrong.** `site_class=ref_shared` means the isolate's reads show an
IS6110 junction at a seam where its matched reference also carries a copy. The
isolate carries the element. The arm writes this as `REF`, which is correct
only in the matched-reference frame (the per-sample `.is6110.vcf` rightly
writes no record). The cohort key, however, is in the H37Rv frame (or a node
frame), and the merged record is `ALT=<INS:ME:IS6110>`. So:

- **At H37Rv-empty keys:** GT=0 says "no insertion". These carriers become
  indistinguishable from the stage-2 non-carriers, whose REF does mean "looked,
  no insertion" (`is6110_p5_stage2.py`, `depth >= 5` and no junction).
- **At H37Rv-occupied keys,** where H37Rv itself has the element, the polarity
  is inverted for everyone except the shared carriers:
  - a carrier whose reference lacks the copy is written ALT;
  - stage-2 REF means "this sample lacks the element", although H37Rv has it.

The design note `is6110/docs/P5_MERGE_PLAN.md` section 1 calls the matched-frame
REF "settled and correct". The error comes from carrying that state into a
record in a different frame.

**Evidence.** Measured on `refbias/<c>/p5/merged.vcf.gz` joined to
`is6110/results/<c>_p1i_cohort_keys.tsv`:

| cell kind | scale200 | gwas1000 |
|---|---:|---:|
| **H37Rv-empty keys** | | |
| carriers written 1/ALT (reference lacks the copy) | 469 | 2,424 |
| carriers written **0/REF** (reference shares the copy) | **1,098** | **5,103** |
| stage-2 non-carriers, 0/REF | 99,310 | 1,375,858 |
| **Node-frame keys** | | |
| carriers written 0/REF | 639 | 2,797 |
| carriers written 1/ALT | 126 | 738 |
| **H37Rv-occupied keys** | | |
| carriers written 0/REF | 97 | 414 |
| carriers written 1/ALT, although H37Rv has the element | 2 | 22 |
| non-carriers written 0/REF (they lack the element H37Rv has) | 2,028 | 11,725 |

In both cohorts, about 70% of carrier cells at H37Rv-empty keys are written as
not carrying.

**Direct contradiction.** Some keys hold both kinds of carrier for the same
insertion. Their genotype differs only by which reference each isolate was
matched to:

| | scale200 | gwas1000 |
|---|---:|---:|
| keys with both kinds of carrier | 38 | 104 |
| REF carriers in those keys | 251 | 2,373 |
| ALT carriers in those keys | 69 | 351 |

Records with at least one shared carrier: 318 of 778 (scale200) and 477 of
2,020 (gwas1000).

**Effect on current outputs:**

- the IS6110 block of both merged VCFs;
- the association chain (`assoc_scan.py --classes` includes `is6110`), and
  `is6110_gene_burden.py`;
- any convergence or event count on IS6110 records.

Because matched references are chosen as close relatives, insertions shared by
a clade are mostly in that clade's reference. Clade-level insertions therefore
appear mostly as REF, and the ALT calls are biased toward isolates whose
reference happens to lack the copy.

**Suggested fix.** State carriage in the frame of the record, and keep the
matched-frame class as information only.

```
# is6110_write_vcf.py, keys_rows
- state="REF" if shared else "ALT",
+ # cohort frame: a carrier carries the element, whatever its reference holds
+ state="ALT",
```

- For `h37rv_state=occupied` keys, either:
  - emit the record as a loss of H37Rv's copy: carriers REF; stage-2 "no
    element" ALT, or a separate `<DEL:ME:IS6110>` record; or
  - exclude these keys from the `<INS>` block.
- Keep `site_class` as a per-sample FORMAT field, not as the genotype.
- Rerun p1iv, p1is and p5vcf for both cohorts.

---

## P3IS-2 (MEDIUM): level-1 accessory presence is emitted for 674 loci the read route cannot see

**Where:**

- `accessory/bin/locus_presence.py:175-232`, the read route: unmapped plus
  clipped reads from the H37Rv CRAM;
- `bin/merge_cohort_vcf.py:866-905`, which emits every catalogue locus, maps
  ABSENT to REF, and whose comment asserts "These loci are sequence present in
  panel genomes and absent from H37Rv".

**What is wrong.** 674 of the 802 catalogue loci are `novelty=copy_number`, and
all 674 have `h37rv_cov >= 0.9`. Their reads place on H37Rv (MAPQ 0 for the
IS6110 loci), so they never enter the unmapped or clipped pool. The read route
is blind to them, and they are written ABSENT (GT 0, "does not carry").

481 of them are 1,340-1,370 bp IS6110 element sequences, nearly identical to
one another (505 distinct sequences among the 674). Even when IS6110 reads do
reach the pool, a read cannot say which of the 481 "loci" it belongs to.

**Evidence** (all `accessory/<c>/*.presence.tsv`):

| IS6110-length copy_number loci (481) | scale200 | gwas1000 |
|---|---:|---:|
| cells | 96,200 | 479,557 |
| PRESENT | **0** | **0** |
| ABSENT | 96,112 | 478,878 |
| ABSENT although `reference_carries=1` | 899 | 4,639 |

All 802 loci are in the merged VCF: `ACCKLASS` `polymorphic` 249 and
`reference_gap` 553 in scale200.

In scale200, 120 cells are REF on these presence records while ALT on an
IS6110 record within 20 bp.

The copy_number class also makes up most of the polymorphic level-1
characters:

| | scale200 | gwas1000 |
|---|---:|---:|
| polymorphic level-1 characters | 57 | 110 |
| of which copy_number class | 31 | 76 |

What a PRESENT call means for a locus whose sequence H37Rv already has is
**UNVERIFIED**: it may be a divergent homolog, or it may be noise.

**Effect on current outputs:**

- The 481 IS6110 loci are monomorphic, so association drops them. The VCF
  deliverable still states about 575,000 (gwas1000) and about 96,000
  (scale200) unmeasured "not carried" cells.
- The 31 and 76 copy_number characters do enter `assoc_scan`
  (`accessory_presence` class).

**Suggested fix.** In `merge_cohort_vcf.py` level 1, skip
`novelty=copy_number` loci, or emit them NOCALL. Copy number belongs to the
`<CNV>` route that `ACCESSORY_IN_VCF_PLAN.md` describes. Correct the comment at
866-868. Validate the remaining copy_number "PRESENT" calls before using them
as characters.

---

## P3IS-3 (MEDIUM): node-frame IS6110 keys on repeated nodes merge unrelated insertions

**Where:**

- `is6110/bin/is6110_project_sites.py:177-210`: the node key comes from
  `odgi position -v` on the junction base;
- `is6110/bin/is6110_write_vcf.py:217-220,254-260`: a node key is used whenever
  flank placement fails;
- `is6110/bin/is6110_p5_stage2.py:95-108`: `carriers()` projects each key from
  its first carrier only.

**What is wrong.** This is pattern 2, a node taken as an identity whatever it
contains. A node key is an identity only if the node occurs once in each path.
Node 46966 is a 1 bp `G` node that every one of the 333 paths traverses 85 to
453 times (`accessory/assets/node_positions.tsv`).

**Evidence:**

| | scale200 | gwas1000 |
|---|---:|---:|
| `node:46966:0` rows | 58 | 239 |
| `node:46966:0` samples | 51 | 190 |
| node-frame rows on a node the carrier's own path visits more than once | 442 of 772 | 1,861 of 3,585 |
| node keys involved | 56 | 98 |

Within one reference, `node:46966:0` collects sites more than 1.1 Mb apart
(for example SAMEA104447025, r_pos 1,994,050 and 2,002,910; SAMEA104679297,
2,002,438 and 3,134,325). All of them become one merged-VCF record.

Stage 2 then projects that key from a single carrier locus and calls REF or
ABSENT for every other sample from it.

The 50 duplicate (sample, key) pairs noted in HANDOFF section 7 are all
node-frame, and 23 or more of them are on node 46966 (see P3IS-5).

**Effect on current outputs:** node-frame IS6110 records in both cohorts. One
record in each cohort is pure artefact; the rest are of uncertain identity.

**Suggested fix:**

- In `is6110_write_vcf.py`, refuse a node key when the node has
  `n_occurrences > 1` in the carrier's path, or is shorter than about 30 bp.
  Fall back to a per-sample unplaced key, or walk to the nearest single-copy
  node and key on (node, signed distance).
- Count and report the sites dropped this way.

---

## P3IS-4 (MEDIUM, 1 record now): DR-array rescue converts the array start with the wrong offset

**Where:** `is6110/bin/is6110_repeat_rescue.py:171-172,180`.

```
near = min(inside, key=lambda r: abs(int(r["clean_pos"]) - L[0]))
off = int(near["orig_pos"]) - int(near["clean_pos"])
agg.update(clean_pos=L[0], orig_pos=L[0] + off, ...)
```

**What is wrong.** The clean-to-original offset is borrowed from the nearest
stack. When the reference's own IS6110 copy is excised between that stack and
the array's first base, the offset is wrong by the element's length.

**Evidence.** Recomputing `clean_to_orig` with each reference's crossmap for
every rescued row:

- scale200: 5 rescued rows, 0 wrong;
- gwas1000: 24 rescued rows, 1 wrong.

The wrong row: SAMEA104679294 (GCF_965121955) has orig_pos 2,385,638 written,
against 2,384,283 from the crossmap. 2,385,638 lies inside that reference's own
excised element (2,384,921-2,386,275), so the seam test fails. The row was
keyed `ref_lacking`, ALT, `node:82943:0`. It is probably the reference's own
DR copy.

**Effect:** one gwas1000 key-table row and VCF cell. It recurs whenever a
reference's DR array contains an excised element.

**Suggested fix:**

```
- off = int(near["orig_pos"]) - int(near["clean_pos"])
- ... orig_pos=L[0] + off
+ ... orig_pos=clean_to_orig(L[0], crossmap_rows_for(R))
```

---

## P3IS-5 (LOW): duplicate (sample, key) rows, last row wins

**Where:**

- `is6110/bin/is6110_p5_merge.py:93`: `observed[(sample, k)] = state`;
- `bin/merge_cohort_vcf.py:560`: `is6110[k][sample] = state`.

**Evidence:**

| | scale200 | gwas1000 |
|---|---:|---:|
| duplicate pairs | 7 | 50 |
| REF/REF | 6 | 49 |
| ALT/REF | 1 | 1 |

All the gwas1000 duplicates are node-frame. They are mostly a symptom of
P3IS-3. The ALT/REF pairs are decided by row order.

**Suggested fix:**

- Resolve explicitly, ALT over REF, and count it.
- Fail if two distinct `r_pos` of one sample share an H37Rv key.

---

## P3IS-6 (LOW): key-level INFO comes from the first carrier's row; orientation is not recorded

**Where:**

- `bin/merge_cohort_vcf.py:561` (`is_meta.setdefault(k, r)`) and 935-940;
- `is6110/bin/is6110_write_vcf.py:262`.

**What is wrong:**

- `EVIDENCE` and `SITECLASS` describe whichever sample came first, although 236
  gwas1000 keys mix evidence classes (57 in scale200) and 104 mix site class
  (38 in scale200). `assoc/bin/write_event_matrix.py:215-216` reads both.
- `MEINFO=...,+` is hard-coded. Review item 4.10 is still open, although the
  `p_*_left/right` pair counts in `elstacks.tsv` carry the orientation.

**Suggested fix:**

- Write per-key summaries, such as `N_TWO_SIDED` and `N_ONE_SIDED`, or move
  these fields to FORMAT.
- Derive polarity from the majority pair class, or write `.`.

---

## P3IS-7 (LOW): p1iv runs reconcile without `--crossmap-dir`

**Where:** `bin/p1i_vcf.sh:101-104`.

**What is wrong.** Without the flag, `is6110_reconcile.py` uses H37Rv's
crossmap (`--crossmap` default) for every sample. Three things follow:

1. The `chrom_side_only` rows' `orig_pos` is in the wrong frame: 166 rows in
   scale200, 762 in gwas1000.
2. The review 4.10 guard ("missing crossmap is fatal") never runs.
3. The ISMapper join, which defaults to the pilot's `refbias/p1f`, compares
   matched-frame positions with H37Rv positions.

**Effect:** none on the VCF. The writer and projection drop rows with no
geometry, and `ismapper` is 0 or blank on every row. The reconcile tables and
summaries are wrong for anyone reading them.

**Fix:** pass `--crossmap-dir is6110/assets/isclean_matched --ismapper-dir ""`.

---

## P3IS-8 (MEDIUM, latent): accessory presence and P3 skip guards ignore the build and catalogue

**Where:**

- `accessory/bin/locus_presence_array.sh:46`: `[[ -s "$OUT" ]] && exit 0`;
- `bin/p3_accessory.sh:137`, which checks `accessory.tsv` and `copynumber.tsv`
  only.

**What is wrong.**

- Presence tables are kept whenever they exist, whatever catalogue they were
  made against.
- HANDOFF section 4.2 moves aside `p3`, `p4`, `p4b` and `p5`, but not
  `accessory/<cohort>`.
- On the new graph planned in HANDOFF section 0l, the catalogue is rebuilt. The
  old tables would then be merged by `locus_id` (`ACC_<H37Rv pos>`) against the
  new catalogue without any warning: an id that has moved means a different
  sequence, and a new id becomes NOCALL for everyone.
- The P3 guard omits `elementdepth.tsv`, which is written last. A task killed
  in the element-depth step is skipped on rerun. That output is unused, so this
  matters less.

**Effect now:** none, because the build is unchanged. It applies at the
rebuild.

**Fix:**

- Write the catalogue's sha256 into each presence table, or a sidecar, and
  redo the table on mismatch.
- Add `accessory/<cohort>` to the move-aside list.
- Add `elementdepth.tsv` to the P3 guard.

---

## P3IS-9 (LOW): P3 arm B and element depth (no downstream consumer)

`refbias/io_contract.tsv` lists `elementdepth.tsv` as "computed and currently
unused". No script reads `copynumber.tsv` or `accessory.tsv` except
`p3_summary.py`.

- **`bin/p3_copynumber.py:62,101-104`.** The review 4.8 fix chooses the window
  direction from column 4 after `frame_convert from-panel`. That column is only
  the global storage flip; odgi's own local strand is moved to column 5 and
  ignored. A locus inside a local inversion is still read in the wrong
  direction. The window is also L+1 bases. Not measured.
- **`p3_copynumber.py:77`, `p3_element_depth.py:64`.** Depth is keyed by
  position with no contig. This is safe today, because all 333 paths are single
  contigs, but it would break on a multi-contig reference.
- **`bin/p3_summary.py:173-180`.** Missing samples are listed but the script
  exits 0 (pattern 3).

**Fix:** use the XOR of both strand columns, or project both ends and sort, as
`p3_element_depth.py` does. Key depth by (contig, pos). Exit nonzero on missing
samples.

---

## P3IS-10 (LOW): documentation against code

- `bin/merge_cohort_vcf.py:14-16` says accessory regions are "NOT emitted yet",
  but level-1 records are emitted (802 per cohort). Lines 866-868 are covered
  under P3IS-2.
- The presence thresholds in code are absent at 0.20 or below and present at
  0.80 or above (`locus_presence.py:138-143`). This matches
  `accessory/TWO_LEVEL_ACCESSORY.md`, which says they are untested.
  `refbias/ACCESSORY_CALIBRATION_RESULTS.md` ships "present at 0.90" for a
  different instrument (all reads against the matched reference plus elements).
  The 0.80 cut has no calibration behind it.

---

## Checked and found sound

**Frames:**

- `graph_frame.to_panel` and `to_refs` are exact inverses, and the flipped
  offset derivation in `graph_frame_offsets.py:96-98` is correct.
- The build frame table resolves all 333 accessions, with `panel_len ==
  refs_len` for every one.
- `frame_detect` complements bases for flipped accessions.

**Coordinate conversions:**

- P3 anchor and GFF conversions are 1-based to 0-based for odgi and back.
- Threaded odgi output is keyed by source, not by order.
- `p1i_matched.sh` and `p1g_isclean.sh` clip scan: SAM POS and end are both
  1-based. `is6110_junctions.py` handles the SA far side and its strand.
- `element_side.junction_of` picks the boundary from the SA CIGAR.
- `is6110_build_isclean`: excision, `clean_to_orig` and `orig_to_clean` are
  round-trip checked. `element_side.to_orig` matches it.
- `discover_elements` converts PAF 0-based half-open to 1-based inclusive, with
  strand-aware extension.

**Seams and flank placement:**

- One 3 bp seam rule (`is6110_seam.py`) is used by promotion, flank placement
  and the writer.
- `place_by_flank` window and gap arithmetic is correct for empty and occupied
  loci, and `h37rv_pos` is correct on both strands.

**Shared-key clustering (HANDOFF section 2.1, `cluster_keys`):**

- The window scan covers every open cluster.
- No cluster contains two sites of one isolate.
- No key mixes empty and occupied H37Rv states.
- Moved distances are 1-6 bp: 443 rows in scale200, 2,306 in gwas1000.

**Earlier review fixes, re-verified:**

- `--all-stacks` is passed (4.3), and an unkeyed site is fatal.
- Stage-2 `short()` keeps the node offset (4.1).
- An odgi or samtools failure is fatal, and an unmeasured cell is NOCALL (4.2).
- The projection store is keyed (target, carrier, position) and build-scoped.
- `is6110_p5_merge` builds keys through `mtb_norm`.
- `locus_presence.py` checks every subprocess, refuses an empty pool, and
  writes atomically (4.7).
- `is6110_junctions.py` always writes a header (4.6), and P1i's done test needs
  all three outputs.
- All 333 crossmaps and matched GFFs exist.

**Completeness:** every refmap sample of scale200 (200) and gwas1000 (997) has
P1i junctions, elstacks and elementdepth; P3 accessory, copynumber and
elementdepth; and a presence table. No stale `.lock.*` directories exist.

**`build_repeat_mask.py`:** the k-mer and tandem flags, BED output (0-based
half-open), and the PE/PPE label routing check.

**`node_path_membership.py`:** per-path offsets and occurrence counts.

## Open questions (not settled)

- **Copy_number PRESENT calls.** Are the 31 (scale200) and 76 (gwas1000)
  copy_number level-1 characters with PRESENT calls real divergent homologs or
  noise? UNVERIFIED.
- **Node offsets on reverse traversal.** For a node a reference traverses in
  reverse, are the offsets from `odgi position -v` used in IS6110 node keys
  oriented the same way across references? Not checked. It would add split keys
  on top of P3IS-3.
- **Local inversions in P3 copy number.** How many P3 copy-number windows fall
  in local inversions? This needs odgi and was not run (P3IS-9).
- **Gene and element coordinates in `build_repeat_mask.py`.** Gene coordinates
  from the snpEff dump and element GFF coordinates are used as 0-based (lines
  227 and 215). That is a possible 1 bp shift that affects labels and gene
  extension edges only. The dump's coordinate base was not verified.
- **Build-time scripts.** `is6110_isclean_summary.py` (P1g summary only),
  `p1i_build_matched.sh` and `p1i_discover_matched.sh` were read only briefly.
  They are not in the chain.
