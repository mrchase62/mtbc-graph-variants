# Review 2: IS6110 and accessory (`is6110_accessory`)

2026-10-06. Branch `audit-fixes` (ffd7bb9) against `main`. Read-only. Scratch:
`/tmp/claude-12043/review2_is6110/`.

**What was run.** I took a `git archive audit-fixes` export (`code/`) and built
a scratch runroot (`rr/`): code from the export, data read-only from
production. On scale200 I ran the whole IS6110 tail with the branch's code:

- rescue, reconcile, project_sites (with `--node-lengths`), place_by_flank,
  write_vcf;
- stage 1, stage 2 (with a copy of the build's projection store), and the
  IS6110 block of `merge_cohort_vcf.py`;
- `write_event_matrix.py` on the scale200 tree.

The log is `scale200.log`; the outputs are in `scale200/`. For the D8 numbers
on gwas1000 I reused the is6110 fix group's chain outputs
(`/tmp/claude-12043/fix_is6110/{scale200,gwas1000}`). Their
repeat_node / verdict / site_class columns do not depend on the later offset
changes, and my scale200 run reproduces their counts exactly (456 repeat_node,
316 node, 1,666 H37Rv).

---

## Findings

### R2-IS-1 (MEDIUM): D6 drops accessory loci the read route does measure, including some clearly variable ones

**Where:**

- `accessory/bin/locus_presence.py:52-62` (`read_route_blind`);
- `bin/merge_cohort_vcf.py:911` (applies it);
- the classification it relies on, `bin/build_accessory_panel.py:118-161`.

**What is wrong.** "Blind" is `novelty == copy_number or h37rv_cov >= 0.9`.
`h37rv_cov` is the union of all blastn HSPs at e-value 1e-10, at **any**
identity. So a divergent homolog (for example 85% identity across its length)
scores 1.0 and is classed as "H37Rv sequence". Reads from such a locus do not
align cleanly to H37Rv. They reach the unmapped and clipped pool, so the read
route can see them. The audit's premise ("their reads place on H37Rv") holds
only for near-identical copies.

**Evidence.** I re-blasted all 802 representatives against the build's H37Rv
and counted the fraction of each locus covered at or above a given identity
(`acc_identity.tsv`). Of the 674 blind loci:

| identity threshold | <50% of the locus covered |
|---|---:|
| ≥99% | 85 loci |
| ≥95% | 60 loci |

The second column is "not H37Rv sequence in any sense reads would see". For a
blind locus, comparing the old read-route calls with the independent
`reference_carries` column shows real measurement
(`blind_poly_{scale200,gwas1000}.tsv`, gwas1000):

| locus | len | identity cov ≥95% | PRESENT if ref carries | ABSENT if ref carries | PRESENT / ABSENT if ref lacks |
|---|---:|---:|---:|---:|---:|
| ACC_2867346 | 656 | 0.00 | 164 | 1 | 245 / 377 |
| ACC_2165937 | 3437 | 0.19 | 130 | 2 | 517 / 217 |
| ACC_0334653 | 5887 | 0.00 | 123 | 32 | 3 / 2 |

In all three, carriers' samples are called PRESENT almost always, so the calls
are not noise. These are polymorphic (klass `polymorphic`), measurable loci that
the branch now writes no record for. Across the 60 divergent blind loci there
are 24 (gwas1000) and 20 (scale200) with both PRESENT and ABSENT calls: 1,918 /
328 PRESENT cells. A further near-identical locus, ACC_3691060 (identity cov
1.0), also tracks its reference (113 of 117 PRESENT have ref carrying). That one
is unexplained (see open questions).

The 481 IS6110-length loci and the bulk of the 674 are correctly blind: 599
are ≥90% covered at ≥95% identity.

**Effect.** A minority of genuinely variable accessory characters leave the
merged VCF and the scan. This is the D6 "no genuinely variable locus dropped"
check: it fails for about 3 strong loci and about 20 weaker ones per cohort.

**Fix.** Define blindness by identity, not homology. For example, call a locus
blind when at least 0.9 of it is covered by H37Rv HSPs at ≥95-97% identity.
Compute this in `build_accessory_panel.py` as a new column and read that column
in `read_route_blind`. Validate with the `reference_carries` concordance above.

---

### R2-IS-2 (MEDIUM): D8 excludes far more than the ambiguous rows, and the lost signal is large and clade-level

**Where:**

- `is6110/bin/is6110_write_vcf.py:217-218` (`node_ok` needs `node_occ == "1"`);
- `:325-327` (repeat_node, no key).

**What was checked.** For the rows D8 excludes, I grouped the old keys by node
(`node_placed`). Using each carrier's own on-path graph H37Rv projection, I
asked whether a key really merged different loci:

- "conflated": spread over 1 kb, or one sample twice;
- "consistent": all carriers within 1 kb, one site each.

| | scale200 | gwas1000 |
|---|---:|---:|
| repeat_node rows (D8) | 456 | 1,945 |
| ... two-sided | 352 | 1,443 |
| ... graph placement on_path | 416 | 1,789 |
| ... ref_shared | 415 | (not split) |
| rows on conflated keys (truly ambiguous) | 257 | 1,227 |
| rows on consistent multi-carrier keys | 167 | 657 |
| singleton keys | 32 | 61 |
| samples losing ≥1 site | 175 / 200 | 788 / 996 |

In scale200 the median affected sample loses 20% of its IS6110 sites. The lost
rows form about 68 loci (50 bp clusters of graph position). They are highly
recurrent: the largest have 43, 38, 34 and 29 carriers. They are mostly copies
the matched reference also holds, sitting in IS6110 hotspots (around 1.99 Mb,
2.63 Mb and 3.12 Mb) where flank placement fails (`flanks_disagree_gap` 226,
`one_flank_unique` 188).

**Effect:**

- About half the excluded rows really were conflated, so excluding them fixes
  P3IS-3 as intended.
- The other half (167 + 32 scale200; 657 + 61 gwas1000) were not shown to be
  ambiguous by this test. They are lost as well, together with every
  clade-level insertion at those hotspots.
- No wrong genotype is introduced: excluded carriers are not written REF at
  their own sites.

This is a decision (D8), not a bug. But the "456 / 1,945 rows" figure in
DECISIONS understates it: this is about 19% of all IS6110 sites, and it removes
a whole class of recurrent loci.

**Suggestion.** Keep D8 for the release if needed, but treat the alternative
(key by the nearest single-copy node, or by the carrier's own flank-anchored
position) as a priority, not an option. At a minimum, report the excluded
fraction per sample in the run summary.

---

### R2-IS-3 (LOW, latent MEDIUM at scale): an on-path, reverse-walked site keyed on a node gets offset 0

**Where:**

- `is6110/bin/is6110_project_sites.py:264-265`: node lengths are requested only
  where `dist != 0`;
- `:312-313`: `noff = ""` for on-path sites with no length;
- `is6110/bin/is6110_write_vcf.py:235` and `:277`: `int(fr.get("node_offset") or 0)`.

**What is wrong.** The writer keys any site whose **flank** placement fails on
its node, including sites that odgi placed on the H37Rv path. For an on-path
site on a node the carrier walks in reverse, project_sites leaves the forward
offset blank. It requests no length, and D36's node table excludes
H37Rv-path nodes anyway. The writer then reads the blank as offset 0. Three
consequences:

- the key and the node-contig POS name the wrong base;
- sites at different bases of one node collapse onto `:0`;
- a carrier walking the same node forward gets a different key (a split).

**Evidence (my scale200 run):**

- 1,178 on-path rows have a blank offset (409 of them with `node_occ = 1`);
- 51 of the 316 node-frame key rows are keyed `node:<id>:0` from a blank
  offset, for example `node:77075:0`, whose walking offset was 12, and
  `node:180962:0`, whose walking offset was 67.

In scale200 every multi-carrier key of this kind happens to hold carriers at
one H37Rv position, and no node mixes a blank-offset key with another key. So
there is no wrong merge today; the offsets are simply wrong. The same path has
not been checked on gwas1000 (UNVERIFIED).

**Fix.**

- Either resolve the forward offset for every site that could become a node
  key: take node lengths for H37Rv-path nodes from the GFA S-lines or
  `odgi`, not from the D36 table.
- Or treat a blank offset in the writer as not keyable (repeat_node-style,
  counted), and never as 0.

Add a test with an on-path reverse-walk row whose flank verdict fails.

---

### R2-IS-4 (LOW): level-2 conditioning is applied to IS6110 node-frame records; on blind loci it silently removes them

**Where:** `assoc/bin/write_event_matrix.py:787-795` and `:826-835`.

**What is wrong.** Every record whose ID starts `node:` is looked up in
`node_locus`, including IS6110 records (ID `node:<nid>:<off>`). Two cases:

- **The locus is blind (D6).** The regenerated presence tables have no
  PRESENT cells for it, so `car_leaf.get(lid)` is None and every leaf becomes
  UNK. The record disappears from the scan without being counted as such.
- **The locus is mosaic.** IS6110 carriers not called PRESENT are masked.

**Evidence.** 10 of the 145 IS6110 node keys in my scale200 run sit on nodes in
`node_locus_scale200.tsv`. 7 of those are on copy_number loci (5 on
ACC_2262192) and 3 on mosaic loci.

D6 also turns level 2 off for the small variants on 147 of 196 locus-assigned
nodes in scale200. 130 of those were already all-UNK, because their locus had
0 PRESENT; 17 change.

**Fix.**

- Restrict the level-2 lookup to `cls == "small"`, or state the intent for
  IS6110.
- For blind loci, leave variants unconditioned (or use the P5 ABSENT cells as
  "inapplicable") rather than all-UNK, and count them.

---

### R2-IS-5 (LOW): stage 2's odgi projection still ignores the strand column

**Where:** `is6110/bin/is6110_p5_stage2.py:244-255`.

**What is wrong.** The leftovers fix (`d721cbd`) corrected odgi's
"one base high at an inverted step" in P4, P4b and `is6110_project_sites.py`.
`run_odgi` here still reads `tp[1]` without checking `f[3]`, and the build
store keeps the uncorrected values.

**Effect.** At inverted steps the depth position and the junction window
centre are 1 bp off. With a 10 bp window and a depth test this is negligible.
I list it only so that "every odgi reader is strand-corrected" is not assumed.

---

### R2-IS-6 (LOW): occupied-key records carry no marker in the merged VCF

**Where:** `bin/merge_cohort_vcf.py:1042-1048`.

**What is wrong.** Under D9 the records at H37Rv-occupied keys hold carriers
GT 0, everyone else `.` or `2`. The INFO has no `H37RV_STATE`, so a VCF reader
cannot tell that GT 0 here means "carries the element H37Rv has", while on
every other `<INS:ME:IS6110>` record it means "no insertion".

**Scale.** 14 records in scale200, all monomorphic, with no tree effect. This
is a deliverable-clarity issue only.

**Fix.** Write `H37RV_STATE=occupied` (and add a header line), or drop these
keys (the D9 alternative).

---

### R2-IS-7 (LOW): shared crossmap and GFF folders are rewritten per build

**Where:**

- `is6110/bin/p1i_discover_matched.sh:69-75`;
- `is6110/bin/p1i_build_matched.sh:99-104`.

**What is wrong.** The new `.ref.sha256` guard correctly redoes an accession
whose reference bytes changed. But `is6110/assets/{matched_gff,isclean_matched}`
are shared by every build, so building build B overwrites build A's crossmap
and isclean FASTA for any accession whose bytes differ. Build A's later stage 2
and writer then read B's crossmap silently.

On the first run after this change, every existing GFF and crossmap (none has
a sidecar) is regenerated, which is cheap.

**Fix.** Make OUTDIR build-scoped (`<build>/is6110/...`), or refuse when the
sha differs and another build is recorded as the owner.

---

## Checked and found sound

**P3IS-1 end to end (scale200, branch code, real data).** Merged-VCF cells by
key class:

| key class | cell | count |
|---|---|---:|
| H37Rv-empty | carriers, reference lacks the copy | 469 ALT |
| H37Rv-empty | carriers, reference shares the copy | **1,098 ALT** |
| H37Rv-empty | stage-2 non-carriers | 99,310 REF / 8,721 ABSENT / 1,802 NOCALL |
| node | carriers | 92 + 224 ALT |
| node | stage-2 non-carriers | 16,480 REF / 10,853 ABSENT / 1,351 `.` |
| H37Rv-occupied | carriers | 2 + 97 REF |
| H37Rv-occupied | non-carriers | 2,296 NOCALL / 405 ABSENT |

- The 1,098 shared-carrier cells were REF before; this is exactly the audit's
  count.
- The occupied non-carriers are D9's NOCALL.
- No record holds carriers of both states, so the "direct contradiction"
  keys are gone.
- No ALT appears at any occupied key.

**Event direction on the tree** (`write_event_matrix`, scale200 tree,
outgroup set, H37Rv tip REF), IS6110 records:

| | old (production) | new (branch) |
|---|---:|---:|
| records with 0 gains | 289 | 66 |
| records with ≥2 gains | 34 | 48 |
| gains | 536 | 719 |
| losses | 1 | 20 |

- Polarity is `presence_is_derived` for all records: ALT is derived, H37Rv's
  REF is ancestral, so gains are insertions.
- The 20 losses are scattered, with 1-2 per record.
- 4 records reconstruct derived at the root. These are sparse calls,
  mostly NOCALL or ABSENT; it is a Fitch artefact, not polarity.
- Occupied keys yield no events, as D9 intends.
- `is6110_gene_burden.py` uses `n_gain >= 1`, so occupied keys simply drop
  out of the burden. Losses of H37Rv's copies (2,028 scale200 cells) are not
  represented anywhere; this is D9's documented cost.

**Writer and readers:**

- **Writer:** the state rules match DECISIONS. The derived H37Rv VCF now has
  shared carriers. The duplicate (sample, key) check is fatal, with none in
  the data. `node_occ` is required.
- **Stage 1** (`p5_merge`): skips repeat_node rows (counted), and treats a
  conflicting duplicate as fatal (D10).
- **Merge:** skips repeat_node rows. EVIDENCE and SITECLASS are written as
  `mixed` when carriers disagree (P3IS-6). MEINFO strand is `.`.
- **Stage 2:** the `occupied()` keys match the writer's `h37rv_state`, and
  the order of the NOCALL branch is right (a junction-candidate check comes
  first). `short()` and `arm_key()` agree for node keys with forward offsets.

**Repeat rescue (P3IS-4).**

- `clean_to_orig` and `load_crossmap` are checked against a real crossmap: the
  clean_junction semantics are right.
- A missing crossmap is counted and read as having no excision.

**Reconcile.** ISMapper blank-versus-0 handling is consistent, and `p1i_vcf.sh`
passes `--crossmap-dir` and `--ismapper-dir ""` (P3IS-7).

**`project_sites`:**

- the strand off-by-one correction mirrors P4's;
- `node_occurrences` parses `odgi paths -H` correctly on real output (its
  values are in my run);
- the forward offset is L-1-off for '-'.

**Wrappers.** `p1i_vcf.sh` passes `--node-lengths`, `--refs`, `--h37rv` and
`--build-id`. `p1i_p5states.sh` passes `--h37rv`, `--graph` and `--paths`.

**Accessory:**

- `locus_presence_array.sh`'s sidecar guard sources the env, which defines
  `mtb_resolve_build` and `mtb_kv`. ACC_CATALOGUE is compared byte-for-byte
  with the build copy.
- The merge's caller-INS deduplication sees only non-node IS6110 keys and
  level-1 records emitted earlier in the same shard.
- `node_path_membership.py`: occurrence counts and offsets are per accession,
  with an adjacency check. `--exclude-path` is correct.
- `node_locus_from_p4.py` keys on node id only, so it is unaffected by the
  forward offsets.
- `build_accessory_panel.py`'s `carrier_frac` is now over `--genomes`.
- The `t2_*` and `merge_catalogues` changes only remove defaults.

---

## Audit findings not addressed

None of P3IS-1 to P3IS-10 is left without a fix or a decision. Partial:

- **P3IS-2:** fixed as D6, but with the over-broad rule in R2-IS-1. The
  audit's own open question (are copy_number PRESENT calls real?) is now
  partly answered: for the divergent ones, yes.
- **P3IS-9** (P3 copy-number strand, L+1 window, exit code on missing
  samples): I found no change in the diff. It has no downstream consumer.
  UNVERIFIED whether it was deliberately deferred.

## Open questions

- ACC_3691060 and similar near-identical loci track `reference_carries`
  (113 / 117). The mechanism is unknown: perhaps a tandem-duplication junction.
  Worth one look before D6 is finalised.
- Some samples are called stage-2 REF at an H37Rv key 11-17 bp from their own
  excluded (repeat_node) site. scale200: 76 REF cells at keys within 20 bp of
  such a site, many offset by exactly 11 or 17 bp. These may be distinct
  hotspot insertions, or the same one with the junction just outside the 10 bp
  window. UNVERIFIED; it predates D8.
- R2-IS-3 on gwas1000 is not run.
- The ABSENT counts in IS6110 node records (10,853 scale200 cells) rely on the
  default `--absent unknown`. Any run with `--absent ref` would turn them into
  "no insertion".
