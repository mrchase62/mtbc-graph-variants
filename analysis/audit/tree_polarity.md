# Audit: tree_polarity (alignments, trees, roots, ancestral alleles, event directions)

Date: 2026-10-05. Read-only. Scratch: `/tmp/claude-12043/audit_tree_polarity/`.
WT = working tree `/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome`;
REPO = this repository.

**Which copy runs in production.**
- `assoc/bin/add_outgroup.py`, `assoc/bin/write_event_matrix.py`,
  `bin/vcf_to_alignment.py`, `bin/build_snp_tree.sh` and `bin/panel_polarity.py`
  exist only in WT. The chain (REPO `assoc/bin/cohort_assoc_tail.sh`, steps 1 to 4
  identical to WT's copy) runs them from the WT root.
- `bin/ancestral_alleles.py`: REPO copy, run by `p0_prepare.sh --step ancestral`.
  WT's copy is the older root-state version.
- **The cohort trees actually used** (`data/trees/<cohort>.rooted.nwk`) were not
  made by chain step 3. They came from
  `analysis/combined_tree/build_alignment.py`, then IQ-TREE via
  `bin/build_snp_tree.sh`, then `analysis/combined_tree/prune_for_cohort.py`.
  Those scripts are in REPO and were run by hand (HANDOFF 0b).

**Summary:** 0 HIGH, 6 MEDIUM, 5 LOW, plus the known Faults A and B, both
confirmed.

## Known faults, confirmed

**Fault A.** `add_outgroup.py:120-121`: `exact = g` inside the loop, so the last
duplicate record wins. I recomputed it from the panel VCF (scratch
`og_check.py`). The outgroup is written REF where some exact duplicate has
GT 1:
- scale200: 386 sites;
- gwas1000: 323 sites.

HANDOFF 0j gives 527 and 461. Those were measured against direct alignment, so
they also include errors from other causes; see TP-1.

Re-running the real writer on scale200 with only Fault A corrected (scratch
`ev_s200_ogA`):
- 261 variants change their gain or loss counts;
- 386 root states change;
- 3 variants newly reach 2 or more gains.

**Fault B.** `bin/ancestral_alleles.py:114` takes the root's child that is not
the outgroup. In `data/trees/cx333.rooted.nwk` that node's two children are:
- a 331-tip MTBC clade;
- `GCF_000253375` (the second canettii, edge 0.025).

So the "ingroup" is the canettii plus MTBC node, as described.

---

## TP-1 (MEDIUM): `add_outgroup.py` makes three more outgroup errors besides Fault A

**Where:** WT `assoc/bin/add_outgroup.py:99-133` (production copy).

**What is wrong.** The script reads the decomposed graph VCF, and three cases
are written REF when they should not be.

1. **A missing genotype is read as REF.**
   - Line 118 sets `g = None` for `.`, so `exact` stays `None`.
   - Lines 127-131 then take the `any_rec` branch and write REF
     (`ref_other_record`).
   - The docstring (line 18-19) says "missing -> N".
2. **Records that start upstream are ignored.**
   - `bcftools -R` returns records that overlap the site, but line 107 keys
     them by their own `POS`.
   - So an outgroup MNP, deletion or complex record that starts before the site
     and covers it is never consulted. The site then falls to "no record → REF"
     (line 133). This is bug pattern 1 ("no record means REF") together with
     pattern 6.
3. Where an upstream MNP carries the cohort's ALT base at the site, the
   outgroup is ALT, but it is written REF.

**Evidence.** I re-queried every cohort site, including overlapping records
(scratch `{scale200,gwas1000}.ogx.q`), and classified each written cell against
`data/trees/<cohort>.og.fasta`:

| written REF, should be | scale200 | gwas1000 |
|---|---|---|
| ALT: an exact duplicate has GT 1 (Fault A) | 386 | 323 |
| N: every exact record has a missing GT | 328 | 299 |
| ALT: an upstream outgroup MNP has the ALT base at the site | 157 | 154 |
| N: an upstream outgroup indel or complex record covers the site | 259 | 355 |
| N: an upstream MNP with a third base | 7 | 8 |

That is 751 (scale200) and 816 (gwas1000) wrong cells beyond Fault A. Example:
- site H37Rv 92231 A>G: no exact SNP record exists;
- the outgroup carries `92229 CGA>AAG` with GT 1, which puts G at 92231;
- the file says A.

**Effect on current outputs.** The writer was re-run on scale200 with all four
corrections (scratch `ev_s200_ogall`) and compared with the Fault-A-only run:
- 407 more variants change gain or loss counts;
- 136 roots go from ancestral to derived;
- 271 roots go from ancestral to unresolved;
- 5 variants fall below 2 gains.

Against production (all corrections): 10 rows of `assoc/scale200/scan.tsv`
change gain counts. None had q < 0.05 under all three nulls before or after.
gwas1000 was not re-run, but the cell counts are similar.

**Suggested fix.** Collect all records with `POS <= site <= POS+len(REF)-1`
using a query that keeps the overlap. Then:
- any exact record with GT 1 → ALT (the Fault A fix);
- all exact records missing → N;
- an overlapping MNP with GT 1 → its base at the offset (ALT if it equals ALT,
  otherwise N);
- an overlapping indel or complex record with GT 1 → N;
- only then "no record → REF".

Apply the same rule in `analysis/combined_tree/build_alignment.py` (TP-6).

## TP-2 (MEDIUM): SNP records with a `*` allele get no outgroup and no outgroup polarity

**Where:**
- WT `bin/vcf_to_alignment.py:41` (`if len(rec.alts) != 1 ...: continue`);
- WT `assoc/bin/write_event_matrix.py:620-626` (outgroup lookup keyed on
  `v["alt"]`);
- `write_event_matrix.py:476` (`--panel-polarity` lookup keyed on `v["alt"]`).

**What is wrong.**
- `merge_cohort_vcf.py` appends `*` to ALT whenever any sample is ABSENT, so
  ALT is `C,*`.
- `vcf_to_alignment.py` drops every such record, so the cohort alignment, the
  sites table and `og.fasta` have no column for it.
- The writer then logs "site not in the alignment", and the outgroup tip is UNK
  for that variant.
- The `panel_polarity` lookup also uses `alt = "C,*"` and can never match.
- The toy example shows this. Record `500 T>C,*` has the outgroup in `s.tsv`
  and `og.fasta`, but its outgroup state is 2 (unknown).

**Evidence** (`assoc/<cohort>/events/variants.tsv`, H37Rv-frame SNPs, share
whose root is unresolved):

| | `*` SNPs, root unresolved | other SNPs, root unresolved |
|---|---|---|
| scale200 | 2,592 of 11,163 (23.2%) | 588 of 47,473 (1.2%) |
| gwas1000 | 3,059 of 29,455 (10.4%) | 565 of 104,415 (0.5%) |

An unresolved root means every branch is undetermined, so the variant
contributes no events.

I also re-ran the real writer on the scale200 `*` SNPs (scratch
`ev_s200_star_og`). The `*` was stripped and GT 2 set to `.`, which is
equivalent under `--absent unknown`. The outgroup allele was read from the
panel VCF. Results:
- 2,388 of the 2,592 unresolved roots resolve;
- 46 variants newly have 2 or more gains.

In addition, 1,578 (scale200) and 1,946 (gwas1000) unpolarised `*` records have
a `panel_polarity` row under the stripped key, and that row is never read.
Of those, 163 and 162 say ALT is ancestral.

DR positive controls are almost untouched: 1 unresolved `*` record in pncA
across rpoB, katG, inhA, embB, pncA, gyrA, rpsL, ethA and gid. The loss is
mostly in lineage-level SNPs.

**Suggested fix.**
- In `vcf_to_alignment.py`, treat `X,*` as biallelic, with GT 2 → N. Only
  `build_alignment.py` does this today.
- In the writer, key both lookups on the first non-`*` ALT.

## TP-3 (MEDIUM): the H37Rv tip is called REF at node-frame records, where H37Rv has no sequence

**Where:** WT `assoc/bin/write_event_matrix.py:596-598` and `768-769`.
`g = np.where(is_ref_tip, ord("0"), g)` is applied to every record.

**What is wrong.**
- A node-frame record sits on sequence the H37Rv path does not traverse, so the
  reference tip has no observation there.
- Level-2 conditioning hides this only for the 254 variants inside a placed
  accessory locus. H37Rv is in no presence table, so it is inapplicable there.
- The rest of the node-frame records are not covered:
  - scale200: 7,349 of its 11,740 node-frame records have ALT calls and **no
    REF call in any isolate**;
  - for those, the fabricated H37Rv REF is the only REF in the tree, and it
    creates gains.

**Evidence.** The real writer was run on the node-frame subset twice:
1. as in production; this reproduces production `variants.tsv` exactly, at
   11,740 of 11,740 rows for scale200 and 22,262 of 22,262 for gwas1000;
2. with H37Rv given an all-missing genotype column, so it is UNK.

| | scale200 | gwas1000 |
|---|---|---|
| 2 or more gains only because of the H37Rv REF | 62 (of 371) | 423 (of 2,090) |
| any gain only because of it | 1,869 | 2,798 |

**Effect on current outputs.**
- Event matrices: as above.
- `scan.tsv`: none of these variants reach it. 39 and 141 node-frame
  non-locus rows are in the scans, all unaffected, because the callability
  floor removes the rest.
- Later stages that read the event matrix (burdens at scale, phyoverlap2) would
  see the spurious gains.

**Suggested fix.** Set the ref-sample tip to UNK (bits 0) for records whose
`FRAME != h37rv`. For symbolic SV, IS6110 and accessory_presence records REF is
correct by construction, so keep it there.

## TP-4 (MEDIUM): skip-if-exists guards will hide the Fault A and B fixes at rerun

**Where:**
- REPO `bin/p0_prepare.sh:334`: `_done ancestral`, a marker only;
- REPO `assoc/bin/cohort_assoc_tail.sh:44,52,60,74`: guards on
  `snps.fasta`, `og.fasta`, `rooted.nwk` and `events/summary.txt`;
- REPO `analysis/combined_tree/assoc_tail.sbatch` (prune skipped if
  `rooted.nwk` exists).

**What is wrong.** Every guard tests only that the file exists. None tests
whether it is newer than its inputs or code.
- `refbias/build/7713a8d71d8e/logs/ancestral.done` exists (2026-10-01). After
  Fault B is fixed, the documented command (HANDOFF 4.1,
  `bash bin/p0_prepare.sh --step ancestral`) prints "already done" and keeps
  the faulty `ancestral.tsv`.
- After TP-1 or Fault A is fixed, `cohort_assoc_tail.sh` keeps the old
  `og.fasta` and the old `events/`.
- Nothing is stale today. I checked the timestamps: each product is newer than
  its inputs.

**Effect.** None on current outputs. At the planned rerun, the fixed code would
silently not run.

**Suggested fix.**
- Make each guard compare against its inputs (`-nt`) and record a hash of the
  code.
- At minimum, add marker and output deletion to the HANDOFF 4.1 and 4.2 rerun
  commands: `rm $B/logs/ancestral.done`,
  `data/trees/<c>.og.fasta` and `assoc/<c>/events/summary.txt`.

## TP-5 (MEDIUM): the production tree path is a hand-run step, and its alignment has the add_outgroup faults for all 332 panel genomes

**Where:**
- REPO `analysis/combined_tree/build_alignment.py:113` (`pan[...] = calls`,
  last duplicate wins);
- `:123-125` ("no record → REF", keyed on record `POS`, so upstream records are
  ignored);
- `:129` (missing in a non-exact record treated as acceptable);
- `prune_for_cohort.py`;
- the chain's step 3, `cohort_assoc_tail.sh:59-65`.

**What is wrong.**
1. Bug pattern 4: the trees used in production come from a route the chain does
   not contain. Run on a new cohort, `cohort_assoc_tail.sh` would build a
   cohort-only tree from `og.fasta`, a different kind of tree, without saying
   so. Only the `analysis/combined_tree/*.sbatch` launchers encode the real
   route.
2. `build_alignment.py` reads the decomposed VCF with the same rules as
   `add_outgroup.py`, so every panel genome's column has the TP-1 and Fault A
   errors. HANDOFF 0j measured about 200 wrong cells per panel genome.

**Effect.** Tree topology input for both cohorts. The effect on topology is
UNVERIFIED. `tree_checks.md` reports every lineage monophyletic, so it is
probably small.

**Suggested fix.**
- Move the combined-tree route (alignment, IQ-TREE, prune) into the chain as
  step 3, with the guard fixed (TP-4).
- Share one corrected allele-from-graph-VCF function between `add_outgroup.py`
  and `build_alignment.py`.

## TP-6 (MEDIUM): `--dedupe first` or `drop` misaligns the outgroup alleles

**Where:** WT `assoc/bin/write_event_matrix.py:622-639` and `679-693`.

**What is wrong.**
- `og_gt` is built per index of `variants` before deduplication.
- The `if dup:` block filters `variants` and `gt_rows` but not `og_gt`.
- So every variant after the first removed record reads its neighbour's
  outgroup allele.

**Evidence.** Toy example (`toy/m_dup.vcf`, record 100 duplicated,
`--dedupe first`): the outgroup state vector shifts by one, from
`[0 0 2 0 2 2 1]` to `[0 1 1 2 0 2 2]`. With the shift:
- record 200 reconstructs derived-at-root with a spurious loss;
- record 600 becomes unresolved.

**Effect on current outputs.** None. The chain passes `--dedupe suffix`, after
which no duplicates remain and the block is skipped. It is MEDIUM because the
help text offers `first` and `drop` as normal options, and the error is silent
and affects every record after the first duplicate.

**Suggested fix.** Filter `og_gt` (`og_gt = og_gt[keep]`) together with
`variants` and `gt_rows`, or build `og_gt` after deduplication.

## TP-7 (LOW): panel multi-allelic sites code the other ALT's carriers as REF in the AA alignment

**Where:**
- WT `bin/vcf_decompose.sh:151-152` (`norm -m +any | norm -m -any`, default
  `--multi-overlaps 0`) produces `graphs/.../snps.vcf.gz`;
- WT `bin/vcf_to_alignment.py` reads it into `data/trees/cx333.snps.fasta`,
  which is the input to `ancestral_alleles.py`.

**What is wrong.** When a multi-allelic site is split, a sample carrying ALT2
gets GT 0 in the ALT1 record, not `.`.
- 449 positions are affected.
- 6,587 sample-record cells are 0 in one record and 1 in its sibling.
- No sample has 1 in both records.

The cohort merged VCF does this correctly: the other ALT's carriers are `.`
(scale200: 61 of 61).

**Effect.** I re-ran the REPO `ancestral_alleles.py` with those cells set to N
(scratch `anc_fixed.tsv`):
- 53 of 72,986 AA calls change;
- 52 go from `ref_ancestral` to tied;
- 1 flips to `alt_ancestral`.

**Suggested fix.** Use `bcftools norm -m -any --multi-overlaps .` in
`vcf_decompose.sh`, or have `vcf_to_alignment.py` set N where a sibling record
at the same position has GT 1.

## TP-8 (LOW): `panel_polarity.py` lets the last duplicate win, and the writer's table lookup ignores chrom

**Where:**
- WT `bin/panel_polarity.py:54-70`: writes one row per decomposed record;
- WT `write_event_matrix.py:469-471`: `pol_tab[(pos, ref, alt)] = ...`, the
  last row wins, and the key has no chrom.

**What is wrong.**
- In `refbias/assets/panel_polarity.tsv`, 5,807 keys repeat. At 1,423 of them
  the outgroup is ALT in one record and REF in another:
  - the last row is REF at 903;
  - the last row is ALT at 520.
- `panel_af` is per record, so it is understated for duplicated alleles.
- The chrom-less key could let a node-frame record (POS=1) match an H37Rv
  position-1 row. There is none today: 0 rows at pos 1.

**Effect.** 41 (scale200) and 39 (gwas1000) records are labelled
`ref_ancestral_outgroup` although the outgroup carries ALT. The summed ALT AF
is below `--min-panel-af` for all of them, so the correct label is
`unpolarised`. Derived stays ALT either way, so no events change.

**Suggested fix.** Collapse duplicates in `panel_polarity.py` (any GT 1 → ALT;
`panel_af` = summed AC / AN) and add chrom to the key.

## TP-9 (LOW): optional polarity table is silently skipped when absent

**Where:** WT `write_event_matrix.py:311-312` and `467`. The default path
`refbias/assets/panel_polarity.tsv` is relative, and
`if ... os.path.exists(...)` skips it with no message.

**What is wrong.** This is bug pattern 3. If the writer runs from any other
directory, every outgroup-polarised record silently becomes `unpolarised`.
That is 1,590 records in scale200 and 1,888 in gwas1000, with 78 and 87 of
them flipped to ALT-ancestral. `summary.txt` does not record whether the table
was read.

**Effect.** None today: production shows `*_outgroup` counts.

**Suggested fix.** Fail if the default path is missing unless `--panel-polarity ""`
is passed explicitly, and print the path used in `summary.txt`.

## TP-10 (LOW): hard-coded counts in VCF headers and docstrings will go stale

**Where:**
- REPO `bin/merge_cohort_vcf.py:192-200` (header text: "371 of 72,986 ... tied",
  "6,567 of 72,986 ... 9.0%");
- WT `write_event_matrix.py:28-29` ("8.9%").

**What is wrong.** These are numbers about the current `ancestral.tsv`. After
the Fault B fix, every merged VCF header will state wrong counts. This is bug
pattern 9.

**Suggested fix.** Compute the counts from the table at merge time, or drop
them.

## TP-11 (LOW): no tests cover the event writer or the outgroup, and the AA test misses Fault B's shape

**Where:** REPO `tests/run_tests.py:662-683`. It is the only test in this area.
- Its tree `(O,(A,(B,C)))` has a single outgroup, so it cannot catch Fault B,
  which needs a second outgroup taxon inside the ingroup.
- `add_outgroup.py` and `write_event_matrix.py` live only in WT and have no
  tests.
- The merge test passes `--ancestral ""`, so the AA and AA_INVERTED join is
  untested.

**Suggested fix.** Add the toy cases from scratch `toy/`:
- duplicates;
- missing GT;
- `X,*`;
- node frame plus the reference tip;
- `AA_INVERTED`;
- an expected per-branch event matrix.

Add an AA test whose tree contains two outgroup taxa.

---

## Checked and found sound

- **Fitch in `write_event_matrix.py:244-285`** (binary masks, UNK as wildcard
  upward, downward `parent & own` otherwise own, ties left AMB): checked by hand
  on a 7-variant, 7-tip toy tree run through the real script (`toy/ev`). The
  per-branch matrix is correct:
  - a gain on the expected clade;
  - `AA_INVERTED` gives a gain where H37Rv's clade carries REF;
  - missing leaves give undetermined branches;
  - the root row is excluded from counts (line 948).
- **Event classes:** gain is A→D and loss is D→A. A branch is undetermined
  unless both ends are resolved. No coin-toss tie breaking.
- **Polarity precedence** (`write_event_matrix.py:199-209`):
  `AA_INVERTED` → derived REF; a resolved AA → derived ALT; `AA=.` (tied) →
  unpolarised.
- **The AA join in `merge_cohort_vcf.py:456-458` and `683-695`:**
  - the key is (pos, REF, primary ALT), so `X,*` records do get AA;
  - `AA_INVERTED` only when the flag is empty and AA == ALT;
  - tied or NODATA → `AA=.` plus `AA_FLAG`.
- **AA values:** always one of REF or ALT. The alignment column holds only
  those two bases.
- **CODE_REVIEW 3.7 claim re-verified:** of the 67,856 previously resolved
  sites, 0 changed; 4,759 ties were resolved and 371 remain.
- **Trees:**
  - `cx333`, `scale200`, `gwas1000`, `*_cx333` `.rooted.nwk` all have 2 root
    children, one of them the tip `GCF_035581225`, and no polytomies;
  - the `labelled.nwk` root `n00001` has children `GCF_035581225` and the
    ingroup (scale200: 201 tips; gwas1000: 998);
  - so the trees are rooted where the code assumes.
- **IQ-TREE** (`build_snp_tree.sh`): `GTR+F+ASC+G4`, `-o` equals the reroot
  taxon, and the logs show "0 constant sites" for all three trees, so +ASC is
  valid. `vcf_to_alignment.py` and `build_alignment.py` both drop constant and
  N-only columns.
- **Outgroup tip:** its terminal branch can never carry a gain or loss by
  construction, and production confirms it: scale200 has 0 such events.
  Pruning it is not needed; the chain does not prune.
- **`prune_for_cohort.py`:** keeps exactly the VCF samples, H37Rv and the
  outgroup. The root stays on the outgroup.
- **Cohort merged VCF multi-ALT positions:** the other ALT's carriers are `.`,
  which is correct.
- **`read_vcf` GT parsing:** the first character of the haploid field; `2`
  (`*`) → UNK under the default `--absent unknown`.
- **Ladder collapse runs before `og_gt` is built,** so indices stay aligned.
  Contrast TP-6.
- **Current products are consistent in time:** each is newer than its inputs.
  The subset re-runs reproduce production `variants.tsv` exactly.

## Open questions

1. TP-5: how much do the duplicate and upstream-record errors in
   `build_alignment.py` move the cohort topologies? A tree rebuild is needed to
   know. UNVERIFIED.
2. Fault A count: my strict "any exact duplicate GT 1" count (386 and 323)
   differs from HANDOFF 0j (527 and 461). The 0j figure was measured against
   direct alignment and probably includes TP-1 categories. I did not reconcile
   it site by site.
3. TP-1 "third allele" (28 and 25 sites): an exact record has GT 0 while another
   record at the same POS (often an indel) has GT 1. The script writes REF.
   Whether the base itself is changed depends on the record. Left
   UNVERIFIED and not counted above.
4. TP-2 at gwas1000: only the unresolved-root count (3,059) was measured. The
   number of newly testable variants was not re-run.
5. Indels and SVs never get an outgroup state; the alignment is SNP-only. This
   is by design, but for indels without a `panel_polarity` row the root is
   decided by the cohort alone. That is a design limit, not a bug.
