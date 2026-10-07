# Review 2: build, graph and panel (area `build_graph_panel`)

Reviewer run 2026-10-06 on `audit-fixes` at ffd7bb9 (`git archive` export).
All runs were local and single-core, in `/tmp/claude-12043/review2_build/`.
No cluster jobs were submitted, and nothing was written to production data or
to the working tree.

## Summary

- **A fresh build works end to end through P0 and up to job submission.**
  - **The graph:** test graph `arm_B_L14_L56`, 52 genomes (L1-4, L7, one L5,
    one L6), treated as a new graph.
  - **The guards:** every graph-identity, build-identity and stale-output
    guard accepted the fresh build and its new output folders.
  - **The references:** every P1 and P2 input resolved from the new build.
    P1 selection and the P1 summary ran on its assets.
- **Three things stop or weaken a fresh build:**
  - R2-BUILD-1 (MEDIUM): the non-dry chain is refused by the I/O contract.
  - R2-BUILD-2 (MEDIUM): a build without `ancestral.tsv` is accepted, and
    the failure only surfaces at the last pass.
  - R2-BUILD-4 (LOW): a wrong outgroup is stamped, and rerunning cannot
    correct it.
- **The stale-output guards held in every case I constructed** (section
  "Checked and found sound").
- **The collapse is exactly to specification on CX333.** 93,214 keys match an
  independent union with 0 mismatching cells.
  - Splitting multiallelic records writes the other ALT's carriers as REF,
    because of bcftools' `--multi-overlaps 0` default: 604 SNP keys, 3,849
    cells.
  - Every reader handles this except `panel_polarity.py` (R2-BUILD-3, LOW in
    effect).

No HIGH findings.

## How the fresh build was run (deviations stated)

| step | what ran | deviation |
|---|---|---|
| graph | `arm_B_L14_L56` `.og` and `.gfa`, copied to `graphs/armB/` | built by `analysis/graph_tests/build_arm.sbatch`, not `pggb_build.sh` (cluster job), so there is no `graph_provenance.tsv` |
| deconstruct | the arm's own `deconstruct.vcf` (the same `vg deconstruct -a -P` command), used as `variants.vcf` | not `vg_deconstruct.sh` |
| decompose | branch `vcf_decompose.sh armB --jobs 2`: vcfbub (23,054 → 18,846 records) and vcfwave on all 4 small chunks finished | **The >5 kb chunk (30 records) had not finished after 40 min on this CPU-starved login node.** I killed it, and step 4 was assembled by hand from the 4 waved chunks plus the 30 raw vcfbub records. The >5 kb arm is UNVERIFIED. |
| collapse / split | branch `vcf_collapse.sh armB` and `vcf_split_classes.sh armB`, unmodified | none: 21,935 keys, 0 duplicates; snps 17,743, indels 2,373, svs 1,819 |
| accessory panel | `docs/PANEL_TREE.md` section 4, literally: `t2_extract_candidates` (236 candidates), 51 `t2_validate_accessory.sh` tasks run serially with `SLURM_ARRAY_TASK_ID`, `t2_summary`, `build_accessory_panel` | none |
| P0 | `stamp`, `paths`, `accessions`, `gff`, `refs`, `frames`, `assets`, `catalogue`, `nodes`, `is6110_intervals`, `panel_polarity`, `ancestral`, `manifest`, with `MTB_OUTGROUP=GCF_022870225` (L5.1) | `gff`: 52 GFFs pre-copied from build 7713a8d71d8e (no network fetch). `refs`: 2 refs indexed by the step, then (1 core, GATK about 2 min per ref) 50 copied byte-identical from 7713a8d71d8e; the step then marked done. |
| panel tree | `vcf_to_alignment.py` on `$B/assets/graph_collapsed.vcf.gz` (16,901 sites, 52 taxa), then `build_snp_tree.sh` with local IQ-TREE 3.0.1 (5 min), rooted on L5 | none. Checks: lineage 1 is the first split of the L1-4,7 clade (11 / 39), and `ancestral` with `ANC_OUTGROUPS=L5,L6` gave 734 ALT-ancestral and 3 unknown |
| chain | `refbias_run.sh --dry-run` for a 2-sample scale200 cohort (`SAMEA111556136`, `SAMN09380073`), plus a non-dry run with `sbatch` stubbed | none |

**Checked against CX333:**

- `frames` gave a table identical to CX333's for the same 52 accessions (34
  same frame, 14 rotated, 4 reverse-complemented).
- `is6110_intervals` produced counts for all 51 non-H37Rv genomes.
- P1 selection with the new build's `panel_snps` chose GCF_000786505 and
  GCF_000153685, both in the right lineage. CX333 had chosen GCF_014900175 and
  GCF_977011555, which are not in arm B.
- The P1 summary wrote `refmap.tsv.build` with the new build ID.

---

## Findings

### R2-BUILD-1 (MEDIUM): a fresh build (and a fresh cohort) is refused by the I/O contract, so the chain cannot be submitted without `IO_CHECK=0`

**Where:**

- `refbias/io_contract.tsv:41` (`p3 pre build yes ${BUILD}/assets/accessory_panel.fasta`);
- `:72` (`p1is pre ... is6110/results/${NAME}_p1i_cohort_keys.tsv`);
- `:77` (`p5svgt pre sample ... ${OUTROOT}/p1g/${S}.isclean.bam`);
- `bin/refbias_run.sh:404-419`.

**What is wrong.**

- **`accessory_panel.fasta`:**
  - No P0 step makes it. `bin/p3_accessory.sh:101-122` builds it lazily, on
    the first P3 task, inside the build.
  - CX333's build has it only because P3 already ran there.
  - On any new build it is absent, and `io_contract.py` reports it as a
    missing required pre-input.
- **The other two rows (fresh cohort, not fresh build):**
  - Both files are made by passes in the same chain: p1iv writes
    `is6110/results/<tag>p1i_cohort_keys.tsv` (`bin/p1i_vcf.sh:74`), and p1g
    writes `<OUTROOT>/p1g/<S>.isclean.bam` (`is6110/bin/p1g_isclean.sh:99`).
  - Neither is declared as a `post` of its pass, so the "produced earlier in
    this chain" exemption cannot apply.
  - These two rows are older than this branch, but every rerun into new
    folders now meets them.

**Evidence.** Dry run, then a non-dry run with `sbatch` stubbed, for a new
cohort on the new build (after `ancestral` was done):

```
    3 required INPUT(S) missing -- do not submit:
      p3            -  refbias/build/26bf359e00df/assets/accessory_panel.fasta
      p1is          -  is6110/results/t2new_p1i_cohort_keys.tsv
      p5svgt       0/2  refbias/t2new/p1g/${S}.isclean.bam
FATAL: I/O contract check failed; fix the inputs or set IO_CHECK=0
```

`_stamp_root` had already written `refbias/t2new/.mtb_build` before the refusal.
That is harmless, since it names the same build.

**Effect.**

- The first rerun on the new graph cannot be submitted as documented.
- The only way through is `IO_CHECK=0`. That also switches off the one
  pre-submission check that catches a build with no `ancestral.tsv`
  (R2-BUILD-2).

**Suggested fix:**

- **In P0** (step `assets`, or a small step after it): build
  `assets/accessory_panel.fasta` from `accessory_{novel,mosaic}.fasta`, with
  the bwa index and faidx. The manifest then hashes it, and P3 finds it
  ready.
- **In `io_contract.tsv`:** add the two `post` rows,
  `p1iv post cohort yes is6110/results/${NAME}_p1i_cohort_keys.tsv` and
  `p1g post sample yes ${OUTROOT}/p1g/${S}.isclean.bam`.

### R2-BUILD-2 (MEDIUM): a build without `ancestral.tsv` completes and is accepted by the chain; and `ancestral` run after `manifest` is not verified at use

**Where:** `bin/p0_prepare.sh:697-698` (manifest's list of required steps
omits `ancestral`); `bin/p5_finish.sh:318-321` (the only hard check, at p5vcf).

**What is wrong.**

- `manifest` is marked done once stamp through panel_polarity are done.
  `ancestral` is not required.
- `refbias_run.sh` accepts any build with `logs/manifest.done`. The chain then
  runs P1 to P5 and the IS6110 and SV passes before `p5_finish.sh --merge`
  stops with "no ancestral allele table".
- If `ancestral` is run after `manifest`, as the docs' order allows, then
  `ancestral.tsv` and `ancestral_inputs/` are not in `manifest.tsv`, and
  `p0_check.py verify` never re-hashes them. The AA source of every merged VCF
  can then change without detection until `manifest` is rerun.

**Evidence.** On the fresh build, `--list` showed `ancestral MISSING`,
`manifest done`, and `refbias_run.sh t2new --dry-run` printed "76 manifest
assets verified" and planned every pass. Only `io_contract.py` flagged
`p5vcf - ancestral.tsv`, and R2-BUILD-1 forces it off. After running
`ancestral` and rerunning `manifest`, 4 ancestral rows were added (81 assets).

**Effect:** a whole cohort's compute is spent before a predictable failure.
A late `ancestral` leaves the immutability guarantee broken for the AA table.

**Suggested fix:**

- Add `ancestral` to `step_manifest`'s required list. A graph that genuinely
  has no tree can be given an explicit opt-out, such as `ANCESTRAL=none`
  recorded in `build_info.tsv`.
- Or have `refbias_run.sh` refuse a build without `assets/ancestral.tsv`.

### R2-BUILD-3 (LOW in effect): `panel_polarity.py` reads a sibling-allele carrier's 0 as REF

**Where:** `bin/panel_polarity.py:94-101`.

**What is wrong.**

- At a position with two base-changing alleles (A>G, A>C), each in its own
  record, the outgroup carrying C is GT 0 in the A>G record.
- This happens in two ways. vcfwave writes per-haplotype records, and
  `vcf_collapse.sh`'s `norm -m -any` uses bcftools' default
  `--multi-overlaps 0`.
- The table then states `ancestral = REF` for A>G, but the outgroup has a
  third base.
- `vcf_to_alignment.py` (TP-7) and `add_outgroup.py` both treat this case as
  unknown. `panel_polarity.py` does not.

**Evidence** (`sib.py` in scratch):

- **CX333**, the branch collapse with outgroup GCF_035581225: 260 of 90,364
  keys are REF-ancestral while the outgroup carries a sibling base change
  (176 SNP keys).
- **The fresh build** (outgroup L5): 148 keys (13 SNP). This count is
  inflated by the unwaved >5 kb stand-in records.

**Effect, checked in `assoc/bin/write_event_matrix.py:518-540`.**

- The table is consulted only where AA left the variant unpolarised. "REF
  ancestral" gives `ref_ancestral_outgroup` with derived = ALT, which is the
  same derived allele as the unpolarised default.
- What changes is the label ("polarity is confirmed"). Two options also
  change, and the production tail passes neither:
  - **pinning** under `--root ancestral` (`cohort_assoc_tail.sh:191-198`
    does not pass it);
  - **eligibility for `--polarise-by-root`** (lines 907-911; not passed
    either).

**Suggested fix:**

- **In `panel_polarity.py`:** treat outgroup GT 0 as no call when the outgroup
  carries a sibling allele that changes the same base, using the same rule as
  `vcf_to_alignment.carried_change`.
- **Optionally, also in `vcf_collapse.sh`:** pass `--multi-overlaps .` to
  `norm -m -any`. This covers the multiallelic subset only, not vcfwave's
  per-haplotype siblings.

### R2-BUILD-4 (LOW): the outgroup is stamped unvalidated, and a stamped outgroup cannot be corrected by rerunning

**Where:**

- `bin/p0_prepare.sh:208-243`: `stamp` keys only `build_id` and writes
  `outgroup` with no membership check;
- `:595-599`: `panel_polarity` refuses any mismatch;
- `config/project_env.sh:192`: the default is `GCF_035581225`.

**What is wrong.**

- `docs/PANEL_TREE.md` recommends keeping canettii out of the next graph.
  The config default is still CX333's canettii.
- **A user who runs P0 before setting `MTB_OUTGROUP`:**
  - gets `build_info.tsv: outgroup GCF_035581225` for a graph that does not
    contain it;
  - then sets `MTB_OUTGROUP` and reruns, and gets `stamp: already done`
    followed by `FATAL: build ... was stamped with outgroup GCF_035581225,
    not 'GCF_022870225'`.
- **The fix is not stated anywhere:** the message does not say to delete
  `logs/stamp.done`, and no doc says it either.

**Evidence:** reproduced on a copy of the fresh build. The first
`panel_polarity` failed with "GCF_035581225 is not a sample". After setting
`MTB_OUTGROUP`, `stamp` printed "already done" and `panel_polarity` refused.

**Effect:** a blocked build. The fix is non-obvious but takes minutes. No
wrong output, because `panel_polarity` and the association tail both refuse.

**Suggested fix:**

- Key `outgroup` in `stamp`'s `_inputs`, so a change gives the standard
  "delete the marker" message.
- Check in `accessions` (or `stamp`, after `paths`) that a non-empty outgroup
  is one of the build's accessions.

### R2-BUILD-5 (LOW): the accessory presence folder is adopted with no evidence, against D38

**Where:** `bin/refbias_run.sh:324-327`.

**What is wrong.** D38 says outputs from before the guards are "adopted only
on evidence ...; everything else refused".

- For `accessory/<cohort>`, a folder with files, no `.mtb_build` and no
  stamped VCF is adopted. The non-dry run writes a `.mtb_build` that names the
  new build.
- The per-table check in `accessory/bin/locus_presence_array.sh:67-82` still
  refuses old tables, because they are older than the new catalogue.
- That refusal happens inside the p3acc array. With `afterok`, p5vcf is then
  cancelled at the very end of the chain.

**Evidence:**

```
mkdir accessory/t2new; echo x > accessory/t2new/SAMEA111556136.presence.tsv
refbias_run.sh t2new --dry-run
  accessory/t2new: no build record; adopting build 26bf359e00df (0 stamped outputs checked)
```

**Effect.** This applies when a registry row keeps its cohort name and only
changes `outroot`. The chain then fails late, and the folder's record is
false. No wrong results.

**Suggested fix:**

- For `ACCDIR`, adopt only if every `*.presence.tsv` has a `.build` sidecar
  naming this build, or is newer than the build's catalogue (the same rule
  the array uses).
- Otherwise refuse.

### R2-BUILD-6 (LOW): the graph-identity checks compare genome sets only

**Where:** `bin/p0_check.py:64-73` (`vcf-samples`), `:76-78`, `:81-86`,
`:103-122`.

**What is wrong.**

- Every P0 asset check proves only that a file's sample or taxon set equals
  the build's.
- Two graphs of the same genomes, such as test arms B and D (`-x auto`), or a
  graph rebuilt in place with `MTB_PGGB_RESUME=1`, pass each other's VCFs.
- `bin/vg_deconstruct.sh` now writes
  `##MTB_deconstruct_graph_sha256=<gfa sha>`, which could tie the VCF to the
  graph. It is not checked.

**Evidence:**

```
p0_check.py vcf-samples --vcf analysis/graph_tests/arm_D_L14_L56_sparse/deconstruct.vcf \
    --accessions <arm B build>/assets/accessions.txt --ref-acc GCF_000195955
  samples of ...arm_D.../deconstruct.vcf: 52 names, identical to the build's     rc=0
```

**Effect:** none on the default path, because `GRAPH_VCF` and `PANEL_SNPS`
default to the `.og`'s own directory. It matters only when the VCF is named
by hand, or after an in-place rebuild.

**Suggested fix:**

- In step `assets`, when the VCF header carries
  `##MTB_deconstruct_graph_sha256`, require it to equal the sha256 of the
  `*.smooth.final.gfa` beside `OG`.
- UNVERIFIED: that vcfbub and vcfwave carry the line through. The test arm's
  VCF predates it.

### R2-BUILD-7 (LOW): docs vs code, commands that do not work as written

**`docs/PANEL_BUILD.md:133`**: `bash bin/mash_nearest.sh data/rotated $Q/mash_nearest.tsv`.

- The script globs `*.fna.gz` (`bin/mash_nearest.sh:38-39`), and `data/rotated`
  holds `*.dnaA_rotated.fasta`.
- Run as written: `FATAL: no *.fna.gz under data/rotated`.
- Use `data/assemblies`, or teach the script a suffix.

**`docs/PANEL_TREE.md`, order of sections.**

- Section 1a reads `$B/assets/graph_collapsed.vcf.gz`, which exists only
  after step `assets`.
- Step `assets` needs `ACCESSORY_DIR`, which is made in section 4.
- Read top to bottom, section 1a fails with a missing file. State the order
  as: P0 cheap steps up to `accessions` → section 4 → `assets` → section 1 →
  `ancestral` → `manifest`.

**`INPUTS.md:103-108`** (the P0 command block):

- **It omits `is6110_intervals` and `ancestral`.** The block's
  `--step manifest`, run straight after `frames`, refuses
  (`steps not done: is6110_intervals`).
- **`INPUTS.md:83` describes `<graph>/snps.vcf.gz` as "`vg deconstruct` output
  over the graph".** It is `vcf_split_classes.sh`'s SNP view of
  `all_variants.collapsed.vcf.gz`.
- **That collapsed file is a required P0 input (`GRAPH_VCF`) and is not
  listed.**

**`--step refs` as an array never marks the step done.**

- Only a non-array run calls `_mark refs` (`bin/p0_prepare.sh:360-369`), and
  `frames` refuses without it.
- Every doc goes straight from the array to `--step frames`
  (`p0_prepare.sh:18`, `:721`, `INPUTS.md:106`). The docs should add
  "`bash bin/p0_prepare.sh --step refs` once after the array (it only counts)".
- This was already true on main.

**`refresh_assemblies.sh` without `--summary-only`** overwrites
`data/ncbi/panel_manifest.tsv` (`bin/refresh_assemblies.sh:132`).

- That file is the input of PANEL_BUILD 1h's "CX333 from its recorded
  decisions" reproduction, which then no longer reproduces 333.
- Say so in 1a, or have the script write the new manifest under a dated name.

**`foreign_insertion_screen.py --lineages`** builds its background only from
`--accessions`.

- On a 1-genome run the background is that genome alone, and all 24 inserts
  read FOREIGN on the first pass.
- PANEL_BUILD 1e should say "all candidates in one call", as 1g does.

**Every other PANEL_BUILD and PANEL_TREE command I ran worked as written:**

- 1a `--summary <20260907> --summary-only`: 516 selected, 2 added;
- 1b on 3 genomes: flags GCF_039770655 at 377 bp T;
- 1d on 3 genomes, with live NCBI;
- 1g on 3 genomes;
- 1h with production decisions: 333, identical to the CX333 FASTA's genomes;
- 1i;
- PANEL_TREE 1a, 1b (as `bash`), 1c, 3 (the P0 step) and 4 (all five
  commands).

### R2-BUILD-8 (LOW): panel QC still treats a missing collinearity or barcode row as a pass (PGB-23 partly open)

**Where:** `bin/build_panel.py:242-249` (`c = col.get(acc)`; no row means not
scrambled) and `:222-229`.

**What is wrong.**

- A genome absent from `--barcode` cannot be excluded as mixed or BCG.
- The provenance rule and `--review` now refuse unscreened genomes. These
  two screens still pass them silently.

**Evidence:** none for CX333; it reproduces 333 exactly.

**Effect:** for the next candidate set, a genome missed by
`collinearity_qc.py` (a minimap2 failure drops it with only a stderr line,
PGB-23) or by tb-profiler enters the panel unscreened.

**Suggested fix:** exit if any candidate reaching those rules has no
collinearity or barcode row, as for provenance.

### R2-BUILD-9 (LOW): `make_pansn_fasta.sh` still drops genomes silently (PGB-18 not addressed)

**Where:** `bin/make_pansn_fasta.sh:57-59`.

**What is wrong.**

- A panel genome with no rotated FASTA, or with more than one contig, is
  skipped with a stderr line.
- The FASTA is still written, and the script exits 0. The graph's genomes
  then differ from `build_panel.py`'s list without an error.

**Suggested fix:** exit non-zero when `missing + multi > 0`, after listing them.

### Asides (outside this area, found while running; pre-existing)

**`bin/p4b_place_sv.sh:89,106`.** A missing `delly.vcf` or `dysgu.vcf` is
skipped silently.

- The sample's `sv_placed.tsv` is then written from inherited graph records
  only, and marked done for every later run.
- In my run, with an empty `P2DIR`, it wrote 78 rows and printed `done`.
- Inside the chain, P2 guarantees the inputs, so this matters only for hand
  runs and pinned arms.

**`bin/p1_summary.py:89`.** `--lineages` still defaults to
`data/tbprof_lineages.csv`, CX333's panel labels.

- It is used only for the report's `inside_label` column, so it is harmless.
- A new panel's genomes missing from it print `unknown`.

**`p0_prepare.sh --list` creates the build directory** (`:141` runs before
the `--list` branch).

- `stamp_build_id.sh`'s "no build at all → leave unstamped" branch then sees
  a directory.
- That branch then fails, as it should once any build exists. Noted for
  clarity only.

---

## Checked and found sound

**The guards accept a fresh build.** Verified by running each one:

- `p0_prepare.sh` with `OG` only;
- every step's `.done` marker, on rerun (the full cheap run is idempotent:
  all steps "already done");
- `p0_check` `fasta-names`, `vcf-samples` ×2, `accessory-panel`, `catalogue`,
  `is6110-intervals` and `taxa`;
- `refbias_run.sh`'s build resolution (the one completed build), `verify`
  (80 assets) and frame-table check;
- `_stamp_root` on a new `OUTROOT`;
- the P1 and P2 per-sample preambles: both passed every build check and
  stopped only at a deliberately fake CRAM;
- `p1_summary.py` with the build's `is6110_intervals.tsv` and `--build-id`;
- the P2 refmap `.build` check;
- the P3 preamble: built `accessory_panel.fasta` in the build;
- the P4 preamble.

**The guards refuse stale outputs.** Each case was constructed and confirmed:

- **P0 graph and build identity:**
  - `OG` with another `GRAPH_DIR`;
  - neither set;
  - `BUILD_ID` forced onto a build of another graph.
- **P0 assets:**
  - CX333 `snps.vcf.gz`: 281 not in the build;
  - CX333 collapsed VCF;
  - CX333 `refbias/panel`: no genome record;
  - no `ACCESSORY_DIR`.
- **P0 markers:** `assets` rerun with a changed `MASK` is refused, with
  instructions.
- **P0 `is6110_intervals`:**
  - a crossmap whose clean length plus excised bases differs by 5 bp from
    `refs`;
  - a missing crossmap.
- **P0 `ancestral`:**
  - CX333's tree and alignment: taxa refused;
  - `ANC_OUTGROUPS` unset: the canettii default refused as not leaves.
- **P0 build verify:**
  - an asset whose content changed;
  - an asset replaced by a link out of the build.
- **`stamp_build_id.sh --graph`** on a graph with no build.
- **P1 per sample:**
  - candidates plus an old-stamped H37Rv VCF;
  - candidates with no marker and no VCF;
  - a marker naming 7713a8d71d8e.
- **P1 summary** with a stale marker.
- **P2:**
  - refmap `.build` of the other build;
  - an old-stamped VCF, with no marker and with a date-only marker;
  - a marker with the right build but another reference.
- **`refbias_run.sh`:**
  - an outroot with old-stamped P2 VCFs;
  - an outroot whose `.mtb_build` names another build;
  - an outroot with files but no record and no VCF.

**Collapse (d).** The branch `vcf_collapse.sh` was run on CX333's
`all_variants.decomposed.vcf.gz` (164,650 records; 322 multiallelic split to
164,976).

- **Output:** 93,214 keys, 4 trimmed, 0 duplicate keys.
- **Independent check:** `verify.py` recomputed the union from the raw file
  (1 over 0 over `.`, any allele index counted as called). Result: 0 missing
  keys, 0 extra keys, 0 mismatching cells.
- **Multiallelic splitting writes other-ALT carriers as 0:** 604 SNP keys,
  3,849 cells. Readers:
  - `vcf_to_alignment.py` (N via siblings): sound;
  - `add_outgroup.py` (covering-record rule, third base → N): sound;
  - `t8_select_reference.py`: one difference per differing allele, which is
    correct;
  - P5 `panel_af`: the per-allele frequency is as intended;
  - `snp_nonredundant.py`: counts a two-ALT difference twice. Negligible.
  - `panel_polarity.py`: see R2-BUILD-3.
- **Zero-carrier records:** 1,898 keys with no ALT carrier, 1,687 of them
  already AC=0 in vcfwave's output, pass into the product. Every reader drops
  them or gives them AF 0. Harmless.

**`t8_select_reference.py`** (allele-keyed selection):

- trimming is minimal and unambiguous for SNPs;
- a repeated key is unioned (ALT over REF over missing);
- the isolate side uses only the GT-called allele and skips `*`;
- the panel side assumes GT-only FORMAT, which the collapse enforces.

**`p1_summary.py` tie-break from build intervals.** The counts equal
crossmap rows. A missing count for a chosen or rank-2 reference is fatal.
Zero-copy genomes have header-only crossmaps (4 in production) and give 0,
not absent.

**`p2_call.sh`, `p2_summary.py`:** the marker records build and reference;
the legacy path needs the stamp and the HaplotypeCaller `--reference`;
`P1WORK` is passed through; line terminators are LF.

**The graph-build scripts:**

- **`pggb_build.sh`** (by reading; a cluster job):
  - the emptiness check precedes the sidecar (PGB-1);
  - the defaults equal CX333's `params.yml`;
  - `-x` is explicit;
  - overrides are folded in once;
  - name tokens are checked;
  - `finished` is appended only on success.
- **`vg_deconstruct.sh`:** temp name, `#CHROM` check, refusal beside a
  `variants.vcf.gz`, no `-K`.
- **`vcf_decompose.sh`:** -I 1000 as in production; the big/small partition
  covers every record once (18,816 + 30 = 18,846); collapse is delegated.
- **`vcf_split_classes.sh`:** reads the collapsed file.
- **`sync_back.sh`:** the includes precede the excludes; decomposed,
  collapsed and provenance files are kept.

**Panel scripts:**

- `build_panel.py` reproduces CX333's 333 genomes exactly: 181 excluded, with
  the same reason counts.
- `refresh_assemblies.sh` selects 516, including both M. orygis.
- `assembly_qc_stats.py`, `assembly_provenance_screen.py`, `variant_counts.py`
  and `snp_outlier_screen.py` (all three columns) ran as documented.
- `qc_master_table.py` ran.

**`stamp`** records `outgroup` and `refs_source`, and writes
`panel_provenance absent:...` when the sidecar is missing.

**Test suite:** 254 tests, OK, 0 skipped (bcftools present).

**GFF frame.** On CX333, all 333 GFFs' seqid and region length equal the
staged `refs/` FASTA (P0 does not check this; see open questions).

## Audit findings not addressed (this area)

- **PGB-18:** `make_pansn_fasta.sh` drops genomes silently (R2-BUILD-9).
- **PGB-23:** missing collinearity or barcode rows still pass (R2-BUILD-8).
  `collinearity_qc.py`'s drop-on-failure was not re-checked.
- **PGB-5 / PGB-4 at scale.** The provenance rule's "483 of 484" and the
  foreign screen's validation table were not re-run (CPU-bound here):
  UNVERIFIED. Only the 3-genome and 1-genome runs above were done.
- **PGB-11** (LV meaning after vcfbub): documentation only. Not re-checked.

## Open questions

1. **The >5 kb vcfwave arm** (30 records in arm B) did not finish in 40 min
   on one starved core.
   - `FALLBACK_I` equals the default `WAVE_I` (1000), so the "never wedge"
     comment protects only when `--wave-i` is lowered.
   - Even with inversion disabled (`-I 1e8`), progress was slow, so the cost
     is the WFA alignment itself.
   - Production had 16 cores and finished. Worth a cluster test on arm B, to
     confirm the full script end to end. UNVERIFIED here.
2. **P0 does not check that each downloaded GFF matches the staged reference**
   (seqid and length). Step `gff` takes the latest version from the
   assembly summary, while `refs` decompresses whatever `data/assemblies`
   holds.
   - They agree today on CX333 (333 of 333).
   - A RefSeq version bump between `refresh_assemblies.sh` and P0 would
     misplace annotation silently.
   - A cheap check in step `gff` would close it.
3. **`ANC_OUTGROUPS` and `MTB_OUTGROUP` are independent** (D44). Nothing
   checks that the build outgroup is among the tree's outgroup leaves. Is
   that intended?
