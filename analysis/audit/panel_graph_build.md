# Audit: panel and graph build (`panel_graph_build`)

Audited 2026-10-05, read-only. Scratch: `/tmp/claude-12043/audit_panel/`.

**Scope.** The code that chose CX333 and will build the next panel and graph:

- assembly QC: `assembly_qc_stats.py`, `rotate_to_dnaa.sh`, `collinearity_qc.py`,
  `assembly_provenance_screen.py`, `foreign_insertion_screen.py`,
  `tbprofiler_collect.py`, `snp_outlier_screen.py`, `mash_nearest.sh`,
  `qc_master_table.py`;
- panel and graph: `refresh_assemblies.sh`, `build_panel.py`, `make_pansn_fasta.sh`,
  `snp_nonredundant.py`, `ppe38_build_panel.py`, `pggb_build.sh`,
  `graph_panel_check.sh`;
- variants out of the graph: `vg_deconstruct.sh`, `vcf_decompose.sh`,
  `vcf_collapse.sh`, `vcfwave_rerun.sh`, `pilot_analyze.sh`, `exclude_samples.sh`,
  `stamp_build_id.sh`.

**Which copy runs.** None of these scripts is in the repository's `bin/`, except
`stamp_build_id.sh`, which is byte-identical in both places. So the working-tree
copy (`MtbPangenome/bin/`) is the one that ran and the one the next build would
run. File:line references below are to the working tree unless marked "repo".

**Severity.** For code that only runs at a rebuild, HIGH means "would materially
change the next panel or graph". Two findings, PGB-6 and PGB-8, also affect
current production outputs.

**Counts:** HIGH 5, MEDIUM 10, LOW 9.

---

## PGB-1 (HIGH, next build): `pggb_build.sh` cannot start a clean build

**Location:** `bin/pggb_build.sh:81-96` and `:146-154`. Working tree only.

**What is wrong.** Commit `6c2ec89` (2026-09-29) added a provenance sidecar,
written at lines 81-96 into the output directory: `graph_provenance.tsv` in
`$OUT_HOST`. The guard at line 150 then refuses any non-empty output directory
unless resume is on. Because the sidecar is already there, every fresh build
stops:

```
_prov="${OUT_HOST}/graph_provenance.tsv"     # line 81, written first
...
elif [[ -n "$(ls -A "$OUT_HOST" 2>/dev/null)" ]]; then   # line 150
    echo "[pggb_build] FATAL: ${OUT_HOST} is not empty and resume is off." >&2
    exit 1
```

The only way past the guard is `MTB_PGGB_RESUME=1`. That passes `-r`, the
resume mode that GRAPH_PROVENANCE.md section 4(c) identifies as the reason the
2025 graphs cannot be verified.

**Evidence.** I ran a copy against a stub `singularity` and a stub
environment, in scratch:

- empty output directory: `FATAL ... is not empty`, exit 1;
- `MTB_PGGB_RESUME=1`: pggb is called with `... -k 51 -K 21 --keep-temp-files -r`.

**Effect on current outputs:** none. CX333 was built on 2026-09-09, before the
sidecar existed; its directory has no `graph_provenance.tsv`. The 0l test
graphs bypass this script (`analysis/graph_tests/build_arm.sbatch` calls pggb
directly).

**Suggested fix:** move the guard (lines 146-154) above the sidecar write, or
write the sidecar into the directory after pggb finishes.

## PGB-2 (HIGH, next build): `pggb_build.sh` defaults are not CX333's settings, and RUNBOOK passes none

**Location:** `bin/pggb_build.sh:26` and `:130-134` (`MIN_MATCH_LEN=51`,
`MASH_KMER=21`); `#SBATCH -c 112` and `-t ${SLURM_CPUS_PER_TASK:-112}` (lines 4
and 162). RUNBOOK.md, Segment 2.

**What is wrong.** The settings disagree across the files:

- **Script defaults:** `-s 10000 -l 30000 -p 95 -k 51 -K 21`, on 112 cores.
- **What CX333 was built with:** `-k 23 -K 15` on 48 threads. This is recorded in
  `graphs/CX333.s10k.k23.K15/*.params.yml` (`min-match-len: 23`,
  `mash-kmer: 15`, `threads: 48`) and in `data/qc/cx333_build.cmd`.
- **What QC_PIPELINE.md stage 3 documents:** the same `-k 23 -K 15`, and
  "Use 48 cores, not 112".
- **RUNBOOK.md Segment 2** runs `sbatch -p intermediate -t 14-00:00
  bin/pggb_build.sh data/fastas/mtb.complex333.fasta.gz CX333.s10k.k23.K15`
  with no pggb arguments. Following it gives a graph with a different seqwish
  match filter and a different mash k-mer, under a directory named `k23.K15`.
- **The RUNBOOK input path is also wrong:** the script prefixes `$MTB_FASTAS`
  (line 58), so `data/fastas/...` resolves to `data/fastas/data/fastas/...`.

**Evidence:** the stub run in PGB-1 shows the arguments pggb would receive.

**Effect:** none now. At the rebuild, the graph would silently differ from
CX333 and from the 0l test arms, which used `-k 23 -K 15`.

**Suggested fix:**

- set the defaults to `MIN_MATCH_LEN=23` and `MASH_KMER=15`, the documented
  settings;
- fail when the output directory name implies settings (`k23.K15`) that
  disagree with the arguments;
- fix the RUNBOOK command: bare FASTA name, explicit arguments, `-c 48`;
- record the container digest in the provenance file, not just
  `pggb_latest.sif`.

## PGB-3 (HIGH, next build): no code produces the CX333 panel; the documented command overwrites it

**Location:** `bin/build_panel.py:66-70`, `:135-142`; RUNBOOK.md Segment 1, step
1i; `data/qc/cx333_build.cmd`.

**What is wrong.**

1. **`build_panel.py` implements only 514 to 484.** It has no input for the
   provenance screen, the foreign screen or the SNP screen. The 484 to 333 step
   (the 150 provenance exclusions, the two engineered genomes, and H37Rv kept as
   the reference path) exists only as a comment in `cx333_build.cmd`.
   - No script writes `mtb.complex333.fasta.gz`. `make_pansn_fasta.sh` defaults
     to `panel.rebuild.tsv` (484) and the name `mtb.complex490`.
   - `qc_master_table.py --panel-fasta` then takes the built FASTA as the
     authority for `in_panel` (lines 95-101). The table records the decision
     after the fact; it does not make it.
2. **`build_panel.py` writes a 3-column TSV** (lines 135-138). RUNBOOK 1i passes
   `--out data/fastas/mtb.complex333.fasta.gz`, so running the runbook would
   overwrite the production panel FASTA with a TSV.
3. **RUNBOOK 1i passes `--manifest data/qc/qc_master_table.tsv`.** The script
   reads every line, header included (lines 67-70), and takes column 3 as the
   organism. In that table, column 3 is `exclusion_reason`.
4. **RUNBOOK 1i omits `--barcode` and `--manual-exclusions`.** Without them, the
   mixed-lineage and manual SNP-excess rules silently do nothing.

**Evidence.** I ran it with outputs in scratch.

- **With the RUNBOOK arguments:** `515 selected -> 492 kept`.
  - the header row `accession` is "kept" as a genome;
  - GCF_000738445 (mixed lineage) and all 6 `manual_exclusions.tsv` genomes are
    kept;
  - BCG is excluded only because column 3 holds the earlier reason
    `BCG_by_barcode`;
  - the output file `rb_out.fasta.gz` is a 493-line TSV.
- **With `panel_manifest.tsv`, `--barcode` and `--manual-exclusions`:** 484 kept,
  identical to `data/ncbi/panel.rebuild.tsv` (diff empty), and it contains all
  333 CX333 genomes. So the code reproduces the first round, not the panel.

**Effect:** none on current outputs. At the rebuild, the panel cannot be
regenerated from code, and the documented command destroys the panel FASTA.

**Suggested fix:**

- make `build_panel.py` read every screen's table (provenance, foreign,
  long-read checks, SNP/indel outliers, retained exceptions);
- have it write the PanSN FASTA itself, or call `make_pansn_fasta.sh` with its
  own output list;
- skip the manifest header;
- make `--barcode` required;
- fix RUNBOOK 1i.

## PGB-4 (HIGH, next build): the provenance screen in code is not the rule that removed 150 genomes

**Location:** `bin/assembly_provenance_screen.py:153-172`;
`data/qc/assembly_provenance.tsv`; QC_PIPELINE.md section 1.3; RUNBOOK 1d.

**What is wrong.** The production table was not written by this script.

- **The columns differ.** Production has `is6110_total, is6110_ref_kept,
  is6110_novel, identical_profile_94, flag`. The code writes `inserted_bp,
  deleted_bp, insertions_ge50, verdict, reasons`.
- **The decision rule differs.** In production:
  - `REFERENCE_STRUCTURE` is exactly `identical_profile_94 = yes` (94 of 94);
  - `SHORT_READ_SUSPECT` is exactly `tech_class = short_read` (56 of 56);
  - `unknown` technology passes (17 genomes, all `ok`).

  The code instead excludes on fewer than 5 insertions of 50 bp or more, plus
  `short_read_only` only when `--require-long-read` is given, plus assembler
  names.
- **The documented threshold was never stored.** The production table has no
  insertion count, so QC_PIPELINE.md's ">= 5 insertions excluded 118 of 118" is
  not the recorded basis for CX333.
- **RUNBOOK 1d omits `--require-long-read`.** The runbook command would
  therefore keep the 56 short-read-only genomes.

**Evidence.** Counts from production `assembly_provenance.tsv`, flag against
`tech_class` and `identical_profile_94`:

| profile identical to H37Rv | technology | flag | n |
|---|---|---|---:|
| yes | short_read or unknown | REFERENCE_STRUCTURE | 94 |
| no | short_read | SHORT_READ_SUSPECT | 56 |
| no | long_read / hybrid / unknown | ok | 334 |

I imported the module's `tech_class()` and checked it against production. It
reproduces production's `tech_class` for all 484 genomes. So the technology
part matches, and the sequence part does not.

**Effect:** none now. At the rebuild, running the code gives a different
excluded set. How different is UNVERIFIED: I did not rerun the screen, which is
about 8 CPU-hours of minimap2.

**Suggested fix:**

- choose one rule and record it in the code, with the insertion count, the
  IS6110 profile and the technology all written per genome;
- make `--require-long-read` (or its replacement, PGB-7) the default;
- make the table the input to `build_panel.py` (PGB-3).

## PGB-5 (HIGH, next build): the foreign-DNA screen can pass engineered genomes

**Location:** `bin/foreign_insertion_screen.py:103-119`, `:77`, `:94`, `:28-48`.

**What is wrong.** There are six problems; (b) and (e) are verified in scratch.

- **(a) The code sums aligned bases over all non-self hits** (up to 50
  secondaries, lines 109-119). The production table's column is
  `best_homologue_bp`, a single best homologue. QC_PIPELINE.md 1.4 documents
  `-N 100 -p 0.05`; the code uses `-N 50 -p 0.1`. A short native-looking piece
  that hits many genomes can push `frac` above 1.
- **(b) Two engineered genomes that share a vector hide each other.** Each
  one's insert finds the other genome as a "homologue".
- **(c) RUNBOOK 1e passes the whole panel as the background.** This is the known
  disagreement with QC_PIPELINE.md 1.4.
- **(d) Self-hits are excluded by the PanSN prefix** (`tgt.split("#")[0]`,
  line 116). A background built from bare contig names would count every
  self-hit, so every insert from a background genome would read as native. The
  recommended one-per-sublineage background must therefore be PanSN-named.
- **(e) Inserts are taken only from `paftools call` indels inside one
  alignment.** An artifact that breaks the alignment is never extracted
  (PGB-6).
- **(f) Skip-if-exists on `inserts.fa` and the PAF** (lines 77 and 94). A partial
  or stale file from an earlier accession list or `--min-len` is reused.

**Evidence for (b).** I ran the script, unchanged, with GCF_044324775 (pJEB) and
GCF_021535155 (attB vector) as the background (`fs/foreign2.tsv`):

| genome | insert at H37Rv 2,765,569 | aligned elsewhere | verdict |
|---|---|---:|---|
| GCF_044324775 | 5,214 bp | 3,493 | **native** (frac 0.67) |
| GCF_021535155 | 4,445 bp | 3,474 | **native** (frac 0.78) |

Production called both FOREIGN, with best homologue 0 and 379, because its
66-genome background held only one of the two.

**Evidence for (e).** GCF_039770655's 1,144 bp artifact (PGB-6) is absent from
`data/qc/foreign_insertions.tsv`. Its 40 kb window aligns to H37Rv as separate
blocks, ending at 2,001,528 and resuming at 2,002,483, so no indel is called.

**Effect:** none on the current panel; both engineered genomes are excluded by
`engineered_exclusions.tsv`. At the rebuild, any two genomes carrying the same
construct, or an artifact that breaks the alignment, would pass.

**Suggested fix:**

- score `best single-genome homologue / insert length`, as production did;
- also align inserts to a vector database (UniVec) and to IS6110, the
  `foreign_recheck.py` approach;
- extract unaligned query gaps between consecutive alignment blocks, not only
  `paftools` indels;
- require a PanSN background, or exclude self-hits by a background-to-accession
  map;
- put a hash of the inputs in the workdir and refuse a mismatch;
- exempt the outgroup in code.

## PGB-6 (MEDIUM, current production): a third artifact genome, the lineage 9 reference, puts false events in both cohorts

**Location:** the panel FASTA, `GCF_039770655#1#NZ_CM129972.1:2,070,500-2,071,700`.
GCF_039770655 is ITM-2021-00504, lineage 9, Flye/Medaka ONT. No screen in
`bin/` tests for this (PGB-5(e) and PGB-7).

**What is wrong.**

- **The artifact.** The genome carries about 490 bp of `(CCATT)n`, a 377 bp pure
  poly-T, then more `(CCATT)n`. Its GC is 29.5% against 65% for the genome.
  It replaces 955 bp of the 3′ end of PE_PGRS31 (Rv1768): H37Rv 2,001,528 to
  2,002,483 becomes 1,144 bp of artifact.
- **Not real lineage 9 sequence.** The closest panel genomes (GCF_022870165,
  GCF_965124535) carry the H37Rv sequence intact: 1,201 of 1,201 bp aligned.
- **A new finding.** HANDOFF 0f lists only GCF_045348265 and GCF_050259585 as
  artifact genomes.

**Evidence.** GCF_039770655 is the matched reference for every lineage 9
isolate: 14 in gwas1000 and 2 in scale200 (`refbias/*/p1/refmap.tsv`). The
merged VCFs carry the artifact as inherited calls:

| record | gwas1000 ALT / REF | scale200 ALT / REF |
|---|---|---|
| `h37rv:2001528:T>TGGCTGG...CATTCC...` (1,144 bp insertion, `CLASS=small;KIND=INDEL`) | 14 / 871 | 2 / 176 |
| `sv:INS:2001528:1143` (inherited) | 14 / 0 | 2 / 0 |
| `h37rv:2001529` 955 bp deletion of PE_PGRS31 | 14 / 876 | 2 / 178 |
| `acc:ACC_2001528` (accessory presence of the 1,143 bp) | 9 / 728 | 5 / 135 |

**Effect.** About 4 to 5 false records per cohort.

- They are carried by exactly the lineage 9 isolates, so they are perfectly
  confounded with lineage 9 in any association or event table.
- The accessory-presence ALT calls (9 and 5 isolates) suggest reads matching the
  low-complexity repeat. Whether that is human-satellite contamination in the
  reads is UNVERIFIED.

**Suggested fix:**

- **at the rebuild:** mask the region with N or replace the assembly;
- **for the current data:** add the three records, and `ACC_2001528`, to a known
  artifact blacklist used by the chain. Treat this like the two 0f genomes under
  the one-rerun rule.

## PGB-7 (MEDIUM, next build): no screen tests long-read error, low-complexity runs or technology from sequence

**Location:**

- `bin/assembly_qc_stats.py:33-42`: only N runs are measured;
- `bin/assembly_provenance_screen.py:48-63`, `:166`: metadata only.

**What is wrong.**

**(a) No single-base-run check.** I scanned all 333 panel sequences
(`scan333.tsv`). Exactly three carry a pure single-base run of 100 bp or more,
and all are artifacts:

| genome | longest run | base | panel position |
|---|---:|---|---:|
| GCF_050259585 | 1,015 | G | 1,585,328 |
| GCF_045348265 | 484 | A | 546,697 |
| GCF_039770655 | 377 | T | 2,070,991 |

The other 330 genomes have no run of 50 bp or more. GCF_045348265's 1,243 bp
poly-A is three pure runs (484 bp longest), so the rule should be windowed
rather than a strict homopolymer test. For example: at least 90% one base over
200 bp, or GC below 35% over 500 bp.

**(b) Technology classification is metadata-only and mis-parses common strings.**
I called `tech_class()` from the module:

| technology string | class | consequence |
|---|---|---|
| `Illumina NovaSeq; ONT` | `short_read` | a hybrid would be excluded |
| `HiSeq 2500` | `unknown` | passes, because only `short_read` is excluded |
| `MiSeq/PacBio RSII` (GCF_000195835, in panel) | `long_read` | should be `hybrid` |

**(c) Assembler names are not checked.** `ALIGNER` lists mapping tools only.
Two "long_read" panel genomes name SOAPdenovo: GCF_050259585 (SY-1, the poly-G
genome) and GCF_000195835.

**(d) Unknown technology passes.** 17 genomes passed with no technology
recorded.

**(e) The read-free long-read checks are not in the pipeline.** The indel
excess and private homopolymer-indel excess checks from 0f exist only in
`analysis/external_assemblies/bin/`.

**Effect.** At the rebuild, error-rich and short-read genomes enter as they did
in CX333. 0h found 82 of 173 GenBank candidates flagged.

**Suggested fix:**

- add a windowed low-complexity and single-base-run rule to
  `assembly_qc_stats.py`;
- port `variant_counts.py`, `repeat_checks.py` and `hybrid_vs_sr.py` into `bin/`
  and the chain;
- classify technology from the sequence (insertion count and IS6110 novel
  copies) and treat metadata as supporting evidence only;
- add `ont`, `miseq`, `hiseq`, `nextseq`, `novaseq`, `iseq` and `dnbseq` to the
  keyword lists;
- add short-read assemblers (SOAPdenovo, SPAdes-only, Velvet, ABySS, SKESA,
  Shovill) as a conflict flag;
- review `unknown` technology by hand, rather than passing it.

## PGB-8 (MEDIUM, current panel): a quarter of CX333 never had the SNP-outlier screen, and its input has no producer

**Location:** `bin/snp_outlier_screen.py`; `data/qc/snp_outliers.tsv`;
`bin/qc_master_table.py:38-40`; RUNBOOK 1g (`--counts data/qc/snp_counts.tsv`).

**What is wrong.**

- **The counts come from the 2025 412-genome graph** (`qc_master_table.py`
  docstring).
- **79 of the 333 panel genomes (24%) have no row** in `snp_outliers.tsv`.
  GCF_045348265 is among them.
- **`snp_counts.tsv` does not exist,** and no script in `bin/` writes it.
- **It is the only error screen beyond provenance,** and it counts SNPs, not
  indels.
- **Groups with MAD 0 get z = 0 and can never flag.** Only one group (lineage 8,
  n = 1) is affected now.

**Evidence:** the set difference between `p333.txt` and `snp_outliers.tsv`
gives 79. `qc_master_table.tsv` has `snps` blank for those 79 panel rows.

**Effect.** 79 panel genomes were never screened for SNP excess.
GCF_045348265's artifact is not SNP-visible, so no wrong exclusion is known.

**Suggested fix:**

- produce the counts in code, from direct alignment to H37Rv (the
  `analysis/panel_checks` states method), not from a graph that does not exist
  yet;
- add indel and homopolymer-indel counts;
- fall back to the major lineage when MAD is 0.

## PGB-9 (MEDIUM): `vcf_decompose.sh` publishes the duplicate-record file that caused fault A; the collapsed file already fixes it

**Location:** `bin/vcf_decompose.sh:130-155`; `bin/vcf_collapse.sh`;
`bin/sync_back.sh:102-105`.

**What decomposition does to records** (CX333, measured from
`all_variants.{decomposed,collapsed}.vcf.gz`):

- `vg deconstruct -a` gives 86,186 records, by level: LV0 15,273, LV1 55,485,
  LV2 12,825, LV3 to LV7 2,602.
- `vcfbub -l 0 -a 100000` gives 64,393 (LV0 15,273, LV1 49,120). No LV1 record's
  parent is present, so nothing is represented twice (see PGB-11).
- `vcfwave -I 1000` gives **164,650 records but only 93,214 distinct
  (pos, ref, alt).**
  - 6,186 keys occur in more than one record, with 71,762 extra records.
  - In 6,129 of those keys, carriers are split across records.
  - This is fault A's mechanism. It comes from vcfwave decomposing each ALT
    allele separately, not from the deconstruct settings. Deconstructing
    without `-a` would not prevent it.
- **The collapse already in the script** (`norm -m +any | norm -m -any |
  fill-tags`) gives 93,214 records, one per key.
  - After allele trimming, carriers are identical to the union over the
    decomposed records, except **64 carrier cells at 12 keys**. These are
    samples carrying two different overlapping alleles at one POS.
  - **41 events remain duplicated** after left-alignment
    (`bcftools norm -f GCF_000195955.pansn.fasta`).
  - **8,934 records are not left-aligned.** The script's comment says
    "canonicalise", but no `norm -f` is run.
  - **LV and ORIGIN in a collapsed record come from one arbitrary input record.**
- 1,675 decomposed records (1,910 after collapse) have no ALT carrier. They are
  dropped without comment by `exclude_samples.sh:72` (`view -c 1`), so 93,214
  becomes 91,304 even though 0 of 6 lab samples were present (job 45740274
  log).

**What is wrong.** Both files are written, and the script says to use the
collapsed one (line 159). Production readers still read the decomposed file:

- `assoc/bin/add_outgroup.py`;
- `bin/sv_intervals.py` (default);
- repo `bin/p5_svgt.sh:174` (default);
- `bin/panel_polarity.py`, `accessory/bin/merge_catalogues.py`,
  `assoc/bin/caller_graph_overlap.py`, `insgt/bin/insertion_contigs.py`,
  `sv2frame/bin/clip_discovery.py`;
- six `analysis/` scripts.

`sync_back.sh` excludes the decomposed file from the durable mirror as
"superseded", even though production depends on it.

**Effect.** Each reader of the decomposed file must handle the duplicates
itself (see the graphvcf audit for each reader).

**Suggested fix:**

- in `vcf_decompose.sh`, write only the collapsed, left-aligned file
  (`bcftools norm -f <H37Rv PanSN> -m +any`, then `-m -any`);
- keep the per-allele file under a name no reader will pick up, such as
  `decompose/waved.vcf.gz`;
- write zero-carrier records and the 64 conflict cells to a side report;
- point every reader at the one file;
- stop excluding it in `sync_back.sh` while anything reads it.

## PGB-10 (MEDIUM, next build): the decomposition code is not what produced the production VCF

**Location:** `bin/vcf_decompose.sh:69` (`WAVE_I=64`); `bin/vcfwave_rerun.sh:16-47`;
log `slurm/analyze_an_CX333_45740274.out`.

**What is wrong.**

- **Production used `-I 1000`.** The CX333 run (2026-09-11) logged
  `vcfwave -I 1000 across 8 parallel jobs`.
- **The code changed afterwards.** Commit `4d46708` (2026-09-13) changed the
  default to 64. The `-I 64` result (`all_variants.inv64.*`) exists but is not
  the file P4 or P5 read (`all_variants.nolab.vcf.gz`, from `-I 1000`).
- **The current script has no fallback for long alleles.**
  `vcfwave_rerun.sh:40-47` records a 23-minute stall on one multi-kilobase
  allele at `-I 64`, and partitions records over 5 kb back to `-I 1000`.
  `vcf_decompose.sh` does not, so a rebuild can hang on the 12-hour limit.
- **No decomposition setting is written into any VCF header.**

**Effect:** none on current files. Rerunning from code gives a different call
set from production, or stalls.

**Suggested fix:**

- move the large-allele partition into `vcf_decompose.sh`;
- decide on `-I` once;
- write `##MTB_decompose=vcfbub -l 0 -a 100000; vcfwave -I N ...` and the vg,
  vcfbub and vcfwave versions into the header.

## PGB-11 (MEDIUM): after vcfbub, LV no longer means nesting depth; one 3.4 Mb snarl explains "most SNPs are nested"

**Location:** `bin/vcf_decompose.sh:63` (comment) and `:100-101`.

**What is wrong.**

- **One top-level snarl covers 77% of the genome.** The deconstruct VCF has a
  single LV0 record at H37Rv 415,369 (`>28744>222850`), with REF 3,388,550 bp
  and an allele of 3,478,275 bp. It spans 415,369 to 3,803,919.
- **`vcfbub -a 100000` drops it and promotes its children.** All 49,120 LV1
  records in `bub.vcf` (76%) have that snarl as their `PS`. They keep `LV=1`
  although nothing in the file nests them.
- **The comment says the opposite.** The script describes `-a` as "alleles
  longer than this stay collapsed"; in fact, those sites are removed.

**Effect:**

- the rearrangement itself is not represented in any variant file;
- any reader that keeps LV0 only, or pools by snarl, loses three quarters of the
  genome. This is the root of the "nested snarls hold most SNPs" misreads in
  pattern 1;
- a lineage 1-4 graph probably will not have this snarl, so LV distributions
  will shift between builds;
- the snarl's identity is UNVERIFIED. It probably reflects the B0/W148 nested
  rearrangement (704,251 to 3,711,742) or another long antiparallel traversal.

**Suggested fix:**

- after vcfbub, rewrite promoted children to `LV=0`, or drop `LV`/`PS` and add
  `POPPED_PARENT=<id>`;
- write the popped >100 kb sites to a side table so rearrangements are kept;
- correct the comment.

## PGB-12 (MEDIUM, next build): `snp_nonredundant.py` gives inflated distances on the decomposed VCF, and its default tie-break differs from its docstring

**Location:** `bin/snp_nonredundant.py:54-62`, `:44-48`, `:116-123`.

**What is wrong.**

- **It accepts any VCF.** With the decomposed file, one ALT split across
  duplicate records counts as several differences.
- **Multi-allelic genotypes are flattened:** GT 1 and GT 2 both become 1.
- **H37Rv is dropped.** The reference path is never a sample in a graph VCF, so
  GCF_000195955 is dropped with a warning, and the H37Rv cluster is formed
  without it.
- **The default representative is the lexicographically largest accession**
  when `--rank-by` is not given. The docstring promises "most insertions
  >= 50 bp".

**Evidence:** 10 CX333 genomes (3 B0/W148 genomes, GCF_026185275 and 6 lineage
4.4.1.1 genomes), SNP records only:

| input | sites | pairwise min / median / max | clusters at <= 50 SNPs |
|---|---:|---|---:|
| decomposed | 128,859 | 16 / 868 / 1,636 | 7 |
| collapsed | 77,056 | 3 / 797 / 1,428 | 4 |

**Effect:** none recorded; the input of the one existing output
(`panel333_lineage2.2.1.nonredundant.tsv`) is not documented. At the rebuild,
"collapse clones within about 50 SNPs" (0i step 3) would under-collapse if fed
the decomposed file.

**Suggested fix:**

- compute distances from direct alignment to H37Rv for all candidates, so the
  reference and non-panel candidates are included;
- otherwise, reject input whose (pos, ref, alt) keys are not unique, and
  require `--rank-by`.

## PGB-13 (MEDIUM, next build): the RefSeq selection regex excludes M. orygis, and GenBank is not searched

**Location:** `bin/refresh_assemblies.sh:57-59`.

**What is wrong.**

- **The organism regex misses M. orygis.** It is
  `^Mycobacterium tuberculosis|^Mycobacterium canetti`. The current RefSeq
  summary has two complete, latest genomes named `Mycobacterium orygis`:
  GCF_015265495 and GCF_033782915. The regex excludes them. This, not
  availability, is why CX333 has no M. orygis (0g lists it as a gap). It is
  the same class of error as the earlier "canettii" spelling.
- **GenBank-only genomes are never considered.** Section 0h found 173 complete
  GenBank genomes not in CX333.

**Suggested fix:**

- select by taxid (MTBC 77643 and its descendants) instead of by name;
- add `assembly_summary_genbank.txt`, de-duplicated against RefSeq by
  BioSample.

## PGB-14 (MEDIUM, next build): `stamp_build_id.sh` would stamp a new graph's VCF with the old build

**Location:** `bin/stamp_build_id.sh:57-71` (repo copy identical).

**What is wrong.**

- **The build directory is not matched to the graph.** The script resolves "the
  only directory under `refbias/build`" and stamps that build's id and graph
  hash into any VCF. It does not check that the VCF came from that graph.
- **RUNBOOK Segment 3 stamps right after `vg deconstruct`,** before a new P0
  exists. A new graph's VCF would therefore be labelled `7713a8d71d8e` (CX333).
- **With no build directory, the script exits 0** and leaves the file unstamped.

**Suggested fix:**

- take the graph path as an argument and refuse when its sha256 differs from
  `build_info.tsv`'s `graph_sha256`;
- exit non-zero when nothing was stamped.

## PGB-15 (MEDIUM, next build): P0 references are not the sequences the graph was built from

**Location:** repo `bin/p0_prepare.sh:51` and `:209-217`, which stage
`data/assemblies/<acc>.fna.gz`. The panel FASTA holds the dnaA-rotated sequences.

**What is wrong.** The two frames differ for 110 of 333 accessions: shifted
origin, 22 of them also on the opposite strand
(`graphframe/docs/GRAPH_FRAME_RESOLUTION.md`). A whole conversion layer
(`graphframe/`) and at least one HIGH bug (CODE_REVIEW 3.2, the node-frame
strand) followed from this.

**Suggested fix:** at the rebuild, stage P0 references by
`samtools faidx <panel.fasta.gz> <acc>#1#<contig>`, so every reference is
byte-identical to its graph path. The frame table then becomes the identity, and
it can be asserted.

## PGB-16 (LOW): mixed-lineage screen cannot see a mixture within one lineage

**Location:** `bin/tbprofiler_collect.py:50-57`.

**What is wrong.** Calls are grouped by root (`lineage4`, `La1`), and the
longest call wins. A lineage 4.1 plus 4.3 chimera has one root and is not
flagged, and most panel genomes are lineage 4.

**Suggested fix:** flag a genome when two calls at or above `--min-fraction`
are not prefixes of each other.

## PGB-17 (LOW): skip-if-exists guards accept partial outputs

**Location:**

- `bin/pilot_analyze.sh:41-48`: `-s variants.vcf` and
  `-s all_variants.decomposed.vcf.gz`;
- `bin/vg_deconstruct.sh:61-64`: writes straight to the final name, under a
  4-hour limit;
- `bin/vcf_decompose.sh:141`: `bcftools sort -o` straight to the final name;
- `bin/rotate_to_dnaa.sh:60`: skips if the FASTA and the log exist, even when
  the log says FAILED.

**What is wrong.** A job killed mid-write leaves a non-empty file, and the next
run skips the step.

**Suggested fix:** write to `*.tmp`, rename on success, and test a completion
marker or `bcftools index` success.

## PGB-18 (LOW): failed rotations are kept, and PanSN assembly drops genomes without a record

**Location:**

- `bin/rotate_to_dnaa.sh:91-104`: `REVERSE_STRAND`, `OFFSET_*` and
  `FAILED_dnaA_not_found` still copy the output to the final name and exit 0;
- `bin/make_pansn_fasta.sh:43-52`: does not read the rotation status, and skips
  missing or multi-contig genomes with only a stderr line.

**Current state:** all 333 panel genomes are `OK +` in `rotation_verify.tsv`.
One non-panel genome is `OFFSET_4397242`.

**Suggested fix:**

- exit non-zero, or write to a quarantine directory, unless the status is OK;
- have `make_pansn_fasta.sh` require OK and fail on any skip.

## PGB-19 (LOW): IUPAC masking is documented but not done

**Location:**

- ASSEMBLY_QC.md section 8 says "`bin/build_panel.py` applies these ... mask to
  N";
- `build_panel.py` writes a list only;
- `make_pansn_fasta.sh:49-50` passes the sequence through untouched.

**Evidence:** `GCF_000572175#1#NZ_CP002883.1` carries 20 IUPAC bases in the
panel FASTA (K 4, M 2, R 7, S 1, W 1, Y 5). No other panel genome has any.

**Suggested fix:** mask non-ACGTN to N in `make_pansn_fasta.sh`. Run the
low-complexity masking from PGB-6 and PGB-7 at the same point.

## PGB-20 (LOW): other RUNBOOK Segment 1 commands disagree with the code

**Location:** RUNBOOK.md lines 59-100.

**What is wrong:**

- **1b:** `assembly_qc_stats.py --dir data/rotated` globs `*.fna.gz`, which that
  directory does not hold. It writes an empty table and reports "0 assemblies".
- **1c:** the same glob problem. `--report` is a flag, but RUNBOOK gives it a
  path, so argparse fails. The docstring also says `--report` "prints the
  distribution"; the code does not.
- **1g:** `mash_nearest.sh data/rotated` stops with "no *.fna.gz".
  `snp_counts.tsv` has no producer (PGB-8).
- **1e runs before the panel exists,** using the previous panel FASTA as its
  background.
- **Missing steps:** rotation (`rotate_to_dnaa.sh`) and PanSN assembly
  (`make_pansn_fasta.sh`) are not in RUNBOOK at all.
- **PIPELINE_SEGMENTS.md** says the rotation check is "in collinearity QC". It is
  written by `rotate_to_dnaa.sh`.

**Suggested fix:** rewrite Segment 1 from a dry run of the actual chain. Better,
make Segment 1 one driver script with an audit of expected products (lesson
10).

## PGB-21 (LOW): `ppe38_build_panel.py` cuts each H37Rv ORF one base too long

**Location:** `bin/ppe38_build_panel.py:21-29` and `:66`.

**What is wrong.** The gene starts are copied from the snpEff dump, which is
0-based (`H37Rv_snpeff_dump.txt`: Rv2354 `2635627 2635954`, CDS
`2635628..2635951`). They are then treated as 1-based (`fetch(s0 - 1, e0)`).

**Evidence** (`data/annotation/ppe38/orf_panel.fasta.fai`):

| ORF | panel length | true length | note |
|---|---:|---:|---|
| Rv2354 | 328 | 327 | starts `CAT...`, the ATG plus one upstream base |
| glyS | 1,393 | 1,392 | |
| PPE40 | 1,849 | 1,848 | |

All six H37Rv-cut ORFs have the extra base.

**Effect:** one flanking base per ORF in the locus profiling. Negligible.

**Suggested fix:** use `fetch(start0, end)` with the dump's 0-based start, or
take the coordinates from a GFF.

## PGB-22 (LOW): conflicted genotypes are silently set to missing

**Location:** `bin/vg_deconstruct.sh:62`; there is no `-K` option.

**Evidence:**

- `variants.vcf` has 3,713 records with `CONFLICT`, naming 134,328 sample cells.
  Without `-K`, vg writes those genotypes as missing. These are genomes with
  more than one traversal (duplications, repeat copies).
- The VCF header records no vg version or command line.

**Suggested fix:**

- record the vg command and version in the header (`stamp_build_id.sh` or
  `bcftools annotate`);
- count CONFLICT cells per sample in the build report;
- any reader that treats missing as REF must treat these as unknown.

## PGB-23 (LOW): screen outputs that are missing for a genome mean "pass"

**Location:**

- `bin/build_panel.py:109-114`: no collinearity row means not scrambled;
- `bin/collinearity_qc.py:100-102`: a minimap2 failure drops the genome with a
  stderr line;
- `bin/assembly_provenance_screen.py:75-78` and `:155-157`: a failed NCBI query
  means technology `unknown`, which passes; a missing assembly is skipped;
- `bin/foreign_insertion_screen.py:82-83`: a missing assembly is skipped;
- `bin/qc_master_table.py:46-50`: a missing table gives `{}`.

**Suggested fix:** `build_panel.py` should require a row from every screen for
every candidate, and refuse to write a panel otherwise.

## PGB-24 (LOW): provenance gaps

- **CX333 has no `graph_provenance.tsv`.**
- **The recorded panel commit cannot be the code that built the panel.**
  `refbias/build/7713a8d71d8e/build_info.tsv` records `panel_commit 1cd160e`,
  inferred as "last commit before graph created 2026-09-17". But:
  - the graph's files date from 2026-09-09 to 2026-09-11 (`params.yml`
    2026-09-09 18:00);
  - the working tree's first commit (`741e894`) is 2026-09-10 16:45, after the
    build began.
- **The container is recorded by path only:** `pggb_latest.sif`, not a version
  or digest. `params.yml` does say v0.7.4.
- **The 2025 selection filter (468 to 416) is unrecorded,** as GRAPH_PROVENANCE
  section 8 already says.

**Suggested fix:**

- at the rebuild, write the sidecar after the build (PGB-1), with the input
  FASTA sha256, the container sha256 and the list of panel accessions;
- correct `build_info.tsv`'s `panel_commit_basis` to "none: built before the
  first commit".

---

## What must change before a new panel build

Grouped by the HANDOFF lessons in sections 0f, 0g, 0j and 0k, with the finding
IDs above.

**Panel selection (0g lessons 3 to 5, 12; PGB-3, -4, -13, -20, -23):**

1. Select candidates by MTBC taxid from RefSeq and GenBank, so M. orygis and
   GenBank-only genomes are included.
2. Make one `build_panel.py` read every screen table and write the PanSN FASTA.
   It should require a row per screen, record panel provenance (accession list,
   sha256s, commit) at build time, and leave no hand step between 484 and 333.
3. Fix RUNBOOK Segment 1, or replace it with a chain script that audits its
   products.

**Assembly QC (0f items 1 to 3; 0g lesson 1; PGB-6, -7, -8, -16, -19):**

4. **Run-length and low-complexity rule:** windowed (at least 90% one base over
   200 bp, or GC below 35% over 500 bp). On CX333 it flags exactly GCF_050259585,
   GCF_045348265 and GCF_039770655, and nothing else. Mask with N or reject.
5. **Long-read error checks in `bin/` and the chain:** indel excess, private
   homopolymer-indel excess, and agreement with the genome's own short-read
   assembly where reads exist. All of these are within sublineage, with a
   robust z, a MAD-0 fallback, and counts from direct alignment, not a graph.
6. **Technology from sequence:** insertion count of 50 bp or more, novel IS6110
   copies, and the homopolymer signature. Metadata only supports the call; fix
   the keyword lists and flag short-read assemblers.
7. **Mixed-lineage test within a root;** mask IUPAC codes.

**Foreign-DNA screen (0f item 4; PGB-5):**

8. A small, PanSN-named, one-per-sublineage background.
9. Best single-genome homologue, not the sum.
10. Re-check with the frequency filter off and against IS6110
    (`foreign_recheck.py`), plus a vector database.
11. Extract gaps between alignment blocks.
12. Exempt the outgroup in code.

**The two known artifact genomes and the third (0f; PGB-6):**

13. Mask or replace GCF_045348265, GCF_050259585 and GCF_039770655. Remember
    that GCF_039770655 is the only lineage 9 panel genome and the reference for
    every lineage 9 isolate.

**Graph and variants (0g lesson 6; PGB-1, -2, -9, -10, -11, -17, -22):**

14. Fix `pggb_build.sh`'s guard order and defaults (`-k 23 -K 15`, 48 cores).
    Record the container digest.
15. Publish one variant file: collapsed and left-aligned, with LV rewritten after
    vcfbub and the decomposition settings in the header. Move every graph-VCF
    reader to it.
16. Put vcfwave's large-allele fallback in `vcf_decompose.sh`. Write each output
    under a temporary name and rename it when the step finishes.

**Outgroup and rooting (0j, 0k):**

17. Keep the canettii and lineage 8 genomes out of the graph and in the kit
    (0k); `add_outgroup.py` then reads the kit, not the graph VCF.
18. Define the ancestral node as the ingroup MRCA (fault B; tree audit).

**Downstream frame (PGB-15):**

19. Stage P0 references from the panel FASTA, so the graph and alignment frames
    are identical.

**Provenance (PGB-14, -24):**

20. `stamp_build_id.sh` must check the graph hash.

---

## Checked and found sound

- **Rotation.** All 333 panel genomes are dnaA at position ≤ 100 on the forward
  strand (`rotation_verify.tsv`); verification maps dnaA back and records strand
  and offset.
- **`build_panel.py` with full inputs** (`panel_manifest.tsv`, `--barcode`,
  `--manual-exclusions`) reproduces `panel.rebuild.tsv` (484) exactly. It keeps
  8 of 8 of the inversion clade, 5 H37Rv and 3 Erdman. The BCG barcode rule
  catches all 10 barcode-BCG genomes.
- **Production `tech_class`** is reproduced by the code's `tech_class()` for all
  484 genomes.
- **`vcfbub`** leaves no parent and child together. All 49,120 promoted LV1
  records lack their parent, so no event is represented twice before vcfwave.
- **The collapse step** is nearly lossless: 93,214 keys, carriers equal to the
  decomposed union except 64 cells at 12 keys.
- **The class-split predicate** (`vcf_split_classes.sh`) matches QC_PIPELINE.md.
  The 50 bp and block-substitution clauses are correct as written.
- **`make_pansn_fasta.sh`** writes PanSN names and checks for duplicate names and
  the reference path. It derives `-n` from the FASTA index.
- **`graph_panel_check.sh`** compares path sets, not counts. It does not compare
  sequence content; a sha256 per path would close that gap.
- **`snp_outlier_screen.py`** uses a robust z within sublineage, with a fallback
  to the major lineage. It reports the excess and deficit tails.
- **`mash_nearest.sh`** records its parameters, takes the nearest non-self
  neighbour, and keys on accession.
- **`collinearity_qc.py`** handles a single full-length block correctly. The
  fixed short-circuit is confirmed: 445 single-block genomes have coverage
  1.0000.
- **The `vcf_decompose.sh` merge** takes its header from a waved chunk, and stale
  `part_*` files are removed before splitting.

## Open questions

1. **What the provenance screen code would exclude on the 484.** Not run here,
   about 8 CPU-hours. It would show how far the coded rule (insertions of 50 bp
   or more, fewer than 5) is from the production rule (identical IS6110 profile
   plus short-read technology).
2. **What the 3.39 Mb top-level snarl at H37Rv 415,369 to 3,803,919 is.**
   UNVERIFIED; most likely the B0/W148 rearrangement.
3. **Whether the `(CCATT)n`/poly-T in GCF_039770655 is human satellite
   contamination or an ONT artifact.** The accessory-presence ALT calls in 9
   gwas1000 and 5 scale200 isolates suggest some reads match it.
4. **Whether downstream readers outside this area rely on `LV`.** For example,
   whether anything filters on LV in the collapsed or nolab files, where LV is
   taken from an arbitrary record. That belongs to the graph-VCF reader audit.
