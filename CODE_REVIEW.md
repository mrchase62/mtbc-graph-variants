# Code review: mtbc-graph-variants

Review date: 2026-09-29. Reviewed at commit `6ab0f8d`.

The tracked code is byte-identical to the working tree at
`/n/boslfs02/LABS/sfortune_lab/Lab/mchase/MtbPangenome` apart from `README.md`,
`.gitignore`, and one line in `bin/sync_back.sh`. So every finding here applies
to the code currently producing results.

## How this review was done

- `bash -n`, `py_compile`, `pyflakes` and `shellcheck` over all 65 files.
- A transitive dependency scan: every script, module and import the code uses,
  checked against what the repository contains.
- A line-by-line read of each pipeline area for correctness defects.
- Where cheap, a finding was confirmed against real outputs in the working tree
  (`refbias/{pilot,scale100,scale200,gwas1000}`) or with a small scratch test.
  These are marked **[measured]**. The others come from reading the code paths
  and are marked **[read]**.

Severity levels:

- **HIGH**: the pipeline cannot run, or it writes wrong genotypes or coordinates.
- **MEDIUM**: wrong results in a definable subset of cases, or failures that go unnoticed.
- **LOW**: edge cases, diagnostics, robustness.

## Summary

The most important conclusions:

1. **A fresh clone cannot run end to end.** The extraction computed the
   dependency closure from shell calls and missed Python imports (section 1).
2. **Merged VCFs already produced are not spec-valid, and some cells are wrong.**
   Findings 3.1, 3.2, 3.3, 4.1 and 4.3 affect gwas1000 and scale200 today.
   Regenerate those cohorts after the fixes, before any association work.
3. **Most defects come from a few repeated patterns:**
   - "not measured" written as ABSENT or REF;
   - strand and base-numbering (0- vs 1-based) handled separately in each script;
   - skip-if-exists "done" checks that accept partial or stale outputs;
   - subprocess exit codes that are never checked.

   Section 8 proposes shared fixes for each.

---

## 1. Repository completeness (HIGH)

### 1.1 Imported modules missing from the repository [measured]

These files are required at runtime and exist only in the working tree:

| missing file | needed by |
|---|---|
| `bin/mtb_norm.py` | `p5_keys.py`, `p5_states.py`, `is6110_p5_merge.py` |
| `bin/p5_states_io.py` | `p5_matrix.py`, `merge_cohort_vcf.py` |
| `graphframe/bin/graph_frame.py` | `frame_convert.py`, `frame_detect.py`, `is6110_project_sites.py`, `is6110_p5_stage2.py` |
| `graphframe/bin/graph_frame_offsets.py` | `graph_frame.py` |
| `is6110/bin/p1i_discover_matched.sh` | named in `refbias_run.sh` as the step that rebuilds matched references |

P5 and the IS6110 arm fail with `ModuleNotFoundError` on a fresh clone.

**Fix:** copy these in, then re-run a closure scan that follows Python
`import`, `from ... import` and `importlib` loads as well as shell calls.

### 1.2 Documented helpers that are not here [measured]

The README and INPUTS.md tell the user to run these, but they are absent:

- `bin/show_config.sh` (setup step 2)
- `bin/stage_in.sh`
- `bin/pggb_build.sh`
- `bin/build_repeat_mask.py`, which INPUTS.md says must generate the repeat mask
- the QC stage scripts listed in INPUTS.md, and `QC_PIPELINE.md`

**Fix:** either include them, or say explicitly that they live in the
companion graph repository and link to it.

### 1.3 Other scripts referenced only in comments [read]

These do not break a run, but readers are sent to files they cannot find:

- `is6110_family_depth.py` and `panisa_score.py` (named in `is6110_junctions.py`)
- `assemble_from_cram.sh`, `write_event_matrix.py`, `vcf_decompose.sh`
- `refbias_t5.sh`, `run_block.sh`, `t16_graph_frame.py`

---

## 2. Configuration (`config/project_env.sh`)

### 2.1 The site file is sourced after the derived paths are set (HIGH) [measured]

`MTB_DATA`, `MTB_GRAPHS`, `MTB_LOGS` and the other working subdirectories are
assigned from `MTB_WORK` *before* `site.local.sh` is sourced. A site file that
sets `MTB_WORK=/elsewhere` therefore gives the following, which is reproduced:

    WORK=/elsewhere  DATA=/n/netscratch/.../MtbPangenome/data  GRAPHS=/n/netscratch/.../graphs

This contradicts the file's own comment that the site file "is sourced before
the defaults below", and it defeats the goal of making the published config
portable.

**Fix:** source the site file at the very top, before the Roots section.

### 2.2 The site-file path is relative to the current directory (HIGH) [measured]

`${MTB_SITE_FILE:-config/site.local.sh}` resolves against `$PWD`, not against
the repository. Sourcing the config from any other directory silently skips the
site file, leaving `MTB_CRAM_ROOT` empty.

**Fix:** resolve it from `$(dirname "${BASH_SOURCE[0]}")/site.local.sh`.

### 2.3 Tool variables used but not declared (LOW) [measured]

Scripts read the following, and fall back to hard-coded defaults when unset:

- `MTB_BWA` and `MTB_GATK_SIF`, which INPUTS.md lists as configuration
- `MTB_WGSIM`, `MTB_TMPBASE`, `MTB_BUILD_DIR` and `MTB_SCRATCH_ROOT`

`show_config` cannot report any of them.

**Fix:** declare them all in `project_env.sh`.

### 2.4 Hard-coded cluster paths outside the config (LOW) [measured]

Despite the rule in the header comment, these scripts carry a literal
`/n/boslfs02/...` default:

- `p0_prepare.sh:55` and `simulate_and_call.sh:36` (GATK container)
- `p2_call.sh:68` (delly env)
- `p3_accessory.sh:79`, `p4_place.sh:66`, `p4b_place_sv.sh:40`,
  `p5_merge.sh:70` and `p5_svgt.sh:79` (odgi)
- `accessory/bin/locus_presence.py:72` (bwa)

`project_env.sh` itself also defaults to lab and home-directory paths for
`MTB_PY`, `MTB_PY_VT`, `MTB_ODGI`, `MTB_K8` and others. For a published
repository, move all of these into `site.local.sh.example`.

---

## 3. P5 merge, genotyping and the VCF writer

### 3.1 GT=2 on records with only one ALT allele (HIGH) [measured]

In `bin/merge_cohort_vcf.py`, only the small-variant block adds a `*` allele
when a sample is ABSENT (line ~332). The interval, caller-SV and IS6110 blocks
(lines ~424, ~451, ~540) map ABSENT to GT=2 but keep a single ALT such as
`<DEL>` or `<INS:ME:IS6110>`.

- In scale200, 976 records (487 SV and 489 IS6110) have a GT index larger than
  their ALT count.
- `bcftools +fill-tags` aborts with `Incorrect allele ("2")`.

**Fix:** add `*` in every block whenever any cell is ABSENT, or encode ABSENT
as `.` and carry the state in FORMAT/ST.

### 3.2 Node-frame REF reads depth at the wrong base (HIGH) [read]

In `bin/p5_states.py:285-288`, the minus-strand position is `start - off`,
which lands outside the node whenever `off > 0`. It should be
`start + len - 1 - off`. 227,401 of 438,234 node positions are on the minus
strand.

`start` also comes from the GFA path walk, which is the panel frame, while the
gVCF is in the refs frame. The H37Rv branch converts between the two with
`frame_convert`; the node branch does not. 110 of 333 accessions are rotated
between the frames.

**Fix:** convert panel to refs, and use the correct strand formula.

**Correction after testing against the data.** The minus-strand formula
proposed above, `start + len - 1 - off`, is also wrong. P4's node offsets
already run along the path, so the base is `start + off` on both strands. The
test read the refs base at the computed position and compared it with the
key's reference allele, over the pilot's references:

| node walked | old code | proposed above | adopted |
|---|---:|---:|---:|
| + | 56.9% | 100.0% (with frame conversion) | 100.0% |
| - | 15.3% | 5.0% | 99.8% |

### 3.3 Any covered position is stated REF (MEDIUM-HIGH) [measured, in part]

In `bin/p5_states.py:234-294`, with `load_gvcf_blocks` at line 49, coverage
counts every gVCF line, including variant lines with GT=1. Only the exact key
is looked up.

- At a multi-allelic site, a sample carrying allele A is written as REF for
  allele G. scale200 has 207 such contradictory cell pairs. [measured]
- A call that P4 dropped also becomes REF: FILTER not PASS, non-ACGT allele, or
  dropped by routing. That contradicts the VCF header's meaning of GT=. [read]

**Fix:**

- Count only GT=0 gVCF lines as reference coverage.
- Index placed records by position. Where the sample has a different allele,
  write NOCALL, or GT=0 with ST=OTHER.

### 3.4 Pilot IS6110 table used as a fallback for any cohort (HIGH) [read]

At `bin/p5_finish.sh:185-190`, a cohort without its own
`is6110/results/<tag>_p1i_cohort_keys.tsv` falls back to the pilot's
`p1i_cohort_keys.tsv`. The value is then passed explicitly, which bypasses the
guard in `merge_cohort_vcf.py:217-225` that was written to stop exactly this.

`merge_cohort_vcf.py:242-244` also adds every sample in the key table to the
sample order. A cohort can therefore gain the pilot's 23 isolates as extra
columns.

**Fix:**

- Remove the fallback and let the merge resolve the table.
- Restrict sample order to the refmap.

### 3.5 Symbolic records: REF is N, POS is one base late, no END (MEDIUM) [measured]

In `merge_cohort_vcf.py` (lines ~426, 452, 498, 541), 7,291 scale200 records
have REF=`N`. That does not match NC_000962.3.

The interval start is the first deleted base (`sv_intervals.py:100, 159`), but
VCF puts the padding base at POS, so every deletion is shifted by +1. Without
END, region queries inside a deletion miss the record.

**Fix:** set POS=start-1, REF to the H37Rv base there, and add INFO/END.

### 3.6 Stale per-sample outputs are skipped rather than recomputed (MEDIUM) [read]

- `p5_merge.sh:130` skips a sample whose states file exists, without comparing
  its `#keys_sha1` to `keys.tsv`. After a key change every task reports "already
  done", and `p5_finish` then fails until files are deleted by hand.
- `p5_svgt.sh:133` has the same skip. Nothing ties `svgt_states.tsv` to
  `sv_matrix.tsv`. SV keys (`sv:TYPE:median_pos:median_len`,
  `p5_sv_matrix.py:316`) shift when cluster membership changes, so a stale
  overlay either fails to join, silently, or joins the wrong cluster.

**Fix:** check the hash or matrix identity before skipping.

### 3.7 Tied ancestral alleles are mostly resolvable (MEDIUM) [measured]

In `bin/ancestral_alleles.py:143-149`, the root is the canettii/MTBC split. The
root set ties whenever the single outgroup differs from the MTBC ancestor.
4,759 of 5,130 TIED sites have an unambiguous state at the MTBC ingroup node,
yet they are written `AA=.`.

The unused locals flagged by pyflakes (lines 46, 51, 70) are harmless. The
parser was tested and is correct.

**Fix:** report the ingroup node's state after a downward Fitch pass from the
outgroup.

### 3.8 Nested genes missed by annotation (LOW-MEDIUM) [measured]

In `bin/p6_annotate.py:64-70`, the backward scan stops at the first gene that
ends before the position. So a longer, earlier gene that contains the position
is never reached. 9,418 genic bp of H37Rv come back as intergenic, and
`pseudogene` features are excluded by `want`.

The unused `fixed` at line 151 is not a dropped fix.

**Fix:** use an interval tree, or scan back by the maximum gene length.

### 3.9 Caller-matrix mode promotes inherited ALT without confirmation (LOW-MEDIUM) [read]

Interval mode in `bin/p5_sv_genotype.py` requires two-frame evidence, after a
measured 26.5% contradiction rate. The caller-matrix fallback (lines 441-445)
still promotes on covered flanks alone.

**Fix:** apply the same two-frame requirement in both modes.

### 3.10 Minor issues

- `p5_sanity.py:158`: `abs(x-med) > 3*med` can never flag a sample below the
  median, such as one with zero REF because its gVCF is broken.
- `merge_cohort_vcf.py:631`: tabix runs with `check=False`, so a failed index
  goes unnoticed.
- `merge_cohort_vcf.py:567`: the CLASS header text says accessory records are
  absent, but `CLASS=accessory_presence` records are written.
- Node-frame records sit at POS=1 on contigs declared with length 1 but carry
  multi-base REF. The comments acknowledge this, but the VCF is still invalid.
  Declare node contigs at their true GFA length and use `node_offset+1`.

---

## 4. IS6110 and accessory arms

### 4.1 Every node-frame non-carrier written as ABSENT (HIGH) [measured]

In `is6110/bin/is6110_p5_stage2.py`, `short()` (lines 98-100) reduces the
stage-1 key to `node:<id>:`. The carrier table (line 93) builds
`node:{x["node"]}:`, where `x["node"]` is already `id:offset` (see
`is6110_write_vcf.py:174`). The two forms never match:

1. No carrier is found.
2. Nothing is queried.
3. The missing projection falls back to `(None, -1)`.
4. That is written as ABSENT, "no projection".

Every node-frame NOCALL became ABSENT, and none became REF or NOCALL:

| cohort | false ABSENT cells |
|---|---:|
| gwas1000 | 446,005 |
| scale200 | 34,362 |
| scale100 | 13,821 |
| pilot | 1,061 |

**Fix:** make `short()` return `node:{f[1]}:{f[2]}:`, and see 4.2.

### 4.2 "Not measured" written as ABSENT (HIGH) [read]

In `is6110_p5_stage2.py:117-120, 133-154, 229-232`, the `odgi position` return
code is never checked. If odgi fails, every key for that target reference
becomes ABSENT. The same happens when the target or carrier reference is
missing from `paths`.

**Fix:** exit on a nonzero return code. Reserve ABSENT for `dist != 0` from a
query that succeeded, and write NOCALL otherwise.

### 4.3 Site projection run without `--all-stacks` (HIGH) [measured]

`bin/p1i_vcf.sh:78-80` calls `is6110_project_sites.py` without `--all-stacks`,
although the output is named `_sites_h37rv_all.tsv`. So only A/B-tier sites are
projected. The others go through two more steps:

1. `is6110_write_vcf.py:169-198` writes them as `unplaced` with key `node:`.
2. `is6110_p5_merge.py:74-76` then drops them silently.

- gwas1000: 1,800 of 11,524 key rows have an empty node key, 816 of them ALT.
- scale200: 319 such rows.
- pilot and scale100: none, because they were run by hand with the flag.

**Fix:** pass `--all-stacks`, and have the writer refuse to emit an empty node
key.

### 4.4 Reference-shared vs reference-lacking decided inconsistently (MEDIUM) [measured]

`is6110_write_vcf.py:168` and `is6110_place_by_flank.py:144-147` require an
exact match on the seam base. `is6110_promote_sites.py` allows 25 bp of slop.
295 gwas1000 sites sit 1-25 bp from a seam and are classified differently by
the two. They are written as ALT with an ALT cohort key, and flank placement
does not step past the reference's own copy of the element.

**Fix:** use one shared nearest-seam helper with the same slop everywhere.

### 4.5 One insertion split across several nearby keys (MEDIUM) [measured]

`is6110_write_vcf.py:197-198` keys on the stack peak, which can sit on either
side of the target-site duplication. In gwas1000, 536 of 1,800 H37Rv keys have
another key 1-10 bp away, 185 of them exactly 1 bp away. This inflates the
convergence counts.

**Fix:** report a canonical position, such as the left seam, or cluster keys
within about 10 bp across samples before the P5 merge.

### 4.6 Samples with no candidate sites crash, and the done marker is written too early (MEDIUM) [read]

`is6110_junctions.py:300-302` returns without writing `--out` when there are no
rows. The following `mv -f "${JUNC}.tmp"` then fails under `set -e`, at
`p1i_matched.sh:269-273` and `p1g_isclean.sh:90-94`. The sample loses its
element-side scan and gets an empty VCF: "no data" is conflated with "nothing
there".

The junction file is also the "already done" marker, but it is written before
`elstacks.tsv` and `elementdepth.tsv`.

**Fix:** always write a header. Write the done marker last.

### 4.7 Accessory presence never checks subprocess exit codes (MEDIUM) [read]

`accessory/bin/locus_presence.py:125-172` runs `samtools view` with stderr
discarded, then `bwa`, `sort` and `depth`, and checks none of their exit codes.
A failed CRAM read gives zero coverage, and all 802 loci are called ABSENT.
There is also no depth floor, so a low-depth sample's loci are ABSENT rather
than UNCERTAIN.

**Fix:**

- Check every return code.
- Mark the sample NOCALL when the query pool is empty or implausibly small.

### 4.8 Copy-number window projected from the start only (MEDIUM) [read]

`bin/p3_copynumber.py:94-95` projects only the locus start and always reads
forward, `rpos .. rpos+L`, which is also L+1 bases. Where the reference runs
reverse to H37Rv, it measures the wrong flank.

**Fix:** project both ends and sort them, as `p3_element_depth.py` does.

### 4.9 Element-depth span not checked (LOW-MEDIUM) [read]

`p3_element_depth.py:202-217` does not check that the projected span is close
to the element's length. The termini sit on collapsed IS6110 nodes, so the two
ends can land in different copies. The window can then be hundreds of kb, and
the depth ratio comes out near 1.

**Fix:** report `no_projection` when the span length differs from the element
length.

### 4.10 Minor issues

- `is6110_write_vcf.py:137`: `hcontig` is read and then discarded. This is
  harmless today, but overriding `--h37rv-contig` would write a mismatched CHROM
  without any error.
- `is6110_write_vcf.py:176` hard-codes MEINFO polarity `+`, although orientation
  is available from the element-side pair counts.
- `is6110_reconcile.py:214` silently uses an identity crossmap when a sample's
  crossmap is missing.

---

## 5. P4/P4b projection and SV placement

### 5.1 Repeat-mask lookup misses nested intervals and is off by one (HIGH) [measured]

`make_hit` in `bin/p4_place.py:80-85`, used at lines 218-223, tests only the
last interval starting at or before the position. The mask has heavily nested
entries: PE_PGRS4 336359-339273 contains 336559-339142, so position 339200 is
missed.

It also compares 1-based VCF positions to 0-based half-open BED with
`s <= p < e`. The correct test is `s < p <= e`.

- 56,926 of 438,885 masked bp (13%) are routed as core.
- That includes 53,429 of 334,640 PE/PPE bp.

In those regions the direct arm is kept and the composed arm is dropped, which
inverts the routing rule.

**Fix:** merge intervals per class before bisecting, and use `s < p <= e`.

A related issue at line 75: the name test `"PE" in name` matches any name that
merely contains "PE". Match on a prefix or on an explicit class column instead.

### 5.2 SV placement ignores strand (HIGH) [read + toy test]

`bin/p4b_place_sv.py:138-169` reads only columns 0-2 of the odgi output. It
never uses the strand. A deletion at POS=1000, END=3000 on a flipped reference
is emitted as `sv:DEL:4002:2002`, with start greater than end. The consequences:

- `p5_sv_matrix.py:260` clusters on `h37rv_pos`, so the event never joins
  forward-strand carriers.
- `sv_intervals.py:159` builds the depth interval entirely outside the real
  deletion.

**Fix:** for the minus strand, take `pos = min-1` and `end = max-1`, and handle
INS the same way. Require `|end-pos|` to be close to SVLEN before a record
counts as "placed".

### 5.3 A VCF with no testable records kills the sample (HIGH) [measured]

`frame_detect.py` exits 2 when it finds no single-base REF records. Under
`set -euo pipefail`, `VFRAME="$(... | cut -f2)"` aborts the script, and the
FATAL branch never prints. This is at `p4_place.sh:116-117` and
`p4b_place_sv.sh:75-78`.

- **P4:** this runs before the legitimate empty-VCF branch at lines 140-146, so
  that branch is unreachable for the 110 of 333 accessions whose refs and
  panel frames differ.
- **P4b:** a header-only delly or dysgu VCF passes `[[ -s ]]` and dies the same
  way.

**Fix:** skip frame detection when there are no data records, and capture the
exit code explicitly.

### 5.4 Samples with no SV calls also lose the inherited SV half (MEDIUM) [read]

`p4b_place_sv.py:121-132` returns early with a header-only table before the
inherited half runs. The reference's own large differences from H37Rv are
dropped for exactly the samples closest to their reference.

**Fix:** remove the early return.

### 5.5 A missing graph VCF silently empties the inherited halves (MEDIUM) [read]

`--graph-vcf` defaults to a relative `graphs/CX333...` path, and neither shell
wrapper passes it. `os.path.exists` then skips the block with no warning. This
is at `p4_place.py:144-145, 327` and `p4b_place_sv.py:113-114, 211`.

**Fix:** pass the path from `$MTB_GRAPHS`, and make a missing file fatal.

### 5.6 sv_twoframe never checks samtools exit status (MEDIUM) [read]

In `sv2frame/bin/sv_twoframe.py:84-97, 114-133`, stderr goes to DEVNULL and the
return code is ignored.

- A CRAM decode error partway through leaves zero depth over the rest of the
  contig, which yields false `depth_absent` ALT calls.
- A failed `view` yields zero clips, so no call is ever marked strong.

**Fix:** check the return code and fail on nonzero.

### 5.7 Minor issues

- `p4b_summary.py:65-67`: ZeroDivisionError when every table is header-only.
  Its "placed" count also includes inherited rows, which were never projected,
  so it overstates the projection rate.
- `p4b_place_sv.sh:32` picks the build with an unsorted `find | head -1`, while
  `p4_place.sh:53-54` refuses when there are several builds.
- `p4_place.py:202`: without `--h37rv`, every minus-strand indel gets the status
  `off_the_end`, not "left as called" as documented.
- `p4_place.py:276`: records with an H37Rv projection but no node projection
  are dropped without being counted.
- `p4b_place_sv.py:174-175`: `k1` and `k2` are dead code. The mask loaded at
  lines 152-157 is never used either, so SVs get no region label.
- `sv_intervals.py:127-129`: after merging, set `svlen = end-start+1`. The
  `svi:` ID currently carries the maximum SVLEN, which can disagree with the
  merged span.

---

## 6. Orchestration and early stages (P0-P3)

### 6.1 A failed sbatch drops downstream dependencies (HIGH) [measured]

`submit()` in `bin/refbias_run.sh:277-279` does `id="$(sbatch ...)"; echo "$id"`.
`set -e` does not apply inside the `$(...)` without `inherit_errexit`, so a
rejected submission returns an empty ID with exit 0. `deps` then drops the
empty ID, and the next job is submitted with no `--dependency`.

- **Scenario:** a QOS limit rejects `p2a`, and `p2sum` runs at once on old
  outputs. P3 and P4 then chain on stale P2 files.

**Fix:** `id=$(sbatch ...) || die`, and require `[[ $id =~ ^[0-9]+ ]]`.

### 6.2 The pre-VCF P5 step does not wait for P4b (HIGH) [read]

At `refbias_run.sh:343-347`, `p5pre` builds `sv_matrix.tsv` from P4b output,
but depends only on `p5states`, and so only on P4. P4 and P4b run in parallel
off P2. The P4b guard in `p5_finish.sh:85-93` only compares file ages with the
script.

**Fix:** submit `p5pre` with a dependency on both `p5states` and `p4b`.

### 6.3 The documented command line is rejected (MEDIUM) [measured]

README.md and INPUTS.md say `bash bin/refbias_run.sh --cohort <name>`.
`refbias_run.sh:98` rejects `--cohort` as an unknown option; the cohort is
positional.

**Fix:** accept `--cohort`, or correct the documentation.

### 6.4 P1 and P2 ignore the documented CRAM reference variable (MEDIUM) [measured]

`p1_select_reference.sh:64` and `p2_call.sh:65` default `CRAMREF` to
`${CRAMROOT}/metadata/reference.fasta`, which is the old collection's layout,
and ignore `MTB_CRAM_REF`. P1g and P1i already use the variable. Neither script
calls `mtb_require_cram_root`.

**Fix:** `CRAMREF="${CRAMREF:-$MTB_CRAM_REF}"`, plus a call to
`mtb_require_cram_root`.

### 6.5 Non-atomic outputs accepted as complete on rerun (MEDIUM) [read]

| where | problem |
|---|---|
| `p2_call.sh:119, 152, 158-160` | "Done" means `vcf.gz` and `delly.vcf` exist. A failed dysgu run keeps its partial VCF. A job killed during dysgu is skipped on rerun. `bcftools view > delly.vcf` is not atomic. |
| `p1_select_reference.sh:109-113`, `p2_call.sh:125-129` | FASTQs are extracted straight into their final names and reused whenever they are non-empty. A timeout leaves truncated reads, which are then aligned on rerun. |
| `p0_prepare.sh:212, 229-245` | `zcat > fa` is not atomic. `refs.done` is written even after `dict_failed` is logged. Leftover `*.bwaindex.PID.bwt` files inflate the count. |

**Fix:** write to a `.tmp` name, `mv` on success, and write a `.done` marker
last.

### 6.6 Real reads can be silently replaced by simulated ones (MEDIUM) [read]

At `simulate_and_call.sh:53`, if FQ2 is missing, the script simulates reads
with wgsim from `$SRC`, which in P1 and P2 is the reference. The isolate's real
reads would be replaced by reads from the reference itself.

**Fix:** refuse to simulate when `SIM_FQ_PREFIX` is set.

### 6.7 Blank lines inflate the array size (MEDIUM) [read]

`refbias_run.sh:114` counts rows with `grep -vc '^#'`, which includes blank
lines. A trailing blank line gives N+1 tasks. The last task dies with "no
cohort row", and afterok cancels the whole chain.

**Fix:** count and index with the same filter, for example
`awk 'NR>1 && $1!=""'`.

### 6.8 The build is chosen by modification time (MEDIUM) [read]

`refbias_run.sh:118` picks the newest build directory (`ls -dt`). Re-running
one `p0_prepare.sh` step on an old build makes that build newest, and the chain
then runs against the old graph.

**Fix:** require `MTB_BUILD_DIR`, or exactly one build with a `manifest.done`.

### 6.9 Concurrent panel indexing in P3 (MEDIUM) [read]

At `p3_accessory.sh:45-55, 100-111`, array tasks race to index the panel, and
`.bwt` is moved into place before `.pac` and `.sa`. Another task can start
`bwa mem` on a half-built index, and one failed task cancels everything
downstream.

**Fix:** build the index in P0, or move `.bwt` last and use `flock`.

### 6.10 Summaries exit 0 when samples are missing (MEDIUM) [read]

`p1_summary.py:130-134, 187` and `p2_summary.py:71-72, 83` list missing samples
but exit 0. Samples then drop out of P5 silently, because P5 derives its
expected count from the refmap.

**Fix:** exit nonzero on any missing sample unless an explicit flag allows it.

### 6.11 Minor issues

- `refbias_run.sh:208-211`: `--status` exits under pipefail at the first empty
  pass directory. `pat` has no default case, so passes p1g..p5vcf reuse p5's
  pattern.
- `sync_back.sh`: there is no `--delete`, so there is no risk of wiping durable
  storage. But:
  - a regressed file on scratch overwrites the durable copy with no backup;
  - a failed `git push` (lines 305-307) still exits 0 and appends to
    `.last_sync`;
  - `|| true` at line 143 hides rsync failures.

  Use `--backup --backup-dir=<date>`, and exit nonzero when the push fails.
- `t8_select_reference.py:64-71`: an uncovered site counts as REF, and
  `argsort` breaks distance ties arbitrarily. Use a covered-sites mask and
  `kind="stable"`.
- `p2_call.sh:169`: quote `$OUTDIR` (SC2046). The SC2140 warnings at lines
  172-173 are false positives.

### 6.12 Checked and found correct

- HaplotypeCaller runs with `-ploidy 1` for both the VCF and the GVCF.
- The mapping from SLURM task ID to table row skips the header correctly, apart
  from the blank-line case in 6.7.
- Per-sample temp names include the sample and PID, so array tasks do not
  collide.
- Per-cohort P4/P4b/P5 directories are passed correctly from the runner.

---

## 7. Lint

pyflakes reports about 80 cosmetic warnings: unused imports and f-strings
without placeholders. Only the unused variables were investigated. They are
covered above in 3.7, 3.8, 4.10 and 5.7, and none hides a dropped computation.

Suggested tooling:

- Add `ruff` and `shellcheck` to a pre-commit hook.
- Pin the Python version, since the login node's default `python3` is 3.6.

---

## 8. Recommendations for the pipeline as a whole

1. **Gate every merged VCF with a validator.** Run the following as a
   `p5_finish.sh` step, and fail the pass on any error:

       bcftools +fill-tags merged.vcf.gz -Ou -- -t AN,AC > /dev/null
       bcftools norm -c e -f H37Rv.fasta merged.vcf.gz -Ou > /dev/null

   Also add a check that every GT index is at most the number of ALT alleles.
   Findings 3.1 and 3.5 would have failed this.

2. **Separate "not measured" from "absent" everywhere.** Findings 3.3, 4.1,
   4.2, 4.6, 4.7 and 5.6 all turn a failed or missing measurement into ABSENT or
   REF.
   - Adopt one rule: ABSENT and REF require a successful measurement, and
     everything else is NOCALL.
   - Make every subprocess call check its return code. A small `run()` helper
     in Python and `set -o pipefail` plus explicit checks in shell cover most
     cases.

3. **Share geometry helpers instead of re-deriving them.**
   - One interval module with the base system stated explicitly and merged
     lookups, tested on nested PE_PGRS intervals.
   - One "sample state at a projected position" helper covering GT-aware
     coverage, frame conversion, strand and multi-copy refusal.
   - One nearest-seam helper for IS6110.
   - One key normaliser.

   Findings 3.2, 4.4, 4.8, 5.1 and 5.2 are each a script getting one of these
   wrong on its own.

4. **Tie "done" markers to their inputs.**
   - Write outputs to `.tmp` and rename them last.
   - Record the build ID, the reference ID, and the keys or matrix hash in each
     per-sample output.
   - Skip a sample only on a match.

   This fixes 3.6, 4.6, 6.5, and the case where P2-P4 keep outputs made against
   an old reference after P1 is rerun.

5. **Make every writer handle empty input.** Each should emit a header-only
   file, not crash, not write `rows[0]`, and not write nothing. Readers should
   treat header-only as "measured, nothing found".

6. **Harden job submission.**
   - `die` on any empty upstream job ID.
   - Declare each pass's inputs next to its dependencies, so 6.2 cannot recur.
   - Use one helper for cohort rows, array size and task-to-sample mapping.
   - Build shared per-build indexes in P0 only.

7. **Add a small regression suite on synthetic data.** Include:
   - a reverse-complemented reference, for P4 indel re-anchoring and P4b
     breakpoint order;
   - nested mask intervals;
   - a round-trip of both IS6110 key frames through stage 2;
   - a sample with zero calls;
   - a multi-allelic site where the sample carries the other allele.

   Run it before any cohort run.

8. **Fix the release artifact.**
   - Add the missing modules from section 1.
   - Move site paths into `site.local.sh.example`.
   - Fix the config ordering (2.1, 2.2).
   - Correct the `--cohort` documentation (6.3).
   - Re-run the dependency closure with Python imports included.

9. **Regenerate affected outputs.** gwas1000 and scale200 are affected by 3.1,
   3.3, 4.1, 4.3 and 5.1, and scale100 and pilot by 4.1. Re-run P4 onward after
   the fixes, and diff the new cell counts against the old ones before any
   downstream association or convergence analysis.

---

## Suggested fix order

| order | items | why |
|---|---|---|
| 1 | 1.1, 2.1, 2.2, 6.3 | makes the repository runnable as published |
| 2 | 4.1, 4.2, 4.3, 3.1, 3.4 | largest wrong-cell counts in existing VCFs, small code changes |
| 3 | 5.1, 5.2, 3.2, 3.3 | coordinate and state correctness in the core arms |
| 4 | 6.1, 6.2, 5.3, 6.4 | orchestration failures that produce stale or missing data |
| 5 | recommendations 1, 2, 4 | stop this class of bug from returning |
| 6 | everything else | medium and low items |

---

## Fix status (2026-09-29)

All changes are uncommitted in this checkout. Test inputs were read from the
working tree at `/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome`, which
was not modified. Test outputs are in `test_runs/`, which is gitignored.

Run the regression suite from the repository root:

    source config/project_env.sh && $MTB_PY tests/run_tests.py

It has 16 tests. All pass on the fixed code, and 14 fail on commit `6ab0f8d`.

### Fixed and verified on real data

| item | change | evidence |
|---|---|---|
| 1.1 | Copied in `mtb_norm.py`, `p5_states_io.py`, `graph_frame.py`, `graph_frame_offsets.py`, `p1i_discover_matched.sh`, `show_config.sh`, `stage_in.sh`, `build_repeat_mask.py` | Tracked code files are byte-identical to the working tree |
| — | New P0 step `--step frames` writes the frame table into the build; the runner exports it as `MTB_GRAPH_FRAMES` | The table was a hidden input that no step produced |
| 2.1, 2.2 | Site file sourced first, found relative to the config file | Derived paths now follow `MTB_WORK`, from any directory (test) |
| 3.1 | `*` added to any record with a GT=2 cell | scale200: 976 records fixed; `bcftools +fill-tags` passes |
| 3.5 | Symbolic records take REF from H37Rv; deletions use the padding base and END | scale200: REF mismatches drop from 7,291 to 0; region queries inside deletions work |
| 3.4 | Pilot key-table fallback removed; sample columns restricted to the P5 set | Record counts are unchanged on scale200 |
| 3.2 | Node-frame REF read at `start + off`, converted panel to refs | Base agreement: 56.9% to 100% on `+`, 15.3% to 99.8% on `-` (see correction in 3.2) |
| 3.3 | Only GT=0 gVCF lines count as coverage; another allele at the key's position blocks REF; SNP reversions stay REF | About 0.6% of H37Rv REF cells per pilot sample become NOCALL; the causes were checked by type |
| 4.1, 4.2 | Stage-2 node keys keep their offset; "not measured" is NOCALL; odgi and samtools failures are fatal | Pilot node frame: 1,061 blanket ABSENT become 700 REF, 46 NOCALL, 315 ABSENT; the H37Rv frame is unchanged |
| 4.3 | `--all-stacks` passed; the writer refuses empty keys | Cohorts run with the flag have 0 empty keys; gwas1000 had 1,800 without it |
| 5.1 | Mask intervals merged; 1-based test; class taken from the prefix | Brute force over all 4.41 Mb: 0 mismatches |
| 5.2 | SV direction from projected geometry; spans checked | On the flipped pilot reference, reversed coordinates on placed events drop from 67 to 0. Forward samples change only through the span check. |
| 4.8 | Copy-number window mirrored on the minus strand | H37Rv flank found on the reverse side for 595 of 651 minus-strand loci, and on the forward side for none |
| 4.7 | Accessory presence checks every subprocess and refuses an empty pool | Output byte-identical to production on a real sample; a bad CRAM now fails |
| 6.1, 6.2, 6.3, 6.7, 6.8 | Runner: sbatch rejection stops the chain; p5pre waits for p4b; `--cohort` accepted; blank-line-safe row count; one completed build or `MTB_BUILD_DIR` | Dry runs and a fake sbatch that rejects the third submission |
| 6.10 | P1 and P2 summaries exit 2 on missing samples (`--allow-missing` overrides) | Tested on the pilot |
| 6.11 | `--status` uses globs, and every pass has a case | Tested on the pilot |
| 3.8 | Nested genes found; pseudogenes included | Brute force over all 4.41 Mb: 0 mismatches |
| rec. 1 | `bin/vcf_gate.sh` runs at the end of `p5_finish.sh` | Fails the existing scale200 VCF on 3 checks; passes the fixed one |

### Fixed, syntax-checked, and not run live

| item | change | why not live |
|---|---|---|
| 5.3 | Undecidable frame with no records is let through in P4 and P4b | Needs a full P4 task |
| 5.4, 5.5 | P4b keeps the inherited half for zero-call samples; the graph VCF is passed and required in P4 and P4b | Unit test only |
| 3.6 | P5 skips only on a matching `keys_sha1`; P5svgt only when newer than the probes | Checksum confirmed on scale200 |
| 4.6 | Junction table always written; P1i and P1g skip only when all outputs exist | Needs an isolate with zero sites |
| 5.6 | `sv_twoframe.py` fails on samtools errors | Needs a live task |
| 6.4 | P1 and P2 use `MTB_CRAM_REF` and require the CRAM root | Same file as before for existing cohorts |
| 6.5, 6.6 | Atomic FASTQ and SV VCF writes; a dysgu failure is fatal; `.p2.done` marker; `SIM_REQUIRE_FQ` guard | Needs a live P2 task |
| 6.9 | `.bwt` installed last in P0, P3 and `simulate_and_call.sh`; P3 waits for a complete index | Race conditions only |
| 4.9, 4.10 | Element-depth span check (no current rows affected); missing crossmap is fatal (every current reference has one) | Checked against existing outputs |
| 5.7 | P4 counts drops; the off-the-end check needs H37Rv; P4b summary rates count called rows only | P4b summary tested on the pilot |
| 3.10 | Sanity outliers are symmetric | — |
| 2.3, 2.4 | Tool variables declared; script fallbacks removed; `config/site.local.sh.example` added | — |
| 6.11 | `sync_back.sh` exits nonzero on a failed push or copy; `SYNC_BACKUP=1` keeps overwritten files | Writes to durable storage, so not run |

### Deliberately not changed

Each of these changes scientific output or existing identifiers, so it is a
decision for the project rather than a mechanical fix.

- **3.7, ancestral alleles at the ingroup root.** This changes the AA method.
- **3.9, two-frame confirmation in caller-matrix mode.** This changes SV genotypes.
- **4.4, seam slop.** The scripts disagree, with 0 bp in two and 25 bp in one, and the right value needs choosing.
- **4.5, clustering IS6110 keys within about 10 bp.** This changes the key space.
- **5.7, SVLEN after interval merging.** This would rename every `svi:` interval ID and orphan the existing genotype tables.
- **6.11, stable tie-breaking in `t8_select_reference.py`.** The current sort is already deterministic, and a stable sort would change the chosen reference for tied samples.

### Decisions on the deliberately unchanged items (2026-10-01)

| item | decision | where |
|---|---|---|
| 3.7 ancestral alleles | fixed: AA is the MTBC ancestor's state; 4,759 of 5,130 ties resolved, no resolved site changed | `bin/ancestral_alleles.py`, P0 `--step ancestral` |
| 3.9 caller-matrix inherited calls | fixed: NOCALL without a second frame, as in interval mode | `bin/p5_sv_genotype.py` |
| 4.4 seam slop | fixed: one 3 bp rule in promotion, flank placement and the writer | `is6110/bin/is6110_seam.py` |
| 4.5 split IS6110 keys | fixed: keys within 6 bp across isolates share a canonical position | `is6110/bin/is6110_write_vcf.py` |
| 5.7 SVLEN after merging | fixed: SVLEN is the merged span | `bin/sv_intervals.py` |
| 6.11 reference tie-break | left unchanged: changing it would force realignment of every cohort | — |

Found while settling them, and fixed: the SV interval arm (catalogue, two-frame
tables, interval-mode genotyping) and accessory presence ran outside the chain,
gwas1000 used scale200's catalogue, and the merge dropped every caller deletion
when a catalogue was present. See HANDOFF.md section 2.

### Behaviour changes to expect on re-run

- **P4.** More masked bases are routed to the composed arm. On two pilot samples, direct core calls fell by 46 to 49, and composed PE/PPE calls rose by a similar number.
- **P4b.** Some events whose projected span disagrees with the called span become breakends. This affected 1 to 25 per pilot sample, mostly in PE_PGRS and IS6110 regions.
- **P5.** About 0.6% of H37Rv-frame REF cells become NOCALL where the sample carries another allele.
- **IS6110 stage 2.** Node-frame cells gain REF and NOCALL calls in place of false ABSENT calls.
- **Merged VCF.** Deletion POS moves back one base, REF is a real base, and END is present.

Existing cohort outputs, including the pilot, scale100, scale200 and
gwas1000, were produced by the old code. Regenerate them from P4 onward, and
from p1iv for the IS6110 arm, before any downstream use.
