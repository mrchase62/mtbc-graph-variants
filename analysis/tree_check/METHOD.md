# Tree check of SV calls: method (draft for review, 2026-10-11)

**Question.** Are the SV calls reliable enough to compare across samples for
convergence testing and GWAS?

**Idea.** A real SV is inherited, so isolates that are very close on the tree
should mostly agree about it. Close isolates do differ by a few real events,
but only a few. If close pairs disagree far more often than their SNP
distance allows, the disagreements are mostly calling noise.

This is report-only. It changes no pipeline code and needs no new truth set.

## Data

The check runs on cohort **scale200_fix**: 200 isolates, the fully fixed
production run. Nothing is rerun. It reads these existing files:

- **Tree:** `runroot/data/trees/scale200_fix.rooted.nwk`. It is built from the
  cohort's SNPs together with the 333 CX333 panel genomes, then pruned to the
  cohort.
- **SNP genotypes:** `scale200_fix/p5/merged.vcf.gz`, the CLASS=small SNPs.
  These give exact pairwise SNP distances.
- **SV calls:** `scale200_fix/p5/sv_matrix.tsv`. It has 5,845 SV keys, all in
  H37Rv coordinates. For each isolate, each key is ALT, REF, ABSENT (the
  sequence is not in that genome) or NOCALL. The INFO fields of the merged VCF
  add each call's class, type, source (graph or caller), tier, component
  (called or inherited), size and IS6110 class.
- **Matched references:** `scale200_fix/p1/refmap.tsv`, which reference each
  isolate was called against.
- **Masks:** H37Rv repeat and callable masks, to sort events into callable
  and masked.

## Steps

### 1. Pairwise SNP distance

For every pair of isolates, count the SNP positions where both have a call
and the bases differ. Positions where either isolate has no call are left
out.

### 2. Close pairs

Close pairs are pairs at most 12 SNPs apart. That is the usual threshold for
recent transmission in TB, so these pairs are, in practice, outbreak-like.
Two wider bands are reported as well: 13–50 and 51–200 SNPs.

From the existing tree, about 23 pairs are under ~45 SNPs. scale200_fix is
not an outbreak study, so the closest band may be small; see **Limits**.

### 3. Close-pair agreement for each SV

Take each close pair and each SV key that is ALT in at least one of the two
isolates:

- **Skipped:** NOCALL in either isolate. This is counted and reported, not
  treated as a disagreement.
- **Agree:** ALT in both.
- **Disagree:** ALT in one isolate, and REF or ABSENT in the other.

For each pair, count the disagreements and the agreements.

### 4. Yardstick: how many disagreements real events allow

How many real SV differences should a pair have? Set it by the SNP distance
and the cohort's own ratio of SVs to SNPs. That ratio comes from the tree: SV
changes per SNP change on the internal branches, which reflect older
inherited history and are less exposed to calling noise.

Plot disagreements per pair against SNP distance:

- **Real differences** grow with distance.
- **Noise** appears as disagreements in pairs 0–5 SNPs apart, roughly
  constant whatever the distance.

The number to report is that baseline, in false SV differences per pair. It
is the noise a GWAS or a convergence test would actually see.

### 5. Split by what might cause the noise

Each split reports both the baseline and the agreement rate:

- **Same matched reference vs different matched references.** About 145 of
  200 isolates share a reference with their nearest neighbour.
  - Disagreement only between isolates on different references points to a
    reference or projection artifact, not to the reads.
  - Events marked inherited (reference-vs-H37Rv differences carried over from
    the matched reference) are reported separately. They are shared by every
    isolate on the same reference, whatever the reads show, so they would
    inflate agreement.
- **Callable vs masked sequence.**
- **Event type and size:** DEL and INS, in bands of <50, 50–500, 500–5k and
  >5k bp.
- **Source:** graph vs caller, the caller's quality band, and support tier
  A/B/C.
- **IS6110 and non-IS6110.** IS6110 calls come from the pipeline's own
  IS6110 caller (CLASS=is6110) and are reported as their own group.
- **The proposed robust tier:** caller DEL/INS in callable sequence, under
  today's agreed definition.

### 6. Whole tree: how often each event arises independently

For each SV key, count the smallest number of separate gains needed to
explain its ALT isolates on the tree. This is parsimony on presence and
absence, with NOCALL treated as unknown.

- **One gain:** the event marks a single clade, as inherited events should.
- **Many gains:** either real convergence or noise.

The two can be told apart by how the ALT isolates sit on the tree:

- **Noise** is mostly single isolates whose close relatives (within 12 or 50
  SNPs) are REF.
- **Real convergence** (for example the tatC–helY element, lost about 30
  times on the CX333 tree) is gained in separate clades. Within each clade,
  close relatives agree.

The count to report is: lone-isolate calls whose relatives within 12 SNPs
all have a definite REF or ABSENT state.

This step matters directly for convergence testing, because a noisy call
looks exactly like a repeated independent gain.

### 7. Controls

- **Positive:** lineage-defining deletions (regions of difference such as
  RD105 and RD207 in lineage 2, RD239 in lineage 1). They should follow the
  tree exactly: one gain each and full agreement within lineages.
- **Noise level from today:** reads mapped back to the isolate's own finished
  genome give ~0.1 false SV calls per isolate in callable sequence, and 0–2
  SNPs in 9 of 11 isolates. The step 4 baseline should be of that order for
  the robust tier. If it is much higher, the extra noise comes from mapping
  to a different reference, from projecting calls to H37Rv, or from the
  merge across samples.
- **Isolates to flag:** any isolate with a large share of the cohort's
  disagreements is reported, not dropped. Examples of the cause are
  M0010874_7 and TB3251 today: likely mixed or a different subculture.

## Outputs (all in `analysis/tree_check/out/`)

- `pairs.tsv`: one row per pair. SNP distance, same reference or not,
  agreements, disagreements, NOCALL skipped, each by stratum.
- `events.tsv`: one row per SV key. Class, type, size, source, tier,
  callable or not, n_alt, n_nocall, minimum gains, lone-isolate calls,
  close-pair agreement.
- `summary.tsv`: the noise baseline and agreement rate for each stratum.
- A plain-language write-up in CHAT_LOG.

## Cost

One Python script reading existing files. It runs as one short job of a few
minutes, well under 0.5 CPU-hours. There are no downloads.

## Limits

- **Few close pairs.** scale200_fix has only about 23 pairs under ~45 SNPs.
  gwas1000 (997 isolates, tree exists) would give many more close pairs, but
  it has not been regenerated with the fixes. The standing rule is to
  validate on fully fixed outputs, so gwas1000 waits until it is.
- **These are production calls.** They come from dysgu and the graph as P4b
  places them, not from the Phase B prototype caller, which is not in
  production yet. The check therefore measures the current production call
  set first. Once the robust tier is in production, the same script
  measures it.
- **Agreement is not correctness.** Two close isolates on the same matched
  reference can share the same reference artifact. This is why step 5 splits
  by same vs different reference. Inherited events must be read with this
  in mind.
- **NOCALL is not REF.** If an event is NOCALL in many isolates, the check
  says little about it. The NOCALL share is reported for every stratum.

## Decisions for the user

1. **Close-pair threshold:** 12 SNPs, with the 13–50 and 51–200 bands as well?
2. **Cohort:** scale200_fix now, gwas1000 once it is regenerated?
3. **Inherited events:** keep them in the main summary or report them only
   separately? The draft reports them separately.
