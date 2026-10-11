# Before scale-up: stricter fold-back QC and call harmonization (plan, 2026-10-11)

User request: "QC the samples with a more rigorous fold-back chimera scan and
harmonize the calls from HaplotypeCaller, SV and IS6110 mapping."
Nothing below has been run.

---------------------------------------------------------------------------
## A. Fold-back chimera scan, version 2

### What a fold-back chimera is
During library preparation a DNA fragment can fold back on itself and be
copied into a molecule that reads forward and then back along the
reverse-complement strand of the same stretch.

Aligned to the reference, such a read runs forward for a while. Then it either
- continues as a second piece on the opposite strand a short distance away
  (a split read with an SA tag), or, if that second piece is too short to
  align,
- ends in a soft clip whose sequence is the reverse complement of nearby
  reference.

### What the current scan (`controls/split_reads.awk`) does and misses
It counts reads whose first SA entry is on the same contig, on the opposite
strand, within 1 kb, per 1,000 primary alignments. Its limits:

1. **Split reads only.** Short fold-backs (second piece under ~20-30 bp) get
   no SA tag and are missed. On the synthetic chimeras (foldsim, 5% of pairs
   replaced) it found 10.9 of 15.8 junction reads per 1,000, about 70%.
2. **Only the first SA entry** is examined.
3. **No mapping-quality filter**, so reads in repeats count, and those
   alignments are unreliable.
4. **The same rate means different things for different read lengths.**
   Longer reads have more room to show a fold-back (160-400 bp reads: 72% of
   isolates above 1 per 1,000; under 110 bp: 41%).
5. **No separation from places where the reference itself folds**
   (inverted repeats and hairpins, e.g. the tatC-helY element). There, the
   reference sequence itself encourages fold-back molecules, and they are
   counted the same way as random library chimeras.

### How it compares with ArtifactsFinder (github.com/lilicai/ArtifactsFinder)
ArtifactsFinder does not measure samples. It scans a reference (hg19/hg38)
for:
- inverted repeats: exact reverse-complement arms of at least 8 bp combined,
  a spacer of at least 5 bp, 50 bp padding;
- palindromes: one mismatch allowed.

It writes a blacklist of positions where library artifacts create false
small variants. That idea is point 5 above: artifacts concentrate at
reference inverted repeats. It does not replace a per-sample rate, and as
written it supports only human genomes. Version 2 therefore keeps the
per-sample rate and measures the reference-driven sites directly from the
data (step 5).

### Method, step by step
1. **Input:** each sample's P2 alignment to its matched reference (BAM or
   CRAM).
2. **Reads counted:** primary alignments, mapped, MAPQ >= 20, not
   duplicates. Their number N is the denominator.
3. **Split fold-backs:** any SA entry, not just the first, with MAPQ >= 20, on
   the same contig, the opposite strand, and starting within 1 kb.
4. **Clip fold-backs:** reads without such an SA, with a soft clip of >= 15 bp
   left after quality trimming. Trimming is the caller's rule: stop at the
   first 5-base window averaging below Q20.
   - The clip counts when the reverse complement of its first 15 bases
     occurs exactly in the reference within 500 bp of the clip point.
   - For context, clips matching the forward strand (tandem-like) and clips
     matching nothing nearby are counted separately.
5. **Sites where the reference folds:** fold-back junctions are tallied by
   position across all samples on the same reference. A position where
   fold-backs occur in more than 20% of samples is a property of the genome,
   not of one library. Those sites are listed, and each sample's rate is also
   reported without them.
6. **Per-sample output:**
   - split, clip and total fold-backs per 1,000 reads;
   - the rate without recurrent sites;
   - duplicate rate, read length, median insert size;
   - the share of all clipped reads that are fold-backs.

### Validation (before any cohort run)
- **Synthetic chimeras** (existing foldsim data; QC_6 with fold-back pairs
  injected at known rates, inverted-repeat and random models, with a list of
  every inserted junction):
  - sensitivity of version 2 vs the current scan at each rate;
  - false fold-backs at rate 0.
  - Goal: >= 90% sensitivity and ~0 at rate 0.
- **Phase B real controls:**
  - separation of the 10 chimeric controls (old rate 4.8-56) from the 15
    clean ones;
  - how well each version predicts a sample's false SV calls (dysgu PASS and
    ours, from the existing Phase B outputs), compared by rank correlation.
- **Read length:** check whether version 2's rate still depends on read length.

### Cutoff
The current cutoff (2 per 1,000) was set on the old metric. A new cutoff
would be proposed from the validation, for the user to decide. Then the scan
runs on the cohorts in use (1,200 isolates) and Phase B.

### Cost
- Validation: the synthetic BAMs plus 62 Phase B BAMs, about 0.5 CPU-hours.
- Cohorts in use: first timed on 5 samples. The estimate is 10-20 CPU-hours,
  since Python is slower than the awk scan; reading the first 2 M reads per
  sample would cut this to ~3 CPU-hours, and the rate is a rate.
- No downloads.

---------------------------------------------------------------------------
## B. Harmonizing HaplotypeCaller, SV and IS6110 calls

### The problem (measured on scale200_fix, `analysis/harmonize/overlap.py`)
Production has three arms that never look at each other.

- **Small variants:** HaplotypeCaller (P2) -> P4 -> P5.
- **SVs:** delly + dysgu (P2) -> P4b -> P5. One deletion can be written up to
  three times: by P4's composed arm (`p4_place.py`, "inherited half"), by
  P4b's inherited records (`p4b_place_sv.py`), and by the graph deletion
  catalogue (`sv_intervals.py`).
- **IS6110:** IS6110 arm -> P5.

The only cross-arm check is in `merge_cohort_vcf.py` `dup_of()`. It drops a
caller SV insertion only when its length is 1,340-1,370 bp and it lies within
20 bp of an IS6110 record. Truncated IS6110 calls (median 234 bp) survive it.

Measured duplicates:
- **787 SV records (13%) repeat a small-class indel.** Mostly these are
  inherited reference-vs-H37Rv differences, written once as exact alleles by
  P4 and again as SV inherited or interval records. Called SV vs
  HaplotypeCaller is only 3%.
- **291 caller SV insertions (14%) are IS6110 insertions** the IS6110 arm
  already has.
- **IS6110 copy losses are not reported by the IS6110 arm.** Its `<INS>`
  allele cannot express a loss; `is6110_p5_stage2.py decide()` leaves those
  isolates NOCALL. The ~1.36 kb SV deletion intervals are the only record of
  a loss.

### Rule: each event is written once, by the arm best placed to call it
1. **Indels under 50 bp belong to HaplotypeCaller.**
   - A caller SV call with a HaplotypeCaller PASS indel of 10-49 bp within
     60 bp in the same sample is dropped. This is the Phase B hand-off rule
     (HC depth >= 50% of the isolate's median, alt >= 90%). On Phase B it
     removed 11 false calls and no true ones.
2. **Events of 50 bp or more belong to the SV block,** so there is one key
   space for them.
   - When a small-class allele of >= 50 bp and an SV record describe the same
     event (start within 60 bp, sizes within 70%, the same carriers), they
     merge into one record. The SV record is kept, and the small-class
     allele's exact sequence goes into its INFO.
   - P4b's inherited records and the graph deletion catalogue describe the
     same reference differences. Only one is kept: the catalogue record,
     since it carries the genotyping.
3. **IS6110 insertions belong to the IS6110 arm.**
   - A caller SV insertion of any size within 100 bp of an IS6110 record in
     the same sample is dropped.
   - This widens `dup_of()`, whose current window is 20 bp and 1,340-1,370
     bp only.
4. **IS6110 copy losses get one clear record.**
   - An SV deletion whose ends lie within 50 bp of the ends of an H37Rv
     IS6110 copy is relabelled CLASS=is6110, as a copy loss.
   - Phase B showed short reads cannot reliably separate a real copy loss from
     a false one: 35 real vs 16 false, and depth does not help. So these
     records are flagged, and kept out of the robust tier until the tree
     check says otherwise.
5. **Every final record lists the arms that saw it** (`SOURCES=`). That is
   cross-support, kept as information instead of duplicate records.

### Where to make the change
In the merge step (`merge_cohort_vcf.py`), plus the sample-level small-vs-SV
check (rule 1), which reads P5's existing states.

- All the needed information is already there: COMPONENT, SITECLASS, the
  IS6110 tables and the states.
- Rerunning only that step on scale200_fix costs about 1 billing-hour.
- Changing P4 or P4b instead would also rerun P5 states (362 billing-hours
  on scale200_fix).

### Verification
- **Unit tests** in `tests/` for each rule. No test covers small-vs-SV
  duplicates or IS6110 copy losses today.
- **Rerun merge on scale200_fix:**
  - overlap.py should report ~0 duplicates;
  - the number of records removed or merged by each rule;
  - the drug-resistance positive controls unchanged;
  - every remaining SV record still reachable from a sample call.
- **Before and after diff** of `p5/merged.vcf.gz`, reported, with no other
  production rerun.

### Not decided by this plan
- **Which SV caller production should use.** The Phase B caller's rules are
  not in production; production uses delly + dysgu. These harmonization rules
  do not depend on which caller is used.
- **The robust-tier filter itself** (callable sequence, DEL/INS, both
  breakpoints). It is a separate step after harmonization.

---------------------------------------------------------------------------
## Suggested order
1. Fold-back version 2: write it, then validate on synthetic and Phase B
   data (~0.5 CPU-hours), and report.
2. Harmonization: code and tests, then the merge-only rerun on scale200_fix
   (~1 billing-hour), and report.

The two are independent and can proceed side by side.

## Decisions for the user
1. Approve fold-back version 2 validation (~0.5 CPU-hours)?
2. Approve the harmonization rules 1-5 and changing production code
   (`merge_cohort_vcf.py`, tests) with a merge-only rerun of scale200_fix
   (~1 billing-hour)?
3. IS6110 copy losses: relabel and flag as in rule 4, or leave them out of
   the output until they can be called reliably?
