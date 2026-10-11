# Fold-back cut-off calibration with panel genomes' own Illumina reads (step 2b, not launched)

Status 2026-10-11: proposed, held by the user to check what it gets us.

## The question it answers
The fold-back cut-off decides which samples are excluded as chimeric. Today
the evidence is two groups of controls with nothing between them:
- clean controls: <= 1.55 fold-back reads per 1,000 (v2 scale), few false calls;
- chimeric controls: >= 6.45, many false calls.

The proposed cut-off of 2.5 sits in that gap without data to support it. In
the 500 random production samples, about 19% fall at 1-5 per 1,000 (old
scale), so the cut-off decides whether roughly 1 in 5 samples is kept.

## What we would do (per isolate)
1. Download the full Illumina run of the same BioSample as the finished panel genome (ENA, md5-checked).
2. Run the standard read preparation (08_qc.sbatch, 09_kraken.sbatch):
   - fastp (the user's settings);
   - TB-Profiler mixed-sample check;
   - rasusa to 100x, with a 60x minimum;
   - kraken2.
3. Run P2 with the reference pinned to the isolate's own panel genome (as marinSyn25): bwa, HaplotypeCaller, dysgu, delly.
4. Run the fold-back scan v2 on that P2 BAM, which holds the QC'd, downsampled reads.
5. Run our SV caller v6 plus local assembly, then self_noise.py. Every call against the isolate's own genome counts as noise, split callable / masked.

## What we get
- Per isolate: a fold-back rate and false-call counts (HC SNPs, HC indels, our typed SV calls, dysgu PASS and AF >= 0.2) in callable sequence.
- A curve of false calls against fold-back rate, filling in the 1-6 gap if isolates land there.
- The own-genome noise estimate, now from 11 isolates of one study, extended to about 45 isolates across about 15 studies and 8 instruments. This is a baseline for the robust tier too.

## Candidates
- 54 paired-end genomes, after removing Phase B controls and a single-end set; 75 Gbases.
- Raw depth:
  - < 70x: 10. These will probably fail the 60x minimum after QC: 5 PRJEB49093 HiSeq 2500 runs, 2 NextSeq 550, 2 MiSeq and 1 NextSeq 500.
  - 70-100x: 4.
  - 100-200x: 13.
  - >= 200x: 27.
- About 44 are expected to pass QC.
- Three runs are very large and would be downsampled after download:
  - KZN MDR: 22 runs, 12 Gb, GA II/HiSeq 2000 (2009);
  - ushi10: 9.2 Gb;
  - TCDC11: 5.2 Gb.

## Limits (what it cannot tell us)
1. **Unknown how many land in the gap.** Instrument mix (MiSeq 22, HiSeq 2500 8) suggests several, but it could be only a handful; then the curve stays thin.
2. **Polishing bias.** Many hybrid assemblies were polished with these same Illumina reads. That forces the genome to agree with the reads at small variants, so HC noise is underestimated. SV-scale noise is less affected but not immune.
3. **Genome and reads from different DNA.** This happened with M0010874_7 and TB3251 in Phase B. Such isolates show many "false" calls that are real differences; they must be flagged (identity check), not read as chimera effects.
4. **Fold-back confounded with library type.** MiSeq/HiSeq 2500 libraries are also older, with different read lengths and inserts. A rise in false calls could come from library era, not chimeras. With ~44 isolates this can only be checked roughly (by study).
5. **Callable sequence only.** Masked-region noise is reported but not used for the cut-off.
6. **Scale.** The rate is measured on the own-genome BAM, while production rates come from H37Rv-aligned CRAMs. Fold-backs are local, so the rates should agree; checkable on the 2 isolates also in production.

## Cost
- Downloads: about 75 Gbases, roughly 35-40 GB compressed, to netscratch (as Phase B).
- Compute, from Phase B measured costs, per isolate: QC 1.2, kraken2 0.2, P2 0.4, callers 0.15.
- About 2 billing-hours per isolate, so about 100-110 billing-hours for 54. Isolates failing QC stop after QC.
- Option: drop the 10 under 70x before download (saves about 2.5 Gb and about 12 billing-hours), and KZN MDR (2009 GA II era).

## Added 2026-10-11 (user): rate plus insert size
Question: does a rule using the fold-back rate and the insert size keep more
samples than the rate alone, for the same false-call level? Prompted by
RW_TB008 (10.9 per 1,000, long inserts, 2 false calls) against N0145 (13.0,
95 calls).

Fixed before looking at the selfcal54 results:
- Points:
  - selfcal54 isolates passing QC;
  - the 11 clean Phase B controls;
  - the 9 chimeric controls (`removal/`, all-reads set).
  All are on QC'd reads, mapped to their own genome.
- Outcome: false calls in callable sequence, counted for our caller's typed
  calls and for dysgu at AF >= 0.2 (the current combination).
- Predictors: fold-back v2 rate. Insert size from the P2 BAM, both as the
  median and as the share of pairs shorter than 2x the read length (an
  overlapping-pair proxy for very short inserts). Library type (adapter call)
  is reported, not fitted.
- Comparison: for each rule, pick its threshold at the largest value where
  the kept isolates average at most the clean-control false-call level
  (mean + 2 SD of the 11 clean controls). Report how many isolates each rule
  keeps.
  1. rate <= c;
  2. rate <= c, or rate <= c2 when the median insert >= i.
  Rule 2's extra thresholds are fitted on few points, so its gain must hold
  in leave-one-study-out checks before it is used.
- Also: Poisson/negative binomial fit of false calls on rate, insert and
  their interaction (descriptive only).
