# What fold-back chimeras do to SV calls: study design (2026-10-11, proposed)

**Question:** how many SV calls do fold-back chimeras cause, of which types
and at which allele fractions? Do they also hide true events? How do the
answers depend on the fold-back rate, apart from library type, insert size
and study?

The user's context:
- Within-host Tn5 samples with very short inserts could not be called below
  10% allele fraction.
- In production, fold-back tracks short inserts within TruSeq libraries
  (`adapters/RESULTS.md`).

The current evidence is correlational: chimeric isolates have more false
calls. The parts below go from cheap and correlational to causal.

## A. Attribution: which false calls sit on fold-back junctions (no new compute)
Data: selfcal54, plus the 11 clean and 8 chimeric Phase B own-genome isolates
(the 8 have reference = own genome in marinB63).

Every call against the isolate's own genome is false. Steps:
1. Fold-back junction positions per isolate, from `foldback_v2.py --sites`
   run on all reads (`--frac 1.0`, changed in `selfcal/03_post.sbatch`).
2. For each false call (our caller, dysgu, delly; type, size, allele
   fraction), test whether a breakpoint lies within 150 bp of a junction
   with >= 1 fold-back read.
3. Chance level: the same test with the calls moved to random positions in
   callable sequence (1,000 shuffles). Calls attributable =
   observed - expected.
4. Report the attributable share by caller, SV type (INV/BND expected; also
   small DEL/DUP/INS?), size, and allele-fraction bin.

## B. Dose-response, with confounders as covariates (no new compute)
1. Per isolate: false calls in callable sequence against fold-back rate, by
   caller and allele-fraction bin (< 0.1, 0.1-0.2, 0.2-0.5, > 0.5).
   - dysgu: FORMAT AF.
   - Our caller: supporting reads / depth.
2. Negative binomial regression:
   false calls ~ fold-back rate + library type (Tn5/TruSeq/unclear) +
   insert-size mean + read length + depth.
   This tells whether fold-back predicts false calls once insert size and
   library are accounted for.
3. Per fold-back band, the lowest allele fraction at which false calls fall
   to the clean-isolate level. This is the minimum allele fraction that band
   can support, and it speaks to the 10% floor in the within-host study.

Limit: still observational. About 65 isolates; the bands in between may be thin.

## C. Read-removal experiment: causal, same library (needs approval)
Remove fold-back reads from a library and rerun everything. Calls that
disappear were caused by those reads. Insert size, Tn5 bias, depth and
genome are unchanged.

1. Isolates: the 8 chimeric Phase B controls (own genome; 6.45-64.7 per
   1,000), plus up to 4 selfcal54 isolates in the 1-6 band.
2. Read QC as the cohort standard (08_qc.sbatch) for the 8. They were never
   QC'd because they failed the fold-back cut-off.
3. `foldback_v2.py` gains `--names`: write the names of fold-back read pairs.
4. Four read sets per isolate, made from the QC'd unaligned BAM:
   - all reads;
   - 50% of fold-back pairs removed;
   - 100% of fold-back pairs removed;
   - the same number of random pairs removed (control for the depth change).
5. Each set goes through P2 (own genome pinned), our caller with assembly,
   and self_noise counts.
6. Report false calls per set, by caller, type and allele fraction. The
   slope over 0/50/100% removal gives calls per 1,000 fold-back reads.
   Calls left at 100% removal come from fold-backs the scan misses
   (sensitivity 0.74-0.88) or from other causes.
- Cost: 8 x QC 1.2 = 10; 12 isolates x 4 sets x (P2 0.4 + callers 0.15)
  = 26. About 36 billing-hours in all.

## D. Does fold-back hide true SVs? Recall (needs approval)
Precision alone misses the cost to sensitivity: chimeric reads could pull
breakpoints or depth signals away from real events.

1. Map each selfcal54 isolate's QC'd reads to the nearest other panel
   genome, which differs by real events.
2. Truth is the assembly-vs-assembly differences (truth.py, as Phase B),
   scored in callable sequence.
3. Score recall and precision against fold-back rate.
4. Optionally repeat on the read-removal sets from C, which gives a causal
   recall effect.
- Cost: P2 + callers for about 44 isolates = about 25 billing-hours, plus
  truth (minutes).
- Limit: few events per pair, so recall is pooled by band.

## Order
A and B run on the selfcal54 outputs already planned. C gives the causal
answer and is the most rigorous. D adds sensitivity. Recommended: A+B as
soon as selfcal54 finishes, C in parallel, D after C.
