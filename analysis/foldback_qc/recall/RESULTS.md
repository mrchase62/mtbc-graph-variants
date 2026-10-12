# Part D: does fold-back hide real SVs? (2026-10-12)

Plan: `../FOLDBACK_SV_EFFECT_PLAN.md` part D.
- Report: `out/analyse_report.txt`.
- Tables: `out/recall_by_group.tsv` and `out/recall_by_isolate.tsv`.

## Design
- **Cohort `fbrecall`:**
  - 39 selfcal54 isolates that passed QC, plus the 9 chimeric controls
    (all reads).
  - Each isolate's QC'd reads are mapped to the nearest OTHER panel genome
    (`make_cohort.py`). That is the closest by panel SNP-profile distance,
    excluding the own genome and genomes within 50 SNPs. Distances are
    51-2,152 SNPs; GCF_035581225 is at 11,696.
- **Cohort `fbrecall100`:** the 9 chimeric controls with every fold-back
  pair removed, on the same reference and against the same truth. This is
  the causal arm.
- **Truth:** every >= 50 bp difference between the own genome and the
  reference (`novel_events/truth.py`), 1,608 events. Scoring uses callable
  sequence only, with strict precision (`score.py`).
- **Groups by v2 rate:**
  - kept: v2 <= 4.9, 27 isolates;
  - excluded: v2 > 4.9, 12 isolates;
  - chimeric: 9.
- **Cost:** about 32 billing-hours.

## Results, pooled by group (callable sequence)

| group | true events | caller | typed recall | typed calls | precision (strict) |
|---|---|---|---|---|---|
| kept | 182 | ours + assembly | 0.797 | 271 | 0.557 |
| | | dysgu AF >= 0.2 | 0.687 | 225 | 0.578 |
| | | combined | 0.841 | 496 | 0.567 |
| excluded | 51 | ours + assembly | 0.784 | 150 | 0.280 |
| | | dysgu AF >= 0.2 | 0.569 | 57 | 0.526 |
| | | combined | 0.863 | 207 | 0.348 |
| chimeric | 78 | ours + assembly | 0.859 | 774 | 0.087 |
| | | dysgu AF >= 0.2 | 0.577 | 1,768 | 0.027 |
| | | combined | 0.872 | 2,542 | 0.045 |
| chimeric, fold-backs removed | 78 | ours + assembly | 0.910 | 105 | 0.676 |
| | | dysgu AF >= 0.2 | 0.705 | 92 | 0.620 |
| | | combined | 0.923 | 197 | 0.650 |

Per isolate, recall does not track the fold-back rate (Spearman -0.15,
p 0.32, 48 isolates).

## What it shows
1. **Fold-backs cost precision, not much recall.**
   - Combined recall is about 0.84-0.87 in every group.
   - Combined precision falls from 0.57 (kept) to 0.35 (excluded) and 0.045
     (chimeric).
   - dysgu at AF >= 0.2 loses the most recall: 0.69 in the kept group,
     0.57-0.58 in the excluded and chimeric groups.
2. **Causal arm: fold-backs hide a few real events.** On the same
   libraries, references and truth, removing the fold-back pairs changes:
   - combined recall 0.872 -> 0.923 (68 -> 72 of 78 events);
   - dysgu AF >= 0.2 recall 0.577 -> 0.705 (45 -> 55);
   - our caller's recall 0.859 -> 0.910;
   - combined precision 0.045 -> 0.65.

   The numbers are small: 78 callable events, 1-35 per isolate. N1176 has 2
   callable events, missed with fold-backs and found without.
3. **How our caller loses events.**
   - Our caller loses them before assembly. Pre-assembly typed recall is
     0.859 with fold-backs and 0.923 without.
   - Breakpoint recall is 0.987 in both cases.
   - So the real breakpoints are still found, but nearby fold-back reads
     change how the event is typed.
   - The 40-window assembly cap is not the cause.
4. **Excluding by the cut-off keeps recall and gains precision.** The
   excluded selfcal isolates have recall similar to the kept ones but half
   the precision.

## Open: precision in the kept group is lower than in Phase B
Kept-group precision is 0.56. Phase B's novel set reached about 0.8-0.83,
and these same isolates gave 1-2 false calls each on their own genome.
- **Distance from truth:** of our caller's ~120 false typed calls, 92 lie
  more than 500 bp from any truth event (54 DEL, 38 INS), so they are not
  placement offsets.
- **Shared references:** 15 of the 92 recur at the same reference position
  in 2 or more isolates mapped to the same reference (9 of 12 on
  GCF_900520315). Those point to the reference or to truth gaps.
- **The rest:** not explained yet. The candidates are:
  - reads mapping to a divergent reference;
  - events missed by assembly-vs-assembly alignment in truth.py;
  - assembly errors in either genome.

This affects the precision of the production SV caller set, so it should be
checked before that set is fixed.

## Not done
- **An AF floor on our caller** (parts A-B suggested 0.2). `score.py` has no
  AF filter for our caller; testing it needs per-call AF (as
  `../fb_effect.py`) applied to the asm tables.
