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

## Kept-group precision explained (2026-10-12, `precheck.py`)
Kept-group precision was 0.567 combined, against about 0.8 in Phase B. Each
false typed call of the kept isolates was checked by aligning its
reference window (+-500 bp) to the isolate's own genome. Output:
`out/precheck.tsv`, `out/precheck_report.txt`.

1. **One isolate is most of it.** GCF_035581225 (L6) has the most distant
   reference (11,696 SNPs; everyone else 51-2,152). It gives 72 of our
   caller's 120 false calls and 64 of dysgu's 95.
   - Most sit where the reference window does not align to the own genome,
     or aligns only in part (absent, split or inverted).
   - These are divergent regions the assembly-vs-assembly truth does not
     describe at the calls' size and position.
2. **The other 26 isolates: 48 false calls of ours.**
   - 27 are no SV-scale difference between the genomes: genuine caller
     false calls, about 1 per isolate.
   - 15 are real 20-49 bp indels, below truth's 50 bp cut-off.
   - 2 are real >= 50 bp differences truth did not report.
   - 4 are complex.
3. **Part D's scoring did not use the HaplotypeCaller hand-off**
   (`score.py --hc-dir`: drop SV calls that are 10-49 bp HC indels), which
   Phase B's numbers used.

| kept group, scoring | combined recall | combined precision | ours | dysgu AF >= 0.2 |
|---|---|---|---|---|
| as in part D | 0.841 | 0.567 | 0.557 | 0.578 |
| + HC hand-off | 0.830 | 0.615 | 0.610 | 0.620 |
| + HC hand-off, without GCF_035581225 | 0.885 | 0.743 | 0.716 | 0.774 |

With those two accounted for, precision is close to Phase B's. A reference
far from the isolate inflates false calls. In production this would hit
lineages the panel covers sparsely, such as L6.
- **A per-isolate reference-distance flag is worth considering.** There is
  only one isolate above 2,152 SNPs here. RW_TB008 at 2,152 had precision
  0.79, so no threshold can be set from these data.

## AF floor on our caller (2026-10-12, `affloor.py`)
Our caller's AF is supporting reads / reads at the start position, as
`../fb_effect.py` computes it. 479 of 4,789 events (assembly-typed, support
0) have no AF. dysgu stays at AF >= 0.2. Output: `out/affloor.tsv`,
`out/affloor_report.txt`.

**Dropping events without an AF costs too much.** Our caller's recall
falls from 0.80 to 0.46 (kept). Any floor must keep them.

Keeping events without an AF, combined (ours + dysgu AF >= 0.2), recall /
precision:

| our floor | kept | excluded | chimeric | chimeric, fold-backs removed |
|---|---|---|---|---|
| none | 0.841 / 0.567 | 0.863 / 0.348 | 0.872 / 0.045 | 0.923 / 0.650 |
| 0.1 | 0.841 / 0.585 | 0.863 / 0.617 | 0.859 / 0.053 | 0.910 / 0.685 |
| 0.15 | 0.835 / 0.586 | 0.863 / 0.664 | 0.859 / 0.058 | 0.910 / 0.689 |
| 0.2 | 0.830 / 0.587 | 0.863 / 0.663 | 0.859 / 0.059 | 0.897 / 0.687 |

Our caller alone at 0.2: kept 0.797 / 0.557 -> 0.742 / 0.595; excluded
0.784 / 0.280 -> 0.745 / 0.830.

- **A floor of 0.1 on our caller costs no combined recall** in the kept or
  excluded groups. It nearly doubles precision in the excluded selfcal
  isolates (0.35 -> 0.62), bringing them close to the kept group.
- **Higher floors** trade recall for little further gain.
- **Chimeric controls barely move.** Their false calls are dysgu's.
- **Implication for the cut-off:** with the floor, the remaining fold-back
  damage in moderate libraries is mostly dysgu's. The cut-off (old <= 2.5)
  still protects against hotspot libraries. Whether a floor would justify a
  looser cut-off rests on only 12 excluded isolates (6 from one 2x300
  study), so it is not tested here.

**Adopted (user, 2026-10-12).** The 0.1 floor, keeping events without an AF,
is now `assemble_candidates.py --min-af 0.1` (function `af_floor`; default 0
= off). It reproduces `out/asm_af/keep_na_0.1/` exactly on all 57 read sets.
