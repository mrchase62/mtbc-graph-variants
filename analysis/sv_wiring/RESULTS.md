# Our caller in production: real-data check (2026-10-12)

Code: commit 6f400f3 (P2 `_ours`, P4b rules and merge, P5 unsized INS,
`bin/sv_rules.py`). Check: `01_fbrecall_check.sbatch` (job 52446924, 23 min
on 8 cores) and `compare.py`. Report: `out/compare_report.txt`.

Cohort fbrecall: 48 isolates with finished P2 and kept BAMs. P2 ran only the
new step on each kept BAM (the path for samples finished before the wiring).
P4b ran into `refbias/fbrecall/run/p4b_test`. 27 of 29 masks were copied
into `refbias/callable_masks/7713a8d71d8e-fix1/` from the evaluation; 2 were
built by P2.

## 1. Reproduction: exact
- Caller tables identical to `recall/out/proto`: 48 of 48.
- Events after assembly and the AF floor identical to
  `recall/out/asm_af/keep_na_0.1`: 48 of 48.
- The 2 masks built by P2 are identical to the evaluation's.
- tatC-helY genotypes: 28 present, 1 absent, 17 reference lacks it,
  2 unresolved.

## 2. P4b rows scored (called, FILTER PASS; callable sequence)

| group | set | typed recall | calls | strict precision | false calls |
|---|---|---|---|---|---|
| kept, without GCF_035581225 (26) | evaluation, both callers summed | 0.885 | 179 | 0.838 | 29 |
| | P4b, callers merged | 0.885 | 117 | 0.761 | 28 |
| excluded (12) | evaluation | 0.863 | 98 | 0.724 | 27 |
| | P4b | 0.863 | 72 | 0.625 | 27 |
| chimeric controls (9) | evaluation | 0.859 | 2,112 | 0.053 | 2,000 |
| | P4b | 0.859 | 2,015 | 0.036 | 1,942 |

- **Recall is unchanged in every group, and so are the false calls.**
- **Precision is lower in P4b only because it merges.** The evaluation
  scored the two callers' calls as a sum, so a true event both callers
  found counted as two correct calls. P4b makes it one row. False calls are
  mostly single-caller and stay as they were.
- **The per-event figure is the one to quote from now on:** strict
  precision 0.76 in the kept group, about 1.1 false calls per isolate.
  Earlier combined precisions (0.785 etc.) counted shared events twice.
- The evaluation's kept-group 0.838 here, against 0.785 before, adds the
  tatC-helY genotype rule; earlier fair comparisons used the HC hand-off
  only.
- Rules applied in P4b over the 48 isolates: a handful of element and HC
  drops per sample at most (`logs/p4b.*.log`).

Not checked: P5 on these rows (the unsized-INS clustering is covered only by
`tests/test_ours_wiring.py`).
