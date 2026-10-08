# Novel structural events, Phase A: first results (2026-10-08)

Plan: `PLAN.md`. The user approved Phase A at 30 billing-hours. Used: **24.5**
(P1 10.1, P2 13.9, simulation 0.1, caller 0.4).

## What was run

1. **Genomes** (`genomes.tsv`): 40 external assemblies, all single-contig, with
   no long-read QC flags and not the same isolate as a CX333 genome.
   - By lineage: L1 8, L2 11, L3 4, L4 15, La3 1, L6 1.
2. **Reads** (`01_simulate.sbatch`): wgsim, 150 bp pairs, fragment 450 +/- 50,
   0.2% error, about 80x.
3. **P1 + P2:** unchanged, through `bin/refbias_run.sh novelA40` on build
   `7713a8d71d8e-fix1`. The registry row is `novelA40`; outputs are in
   `refbias/novelA40/run`. Matched references are 29-944 SNPs from the
   assemblies (median about 130).
4. **Truth** (`truth.py` -> `out/truth.tsv`): each assembly aligned to its
   matched reference with minimap2 `-cx asm5`. 993 events of 50 bp or more:
   - indels inside alignments, from paftools call;
   - events between alignments.

   | type | count |
   |---|---:|
   | DEL | 441 |
   | INS | 471 |
   | REPL | 73 |
   | INV | 6 |
   | REARR | 2 |

   88 insertions are IS6110-sized. Events at the circular origin are not
   included.
5. **Callers, scored with `score.py`:** a call matches within 50 bp, or by
   50% reciprocal overlap for deletions.
   - the prototype `breakpoint_caller.py` (`02_call.sbatch`);
   - dysgu and delly from P2, all calls and PASS only.

## Scores (all 993 events)

| caller | typed recall | breakpoint found | precision |
|---|---:|---:|---:|
| prototype v1 | 0.30 | 0.59 | 0.63 |
| prototype v2 | 0.33 | 0.56 | 0.62 |
| dysgu PASS | 0.51 | 0.53 | 0.65 |
| delly PASS | 0.26 | 0.28 | 0.49 |
| dysgu + delly PASS | 0.58 | 0.59 | 0.59 |
| prototype v2 + dysgu PASS | 0.58 | 0.65 | 0.64 |

- **v2** adds three rules: tandem duplications typed as insertions, deletions
  joined by depth alone, and clusters within 200 bp of the circular origin
  dropped.
- **Precision is a lower bound.** A typed call is counted false when no truth
  event is within 50 bp, and the truth set misses events at the origin.

## By class

- **IS6110-sized insertions (88):**
  - dysgu PASS 0.90, delly 0.05;
  - the prototype finds the breakpoint in 0.94 but types only 0.46.
- **Insertions (471):**
  - dysgu 0.54, delly 0.09;
  - the prototype finds 0.62 of the breakpoints but cannot size them.
    Local assembly is the step that would.
- **Deletions (441):** dysgu 0.50, delly 0.47, prototype 0.42; together 0.60.

## What every caller misses (v1: 313 events, `out/missed_v1.tsv`)

| class | events |
|---|---:|
| repeat: 30% or more of reads at MAPQ < 20 | 136 |
| no clipped reads, depth changed | 65 |
| no clipped reads, depth normal at the event point | 74 |
| clipped reads present but not called | 38 |

- **The "depth normal" group is mostly tandem copy-number change measured at
  the wrong point.** Example: nov_R30215, a 1,728-bp insertion at 3,725,988.
  There are no clipped reads within 3 kb, but depth doubles 500-2,000 bp to
  one side (169x against 84x), because the copied unit is a tandem repeat
  whose junction reads align perfectly to the next unit.
- **Implication:** about a third of the misses are tandem-repeat copy-number
  changes. Only a read-depth method can see those, not a breakpoint caller.
  Another 40% are in repeats where MAPQ cannot place reads.

## Not yet done

- Local assembly of the prototype's candidates, to type and size insertions
  and BND calls.
- A read-depth copy-number scan in the matched reference's frame, for
  tandem-repeat expansions and contractions.
- Phase B with real reads. Simulated reads are cleaner than real ones, and
  dysgu did worse on real ctpV reads than here.
