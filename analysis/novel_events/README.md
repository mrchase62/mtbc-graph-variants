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

## Steps 1 and 2 (2026-10-08; job 51415882, 41 min)

Phase A's total is now **30.0 billing-hours**, the approved budget.

**Step 1: local assembly of each prototype v2 candidate**
(`assemble_candidates.py` -> `out/proto_v2_asm`).
- Candidates within 500 bp share a window.
- Reads within 1.5 kb of the window, plus their mates, come from the P2 BAM
  and are assembled with SPAdes `--isolate`.
- Contigs of 300 bp or more are aligned to the matched reference with
  minimap2 `asm5`.
- Events are read off the contigs with the truth-set rules.
- Windows that give no events keep the prototype's calls.

**Step 2: read-depth scan** (`depth_scan.py` -> `out/depth`).
- samtools depth with all mapping qualities, in 100 bp windows, against the
  sample's median window depth.
- 2 or more windows at 1.4x or above is a gain (INS); 2 or more at 0.6x or
  below is a loss (DEL).
- A true event matches when it lies within the run +-200 bp. Long runs do not
  inflate this: 20 runs over 5 kb hold 14 events, and recall without them is
  0.373 against 0.386.

| caller | typed recall | breakpoint found | precision |
|---|---:|---:|---:|
| prototype v2 | 0.33 | 0.56 | 0.62 |
| prototype v2 + assembly | 0.38 | 0.58 | 0.67 |
| depth scan | 0.36 | 0.39 | 0.44 |
| prototype + assembly + depth | 0.61 | 0.76 | 0.54 |
| dysgu PASS | 0.51 | 0.53 | 0.65 |
| prototype + assembly + dysgu PASS | 0.60 | 0.66 | 0.66 |
| **prototype + assembly + depth + dysgu PASS** | **0.77** | **0.81** | **0.58** |

**By size** (all four together; in brackets, dysgu PASS alone):

| size | typed recall |
|---|---|
| 50-150 bp | 0.64 (0.48) |
| 150-500 bp | 0.80 (0.45) |
| 500-2,000 bp | 0.91 (0.60) |
| > 2 kb | 0.91 (0.55) |

**What each step adds:**
- **The depth scan** carries the 150-500 bp and > 2 kb classes, which are the
  tandem copy-number changes: recall 0.64 and 0.68 alone. It is weak below
  150 bp (0.13).
- **Assembly** raises precision (0.62 -> 0.67) and types more 50-150 bp events
  (0.24 -> 0.36). It does not yet type IS6110-sized insertions: their
  contigs end inside the element, so the breakpoint is found (0.93) but the
  event is not sized.
  - dysgu types those (0.90), so the union is 0.92.

**Precision:**
- The depth scan's 0.44 is the main false-call source: 488 of 866 runs match
  no true event. Some may be events the truth set misses, but most are
  probably depth noise.
- Real reads add GC bias, which makes depth calls noisier.

**Next:** Phase B, scored with real reads. The 63 Marin Illumina runs are
downloading. Running P1/P2 on them needs the user's approval: about 30
billing-hours at the novelA40 rate (P1 + P2 about 0.6 per sample), plus about
6 for the callers.
