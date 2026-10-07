# Step 4: independent second review of the audit fixes (2026-10-06)

You are reviewing the fixes on branch `audit-fixes` (ffd7bb9) against `main`.
The diff is about 95 files and 11,500 lines. It was written by eleven separate
fix groups, each working in isolation, and merged afterwards. Your job is to
catch errors in the fixes themselves before they reach main and before a new
graph is built.

**The user's priority:** right answers over speed. Downstream analysis is
almost all tree-based, and ancestral reconstruction is crucial.

## Where things are

- **Repository:** `/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants`.
  - To read the branch: `git -C <repo> show audit-fixes:<path>`.
  - To read the diff: `git -C <repo> diff main...audit-fixes -- <paths>`.
  - A ready checkout is at `<repo>/.claude/worktrees/integration`. It is
    **read-only for you**; do not commit or edit there.
- **Audit and fix records** (under `<repo>/analysis/audit/`):
  - `CONSOLIDATED.md`: the original findings;
  - the per-area reports;
  - `FIX_BRIEF.md`: what the fix groups were told;
  - `DECISIONS.md`: 44 choices said to be implemented.
- **Production data, read-only:** `/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome`,
  also reachable through `<repo>/runroot/`.

## Rules (strict)

1. **Read-only for code and data.** Write only:
   - your report, `analysis/audit/review2/<area>.md` in the main checkout;
   - scratch files in `/tmp/claude-12043/review2_<area>/`.
2. **No cluster jobs** (no sbatch or srun). Local runs of a few minutes are
   fine. Run the branch's code from the integration checkout or from a
   `git archive audit-fixes` export in your scratch folder.
3. **Never write into the working tree.**
4. **Do not trust the fix groups' reports; verify.** Reproduce the key claims
   you rely on.

## What to check

1. **Correctness of each fix:**
   - Does it fix the finding fully, in every code path, including the sharded
     and assembled paths, per-sample and cohort modes, and both frames?
   - Is it what `DECISIONS.md` says was implemented?
2. **New bugs introduced:**
   - off-by-one and frame errors (0- vs 1-based, H37Rv vs reference vs node
     frame, strand);
   - allele orientation;
   - key formats that changed in one producer but not in every reader;
   - silently dropped records;
   - new required arguments not passed by callers (the chain, `refbias_run.sh`
     and the `.sh` wrappers);
   - performance traps at 10,000 samples.
3. **Interactions between groups.** Several groups edited the same files
   (`merge_cohort_vcf.py`, `p5_states.py`, `p4_place.py`,
   `cohort_assoc_tail.sh`, `p0_prepare.sh`, `add_outgroup.py`,
   `write_event_matrix.py`). Check that the merged result is coherent:
   - one group's assumption is not broken by another's change;
   - node-offset normalisation is applied in every producer and reader of
     node keys;
   - the build-ID guards do not reject legitimate fresh runs.
4. **Tests:**
   - Do they test the behaviour, or only that a function runs?
   - Would they catch a regression?
   - Are any skipped silently where tools are absent?
5. **Coverage against the audit.** Does every HIGH and MEDIUM finding in
   `CONSOLIDATED.md` and the per-area reports have a fix, a documented
   decision, or a stated reason it is out of scope? List any that slipped
   through.
6. **End-to-end coherence.** Trace a variant, an IS6110 insertion, an SV
   deletion and an accessory locus from the placement output to the merged
   VCF to the event matrix to the scan. Do the formats and semantics agree at
   every hand-off? Where cheap, do it on real data (a few scale200 samples).

## Report format

**One section per finding:**

- an ID (`R2-<AREA>-<n>`) and a severity:
  - **HIGH:** wrong results;
  - **MEDIUM:** wrong in a minority of cases, or would be at scale;
  - **LOW:** robustness or clarity;
- `file:line` on `audit-fixes`;
- what is wrong;
- the evidence (code excerpt, and a reproduction or count wherever possible);
- the effect;
- a suggested fix.

**Then:**

- "checked and found sound", with what was checked, so coverage is visible;
- "audit findings not addressed";
- open questions.

Mark anything you could not verify as UNVERIFIED. False alarms cost the user
time; missed bugs cost a rerun.
