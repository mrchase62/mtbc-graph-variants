# Fix brief: section A of the audit (2026-10-05)

You are fixing a group of audit findings in the mtbc-graph-variants
repository. You work in your own git worktree (the Agent tool created it); your
branch will be merged by the coordinating session.

**Read first:**

- `analysis/audit/CONSOLIDATED.md` and the per-area report(s) named in your
  task, in the main repository at
  `/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants/analysis/audit/`.
  Your worktree has no `analysis/`, so read them there.
- `CODE_REVIEW.md` and `HANDOFF.md` section 0m, in your worktree.

## Rules (strict)

1. **Change only the files your task assigns to you.** If a fix needs a change
   elsewhere, describe it in your report instead.
2. **No cluster jobs** (no sbatch or srun). Local runs of a few minutes are
   fine.
3. **Never write in the working tree,** `/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome`
   (it is read-only test data), **nor in the main repository checkout.**
   - Write only in your worktree and in `/tmp/claude-12043/fix_<group>/`.
   - Real inputs can be read through
     `/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants/runroot/`,
     which links the working tree's data. Run YOUR worktree's scripts by path
     against those inputs, with outputs in your scratch folder.
4. **Fix the finding and nothing else.** No refactors, no new methods, no
   changes in behaviour beyond what the finding requires. Keep the
   surrounding code's style and comment density.
5. **Where the correct behaviour needs a scientific choice,** implement the
   conservative option and mark it `DECISION` in your report with the
   alternatives. Conservative means unknown rather than a guess: write
   missing, not REF or ALT.
6. **Do not touch the canettii-isolate question (ASSOC-1).** The user deferred
   it.
7. Commit your work on your worktree branch with clear messages, ending with
   `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Do not push.

## Done means

- **A regression test for each fixed finding,** in a NEW file
  `tests/test_audit_<group>.py`, following `tests/run_tests.py`'s pattern:
  - unittest;
  - small synthetic inputs;
  - standard library plus the pipeline's own scripts;
  - seconds to run.

  The test must fail on the old code and pass on the new. State that you
  checked this.
- **The existing suite still passes:** `$MTB_PY tests/run_tests.py`, after
  `source config/project_env.sh`.
- **Validation on real data:**
  - run the fixed code on copies of scale200 inputs (and gwas1000 where
    cheap);
  - report before and after counts against the audit's numbers, so the
    coordinating session can predict what the rerun should show.

## Report (your final message)

- the branch name and the commit list;
- per finding: fixed / partly / not fixed, the test name, and the before and
  after counts;
- every `DECISION` and its alternatives;
- anything found along the way that is not in the audit.
