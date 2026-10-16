"""Checks against lost run time (2026-10-15): the P5 inline-projection guard,
the watchdog, and refbias_run's preflight and watchdog submission."""
import os, shutil, subprocess, unittest

R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class ChainChecks(unittest.TestCase):
    def test_scripts_parse(self):
        for s in ("bin/chain_watch.sh", "bin/p5_merge.sh", "bin/refbias_run.sh"):
            r = subprocess.run(["bash", "-n", os.path.join(R, s)], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, s + r.stderr)

    def test_p5_states_refuses_large_inline_projection(self):
        src = open(os.path.join(R, "bin/p5_merge.sh")).read()
        self.assertIn('if [[ "$STEP" == --states && -s "$MISS" ]]', src)
        self.assertIn("P5_INLINE_PROJ_MAX", src)
        # the guard sits before the inline odgi call
        self.assertLess(src.index("P5_INLINE_PROJ_MAX"), src.index('"$ODGI" position -i "$OG" -F "$MISS"'))

    def test_project_step_precedes_states_in_chain(self):
        src = open(os.path.join(R, "bin/refbias_run.sh")).read()
        self.assertIn('P5STEP=--project', src)
        self.assertIn('"${JOB[p5proj]:-}" "-c 2 --mem=4G"', src)

    def test_watchdog_is_default_and_optional(self):
        src = open(os.path.join(R, "bin/refbias_run.sh")).read()
        self.assertIn("WATCH=1", src)
        self.assertIn("--no-watch) WATCH=0", src)
        self.assertIn("bin/chain_watch.sh", src)

    @unittest.skipUnless(shutil.which("squeue") and shutil.which("sacct"), "no Slurm client")
    def test_watchdog_once_on_unknown_cohort_is_quiet(self):
        r = subprocess.run(["bash", os.path.join(R, "bin/chain_watch.sh"), "no_such_cohort_xyz", "--once"],
                           capture_output=True, text=True, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertNotIn("FAILURE", r.stdout)


class P0Projections(unittest.TestCase):
    def test_step_exists_and_is_not_in_the_default_run(self):
        src = open(os.path.join(R, "bin/p0_prepare.sh")).read()
        self.assertIn("step_projections()", src)
        default = src[src.rindex('if [[ -n "$STEP" ]]; then'):]
        self.assertNotIn("_run projections", default)
        self.assertIn("PROJ_EXTRA", src)
        # per-genome completion is tied to the position list it was made for
        self.assertIn("positions_sha", src)


if __name__ == "__main__":
    unittest.main()
