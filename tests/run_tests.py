#!/usr/bin/env python3
"""Regression tests for defects found in the 2026-09-29 review (CODE_REVIEW.md).

Every test here pins one bug that produced plausible-looking wrong output
rather than an error. They run on small synthetic inputs, need nothing beyond
the standard library and the pipeline's own scripts, and take a few seconds:

    $MTB_PY tests/run_tests.py            # from the repository root
    $MTB_PY tests/run_tests.py -v

Run them before any cohort run that follows a code change.
"""
import csv
import importlib.util
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
import warnings

# the scripts under test open files without context managers; that is not
# what these tests are about
warnings.simplefilter("ignore", ResourceWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)                      # several scripts load siblings by relative path
sys.path.insert(0, os.path.join(ROOT, "bin"))


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, rel))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def write(path, text):
    with open(path, "w") as fh:
        fh.write(textwrap.dedent(text).lstrip("\n"))
    return path


class RepeatMaskLookup(unittest.TestCase):
    """CODE_REVIEW 5.1: nested intervals missed, and a 0/1-based slip."""

    def setUp(self):
        self.p4 = load("p4_place", "bin/p4_place.py")

    def test_nested_interval_outer_tail_is_masked(self):
        # PE_PGRS4 336359-339273 contains 336559-339142 (BED, 0-based half-open)
        iv = self.p4.merge_intervals([(336359, 339273), (336559, 339142)])
        hit = self.p4.make_hit(iv)
        self.assertTrue(hit(339200))          # past the inner end: was "core"

    def test_bed_edges_against_1_based_positions(self):
        hit = self.p4.make_hit(self.p4.merge_intervals([(100, 200)]))
        # BED [100, 200) covers 1-based 101..200
        self.assertFalse(hit(100))
        self.assertTrue(hit(101))
        self.assertTrue(hit(200))
        self.assertFalse(hit(201))

    def test_class_prefix_not_substring(self):
        with tempfile.TemporaryDirectory() as d:
            bed = write(os.path.join(d, "m.bed"), """
                c\t10\t20\tpe_ppe|PPE1|named
                c\t30\t40\ttandem|pepA|tandem100
                """)
            pe, other = self.p4.load_intervals(bed)
        self.assertEqual(pe, [(10, 20)])
        self.assertEqual(other, [(30, 40)])   # "pepA" contains PE; not PE/PPE


class GeneLookup(unittest.TestCase):
    """CODE_REVIEW 3.8: a gene containing a shorter nested gene was missed."""

    def test_nested_gene(self):
        p6 = load("p6_annotate", "bin/p6_annotate.py")
        rows = [(100, 1000, "outer", "", ""), (200, 300, "inner", "", "")]
        hit = p6.make_lookup(rows)
        self.assertEqual(hit(250)[2], "inner")
        self.assertEqual(hit(500)[2], "outer")   # was None: scan stopped at inner
        self.assertIsNone(hit(1001))


class Is6110Stage2Keys(unittest.TestCase):
    """CODE_REVIEW 4.1: node keys lost their offset and never found a carrier."""

    def test_node_key_keeps_offset(self):
        st2 = load("is6110_p5_stage2", "is6110/bin/is6110_p5_stage2.py")
        self.assertEqual(st2.short("node:46502:0:N><INS>"), "node:46502:0:")
        self.assertEqual(st2.short("h37rv:932204:C><INS>"), "h37rv:932204:")
        # the carrier table builds f'node:{node}:' from node = "<id>:<off>"
        self.assertEqual(st2.short("node:46502:0:N><INS>"), f'node:{"46502:0"}:')


class GvcfEvidence(unittest.TestCase):
    """CODE_REVIEW 3.3: any gVCF line counted as REF evidence."""

    def test_variant_lines_are_not_reference_coverage(self):
        p5 = load("p5_states", "bin/p5_states.py")
        with tempfile.TemporaryDirectory() as d:
            g = write(os.path.join(d, "s.g.vcf"), """
                ##fileformat=VCFv4.2
                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts
                c\t1\t.\tA\t<NON_REF>\t.\t.\tEND=99\tGT:DP\t0:30
                c\t100\t.\tC\tA,<NON_REF>\t.\t.\t.\tGT:DP\t1:30
                c\t101\t.\tGTT\tG,<NON_REF>\t.\t.\t.\tGT:DP\t1:30
                c\t102\t.\tT\t<NON_REF>\t.\t.\tEND=200\tGT:DP\t0:30
                """)
            blocks, spans, alt_at = p5.load_gvcf(g)
        cov = p5.make_cov(blocks, 5)
        inside = p5.make_in_spans(spans)
        self.assertTrue(cov(50))
        self.assertFalse(cov(100))            # a SNP call is not a reference block
        self.assertTrue(inside(100))
        self.assertTrue(inside(103))          # deleted base, though a block covers it
        self.assertFalse(inside(104))
        self.assertEqual(alt_at[100][:2], ("C", "A"))

    def test_spans_with_long_early_span(self):
        p5 = load("p5_states", "bin/p5_states.py")
        inside = p5.make_in_spans([(1, 100), (5, 6), (50, 55), (200, 201)])
        self.assertTrue(inside(99))           # reached only through the first span
        self.assertFalse(inside(150))


class SvReverseStrand(unittest.TestCase):
    """CODE_REVIEW 5.2: reverse-strand SVs came out start > end, shifted."""

    def run_p4b(self, d, proj_lines, svtype="DEL", pos=1000, end=3000):
        vcf = write(os.path.join(d, "s.delly.vcf"), f"""
            ##fileformat=VCFv4.2
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts
            R\t{pos}\t.\tN\t<{svtype}>\t60\tPASS\tSVTYPE={svtype};END={end};SVLEN={end - pos}\tGT\t1
            """)
        posf = write(os.path.join(d, "bp.pos"),
                     "#source.path.pos\ttarget.path.pos\tdist.to.ref\tstrand.vs.ref\n"
                     + "".join(proj_lines))
        mask = write(os.path.join(d, "m.bed"), "c\t1\t2\ttandem|x|y\n")
        out = os.path.join(d, "out.tsv")
        subprocess.run([sys.executable, "bin/p4b_place_sv.py", "--sample", "s",
                        "--reference", "R", "--build-id", "b", "--delly", vcf,
                        "--positions", posf, "--mask", mask, "--no-inherited",
                        "--out", out], check=True, capture_output=True)
        with open(out) as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    def test_reverse_deletion(self):
        # R bases 1001..3000 deleted; R runs reverse to H37Rv, R x <-> H 5001-x
        # (0-based in the file: R 999 -> H 4001, R 2999 -> H 2001)
        with tempfile.TemporaryDirectory() as d:
            rows = self.run_p4b(d, [
                "R#1#c,999,+\tH#1#h,4001,+\t0\t-\t+\n",
                "R#1#c,2999,+\tH#1#h,2001,+\t0\t-\t+\n"])
        r = rows[0]
        self.assertEqual(r["frame"], "h37rv")
        # deleted in H: f(3000)=2002 .. f(1001)=4001; padding base 2001
        self.assertEqual((int(r["h37rv_pos"]), int(r["h37rv_end"])), (2001, 4001))
        self.assertEqual(r["key"], "sv:DEL:2001:4001")

    def test_forward_unchanged(self):
        with tempfile.TemporaryDirectory() as d:
            rows = self.run_p4b(d, [
                "R#1#c,999,+\tH#1#h,1099,+\t0\t+\t+\n",
                "R#1#c,2999,+\tH#1#h,3099,+\t0\t+\t+\n"])
        self.assertEqual(rows[0]["key"], "sv:DEL:1100:3100")

    def test_span_mismatch_is_breakend(self):
        with tempfile.TemporaryDirectory() as d:
            rows = self.run_p4b(d, [
                "R#1#c,999,+\tH#1#h,1099,+\t0\t+\t+\n",
                "R#1#c,2999,+\tH#1#h,90099,+\t0\t+\t+\n"])
        self.assertEqual(rows[0]["frame"], "bnd")

    def test_no_calls_still_writes_header(self):
        with tempfile.TemporaryDirectory() as d:
            empty = write(os.path.join(d, "e.vcf"), "##fileformat=VCFv4.2\n")
            mask = write(os.path.join(d, "m.bed"), "c\t1\t2\ttandem|x|y\n")
            posf = write(os.path.join(d, "bp.pos"), "")
            out = os.path.join(d, "out.tsv")
            subprocess.run([sys.executable, "bin/p4b_place_sv.py", "--sample", "s",
                            "--reference", "R", "--build-id", "b", "--delly", empty,
                            "--positions", posf, "--mask", mask, "--no-inherited",
                            "--out", out], check=True, capture_output=True)
            with open(out) as fh:
                self.assertTrue(fh.readline().startswith("sample\treference"))

    def test_missing_graph_vcf_is_fatal(self):
        with tempfile.TemporaryDirectory() as d:
            mask = write(os.path.join(d, "m.bed"), "c\t1\t2\ttandem|x|y\n")
            posf = write(os.path.join(d, "bp.pos"), "")
            r = subprocess.run([sys.executable, "bin/p4b_place_sv.py", "--sample", "s",
                                "--reference", "R", "--build-id", "b",
                                "--positions", posf, "--mask", mask,
                                "--graph-vcf", os.path.join(d, "nope.vcf.gz"),
                                "--out", os.path.join(d, "o.tsv")],
                               capture_output=True, text=True)
            self.assertNotEqual(r.returncode, 0)


class RunnerCohortTable(unittest.TestCase):
    """CODE_REVIEW 6.7 and 6.3: blank lines inflated N; --cohort was rejected."""

    def runner(self, table_text, *extra):
        d = tempfile.mkdtemp()
        tab = write(os.path.join(d, "cohort.tsv"), table_text)
        crams = write(os.path.join(d, "crams.tsv"), "sample\tcram_relpath\n")
        reg = write(os.path.join(d, "reg.tsv"),
                    f"t\t{tab}\t{crams}\t{d}/out\tp1\tnote\t{d}/work/\n")
        # a minimal completed build, so the runner gets past build selection
        b = os.path.join(d, "build", "b0")
        os.makedirs(os.path.join(b, "logs"))
        og = write(os.path.join(d, "g.og"), "x\n")
        write(os.path.join(b, "logs", "manifest.done"), "done\n")
        write(os.path.join(b, "build_info.tsv"), f"build_id\tb0\ngraph\t{og}\n")
        frames = write(os.path.join(d, "frames.tsv"), "accession\n")
        env = dict(os.environ, REGISTRY=reg, BUILD_ROOT=os.path.join(d, "build"),
                   MTB_GRAPH_FRAMES=frames)
        return subprocess.run(["bash", "bin/refbias_run.sh", "--cohort", "t",
                               "--dry-run", *extra], capture_output=True,
                              text=True, env=env)

    def test_cohort_flag_accepted(self):
        r = self.runner("sample\tlineage\nA\t4\nB\t4\n")
        self.assertNotIn("unknown option", r.stderr)
        self.assertEqual(r.returncode, 0, r.stderr)

    def test_trailing_blank_line_not_counted(self):
        r = self.runner("sample\tlineage\nA\t4\nB\t4\n\n")
        self.assertRegex(r.stdout + r.stderr, r"cohort\s*:\s*t\s+2 isolates")

    def test_blank_line_between_rows_refused(self):
        r = self.runner("sample\tlineage\nA\t4\n\nB\t4\n")
        self.assertIn("between sample rows", r.stderr)


class ConfigSiteFile(unittest.TestCase):
    """CODE_REVIEW 2.1 and 2.2: site file order and working directory."""

    def test_site_file_moves_derived_paths_from_any_directory(self):
        d = tempfile.mkdtemp()
        cfg = os.path.join(d, "config")
        os.makedirs(cfg)
        with open(os.path.join(ROOT, "config", "project_env.sh")) as src, \
                open(os.path.join(cfg, "project_env.sh"), "w") as dst:
            dst.write(src.read())
        write(os.path.join(cfg, "site.local.sh"), "export MTB_WORK=/elsewhere\n")
        r = subprocess.run(
            ["env", "-i", f"HOME={os.environ.get('HOME', '/tmp')}", "bash", "-c",
             f"cd / && source {cfg}/project_env.sh && echo $MTB_DATA"],
            capture_output=True, text=True)
        self.assertEqual(r.stdout.strip(), "/elsewhere/data")


if __name__ == "__main__":
    unittest.main(warnings="ignore")
