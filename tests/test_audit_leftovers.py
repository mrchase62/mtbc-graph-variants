#!/usr/bin/env python3
"""Regression tests for the items left open after the section-A audit fixes
(analysis/audit/DECISIONS.md, "Found during the fixes, still open").

Same pattern as tests/run_tests.py: unittest, small synthetic inputs, the
standard library plus the pipeline's own scripts, seconds to run.

    $MTB_PY tests/test_audit_leftovers.py -v
"""
import csv
import importlib.util
import io
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
import warnings

warnings.simplefilter("ignore", ResourceWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
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


def rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


H = "GCF_000195955#1#NC_000962.3"
VCFHDR = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t"
          "INFO\tFORMAT\tS\n")


# --------------------------------------------------------------------------
class OdgiInvertedStepP4(unittest.TestCase):
    """At an inverted step (odgi's own flag, column 5 after frame_convert,
    `-`) the reported target is one base high and R reads complemented.
    p4_place used column 4 alone and placed the record one base off with the
    uncomplemented allele."""

    def test_parse_pos_file_corrects_inverted_steps(self):
        p4 = load("p4_place_lo", "bin/p4_place.py")
        with tempfile.TemporaryDirectory() as d:
            f = write(f"{d}/h.pos", f"""
                #source.path.pos\ttarget.path.pos\tdist.to.ref\tstrand.vs.ref
                R#1#c,99,+\t{H},499,+\t0\t+\t+
                R#1#c,199,+\t{H},599,+\t0\t+\t-
                R#1#c,299,+\t{H},699,+\t0\t-\t-
                R#1#c,399,+\t{H},799,+\t0\t-\t+
                """)
            got = p4.parse_pos_file(f)
        self.assertEqual(got[100], (500, 0, "+"))
        self.assertEqual(got[200], (599, 0, "-"))     # was (600, 0, "+")
        self.assertEqual(got[300], (699, 0, "+"))     # was (700, 0, "-")
        self.assertEqual(got[400], (800, 0, "-"))

    def test_snp_at_inverted_step_placed_at_homolog_complemented(self):
        with tempfile.TemporaryDirectory() as d:
            write(f"{d}/direct.vcf", VCFHDR)
            # R position 200 is A>G; odgi says H37Rv 600 (+,-): the homolog
            # is 599 and H37Rv reads T there
            write(f"{d}/matched.vcf", VCFHDR + "c\t200\t.\tA\tG\t50\t.\t.\tGT\t1\n")
            write(f"{d}/hpos.tsv", f"#src\ttgt\tdist\n"
                  f"R#1#c,199,+\t{H},599,+\t0\t+\t-\n")
            write(f"{d}/npos.tsv", "#src\tnode\nR#1#c,199,+\t77,3,+\n")
            write(f"{d}/mask.bed", "c\t500\t700\tpe_ppe|PPE1|named\n")
            write(f"{d}/loci.tsv", "pos\tlocus_id\n")
            write(f"{d}/graph.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\t"
                  "REF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tOTHER\n")
            write(f"{d}/h37.fa", ">c\n" + "C" * 598 + "T" + "C" * 401 + "\n")
            r = subprocess.run(
                [sys.executable, "bin/p4_place.py", "--sample", "S",
                 "--reference", "RREF", "--build-id", "b",
                 "--direct", f"{d}/direct.vcf", "--matched", f"{d}/matched.vcf",
                 "--h37rv-pos", f"{d}/hpos.tsv", "--node-pos", f"{d}/npos.tsv",
                 "--mask", f"{d}/mask.bed", "--loci", f"{d}/loci.tsv",
                 "--graph-vcf", f"{d}/graph.vcf", "--h37rv", f"{d}/h37.fa",
                 "--out", f"{d}/placed.tsv"], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            got = [(q["h37rv_pos"], q["ref"], q["alt"], q["frame_strand"])
                   for q in rows(f"{d}/placed.tsv")]
        self.assertEqual(got, [("599", "T", "C", "-")])

    def test_p4b_breakpoints_at_inverted_steps(self):
        # R deletes 1001-1500; both breakpoints sit on inverted steps, odgi
        # reporting 3001 for R 1000 and 2501 for R 1500 (one high). The
        # homologs are 3000 and 2500, so H37Rv loses 2500-2999 (pad 2499).
        with tempfile.TemporaryDirectory() as d:
            write(f"{d}/delly.vcf", VCFHDR +
                  "c\t1000\ta\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=1500\tGT\t1/1\n")
            write(f"{d}/bp.pos", f"#src\ttgt\tdist\n"
                  f"R#1#c,999,+\t{H},3000,+\t0\t+\t-\n"
                  f"R#1#c,1499,+\t{H},2500,+\t0\t+\t-\n")
            write(f"{d}/mask.bed", "c\t1\t2\tx\n")
            r = subprocess.run(
                [sys.executable, "bin/p4b_place_sv.py", "--sample", "S",
                 "--reference", "RREF", "--build-id", "b",
                 "--delly", f"{d}/delly.vcf", "--positions", f"{d}/bp.pos",
                 "--mask", f"{d}/mask.bed", "--no-inherited",
                 "--out", f"{d}/sv.tsv"], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            got = [(q["h37rv_pos"], q["h37rv_end"]) for q in rows(f"{d}/sv.tsv")]
        self.assertEqual(got, [("2499", "2999")])        # was 2500, 3000


class OdgiInvertedStepIS6110(unittest.TestCase):
    """is6110_project_sites.py read raw odgi output and ignored odgi's own
    flag (column 4 there); its ISMapper directory defaulted to the pilot's
    and wrote 0 for every isolate without a table there."""

    def run_sites(self, d, extra=()):
        write(f"{d}/frames.tsv",
              "accession\tpanel_len\trefs_len\tstrand\toffset\tagree\n"
              "R1\t100000\t100000\t+\t0\tok\n"
              "GCF_000195955\t100000\t100000\t+\t0\tok\n")
        odgi = write(f"{d}/odgi", textwrap.dedent(f"""\
            #!{sys.executable}
            import sys
            a = sys.argv[1:]
            if a[0] == "paths" and "-L" in a:
                print("R1#1#c"); print("{H}")
            elif a[0] == "paths" and "-H" in a:
                print("path.name\\tpath.length\\tpath.step.count\\tnode.1\\tnode.2\\tnode.3\\tnode.4\\tnode.5")
                print("R1#1#c\\t10\\t5\\t0\\t0\\t0\\t0\\t1")
            elif a[0] == "position":
                q = [l.strip() for l in open(a[a.index("-F") + 1]) if l.strip()]
                for l in q:
                    p, pos, _ = l.rsplit(",", 2)
                    if "-v" in a:
                        print(l + "\\t5,3,+")
                    else:
                        print(l + "\\t{H}," + str(int(pos) + 1000) + ",+\\t0\\t-")
            """))
        os.chmod(odgi, 0o755)
        write(f"{d}/recon.tsv",
              "sample\tgeometry\tchrom_side\tclean_pos\torig_pos\treads\treads_q\n"
              "A\ttsd\tagreed\t5000\t5000\t9\t9\n")
        write(f"{d}/refmap.tsv", "sample\treference\nA\tR1\n")
        env = dict(os.environ, MTB_GRAPH_FRAMES=f"{d}/frames.tsv")
        r = subprocess.run(
            [sys.executable, "is6110/bin/is6110_project_sites.py",
             "--reconcile", f"{d}/recon.tsv", "--refmap", f"{d}/refmap.tsv",
             "--graph", "g.og", "--odgi", odgi, "--workdir", d,
             "--out", f"{d}/out.tsv", "--threads", "1", *extra],
            capture_output=True, text=True, env=env)
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        return rows(f"{d}/out.tsv")[0]

    def test_inverted_step_homolog(self):
        with tempfile.TemporaryDirectory() as d:
            r = self.run_sites(d)
        # query R 0-based 4999 -> odgi H37Rv 0-based 5999 with flag `-`:
        # homolog 0-based 5998, 1-based 5999; was 6000
        self.assertEqual((r["h37rv_pos"], r["frame_strand"]), ("5999", "-"))
        self.assertEqual(r["key"], "h37rv:5999")

    def test_ismapper_join_off_by_default_is_blank_not_zero(self):
        with tempfile.TemporaryDirectory() as d:
            r = self.run_sites(d)
            self.assertEqual(r["ismapper"], "")
            # a directory with no table for this isolate: blank too
            os.makedirs(f"{d}/ism")
            r = self.run_sites(d, ["--ismapper-dir", f"{d}/ism"])
            self.assertEqual(r["ismapper"], "")


# --------------------------------------------------------------------------
class RetierOverlap(unittest.TestCase):
    """retier_intervals.overlapper: only the 8 spans before the interval were
    looked at, and 1-based intervals were compared as if 0-based."""

    def setUp(self):
        self.rt = load("retier_intervals_lo", "bin/retier_intervals.py")

    def test_long_span_more_than_eight_back(self):
        spans = [(0, 1000, "long")] + [(100 + 10 * i, 105 + 10 * i, f"s{i}")
                                       for i in range(10)]
        hit = self.rt.overlapper(sorted(spans))
        self.assertEqual(hit(900, 950), "long")          # was ""

    def test_bed_half_open_against_1_based(self):
        hit = self.rt.overlapper([(99, 100, "a"), (120, 130, "b")])
        # BED [99,100) is 1-based base 100; [120,130) is 121..130
        self.assertEqual(hit(100, 110), "a")             # was ""
        self.assertEqual(hit(101, 120), "")
        self.assertEqual(hit(121, 121), "b")
        self.assertEqual(hit(90, 99), "")


# --------------------------------------------------------------------------
class MergedHeaderMixed(unittest.TestCase):
    """EVIDENCE and SITECLASS can be `mixed`; the header must say so."""

    def test_header_lists_mixed(self):
        m = load("merge_cohort_vcf_lo", "bin/merge_cohort_vcf.py")
        buf = io.StringIO()
        m._header_body(buf.write)
        lines = {l.split("ID=")[1].split(",")[0]: l
                 for l in buf.getvalue().splitlines() if l.startswith("##INFO")}
        self.assertIn("mixed", lines["EVIDENCE"])
        self.assertIn("mixed", lines["SITECLASS"])


# --------------------------------------------------------------------------
class SanityMissingH37Rv(unittest.TestCase):
    """p5_sanity check 3 wrote 0, and took it into the median, for a sample
    with no H37Rv count."""

    def test_na_not_zero(self):
        with tempfile.TemporaryDirectory() as d:
            hdr = ["key", "frame", "region", "kind", "h37rv_pos", "n_alt",
                   "n_nocall", "A", "B", "C"]
            with open(f"{d}/m.tsv", "w") as fh:
                fh.write("\t".join(hdr) + "\n")
                for i in range(20):
                    st = ["ALT", "ALT", "ALT"] if i < 10 else \
                         ["REF", "ALT", "REF"]
                    fh.write("\t".join([f"h37rv:{i+1}", "h37rv", "core", "SNP",
                                        str(i + 1), "1", "0"] + st) + "\n")
            write(f"{d}/refmap.tsv", "sample\tsnp_distance\nA\t1\nB\t2\nC\t3\n")
            write(f"{d}/cohort.tsv", "sample\tmeandepth\nA\t50\nB\t60\nC\t70\n")
            # A: 10 ALT / 10 = 1.0; B: 20 / 10 = 2.0; C: no count
            write(f"{d}/p2.tsv", "sample\th37rv_small\nA\t10\nB\t10\n")
            r = subprocess.run(
                [sys.executable, "bin/p5_sanity.py", "--matrix", f"{d}/m.tsv",
                 "--refmap", f"{d}/refmap.tsv", "--p2-summary", f"{d}/p2.tsv",
                 "--cohort", f"{d}/cohort.tsv", "--graph-vcf", f"{d}/none.vcf.gz",
                 "--out", f"{d}/s.tsv"], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            got = {q["check"]: q for q in rows(f"{d}/s.tsv")}
        self.assertEqual(got["alt_vs_p2"]["value"], "1.5")   # was 1.0
        self.assertIn("1 samples NA", got["alt_vs_p2"]["detail"])
        self.assertIn("NA", r.stdout)


# --------------------------------------------------------------------------
_assoc = None


def assoc():
    global _assoc
    if _assoc is None:
        _assoc = load("test_audit_assoc_lo", "tests/test_audit_assoc.py")
    return _assoc


class SmallDeletionSpan(unittest.TestCase):
    """A small deletion was credited to the unit at its anchor base, not the
    genes whose bases it removes (the SV rule of ASSOC-5)."""

    def setUp(self):
        a = assoc()
        self.base = a.BurdenBase("units")
        self.base.setUp()

    def test_deletion_crossing_gene_boundary_credits_both(self):
        # G1 1001-1100, G2 1101-1200; anchor 1099 deletes 1100-1102
        u = self.base.units([("1099|AGGG|A", "h37rv:1099", 1099, "small",
                              "core", ["L01", "L05"], 0)], "small")
        self.assertIn("G1", u)
        self.assertIn("G2", u)                           # was G1 only

    def test_insertion_and_snp_stay_at_anchor(self):
        u = self.base.units([("1100|A|AT", "h37rv:1100", 1100, "small",
                              "core", ["L01", "L05"], 0),
                             ("1100|A|G", "h37rv:1100s", 1100, "small",
                              "core", ["L02", "L06"], 0)], "small")
        self.assertEqual(set(u), {"G1"})

    def test_span_helper(self):
        bur = self.base.bur
        self.assertEqual(bur.small_deleted_span(
            {"pos": "1099", "ref": "AGGG", "alt": "A,*"}), (1100, 1102))
        self.assertIsNone(bur.small_deleted_span(
            {"pos": "1099", "ref": "AG", "alt": "TC"}))
        self.assertIsNone(bur.small_deleted_span(
            {"pos": "1099", "ref": "A", "alt": "AT"}))


if __name__ == "__main__":
    unittest.main()
