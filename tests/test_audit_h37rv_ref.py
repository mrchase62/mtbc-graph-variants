#!/usr/bin/env python3
"""H37Rv as a candidate matched reference (the user's decision D21,
2026-10-07), and what P2-P5 need for R = H37Rv.

H37Rv was never a candidate (audit P0P2-7): it is the deconstruct reference,
so it has no panel VCF column. It is now an all-REF column. Running the
H37Rv-pinned arm's inputs through today's P4 found that projecting H37Rv
onto itself through the graph is not the identity where its path passes a
node more than once (729 of 51,139 composed records moved, all in PE/PPE
tandem repeats). A path projected onto itself is now the identity, in
frame_convert.py (every shell call site) and in the two IS6110 scripts that
read odgi directly.
"""
import csv
import importlib.util
import os
import subprocess
import sys
import tempfile
import textwrap
import types
import unittest
import warnings

warnings.simplefilter("ignore", ResourceWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))

import test_audit_selection as sel                               # noqa: E402
import test_audit_leftovers as lo                                # noqa: E402

vcf, py, H = sel.vcf, sel.py, lo.H

FRAMES = ("accession\tpanel_len\trefs_len\tstrand\toffset\tagree\n"
          "R1\t100000\t100000\t+\t0\tok\n"
          "GCF_000195955\t100000\t100000\t+\t0\tok\n")


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, rel))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def write(path, text):
    with open(path, "w") as fh:
        fh.write(textwrap.dedent(text).lstrip("\n"))
    return path


# --------------------------------------------------------------------------
class H37RvIsACandidate(unittest.TestCase):
    """t8_select_reference.py adds H37Rv as a column REF at every panel site."""

    H = "GCF_000195955"
    SAMPLES = ["gA", "gG"]
    select = sel.SelectionByAllele.select

    def test_isolate_nearest_h37rv_picks_it(self):
        # gA and gG each carry one ALT the isolate lacks; H37Rv carries none
        got = self.select([(100, "C", "T", "", ["1", "0"]),
                           (200, "C", "T", "", ["0", "1"])], [],
                          full=True, top=3, h37rv=self.H)
        self.assertEqual(got[0]["reference"], self.H)
        self.assertEqual(got[0]["snp_distance"], "0")

    def test_h37rv_distance_is_the_isolates_alt_calls_at_covered_sites(self):
        # the isolate carries 100 and 200; 200 is uncovered (masked by depth)
        d = self.select([(100, "C", "T", "", ["1", "1"]),
                         (200, "C", "T", "", ["1", "1"]),
                         (300, "C", "T", "", ["0", "1"])],
                        [(100, "C", "T", "1"), (200, "C", "T", "1")],
                        depth={100: 30, 200: 0, 300: 30}, top=3, h37rv=self.H)
        self.assertEqual(d, {"gA": 0, "gG": 1, self.H: 1})

    def test_compared_over_every_covered_site(self):
        # gG has a missing cell; H37Rv has none
        got = self.select([(100, "C", "T", "", ["1", "."]),
                           (200, "C", "T", "", ["1", "1"])],
                          [(100, "C", "T", "1"), (200, "C", "T", "1")],
                          full=True, top=3, h37rv=self.H)
        n = {x["reference"]: x["n_compared"] for x in got}
        self.assertEqual(n[self.H], "2")
        self.assertEqual(got[-1]["reference"], self.H)

    def test_a_panel_column_named_h37rv_is_refused(self):
        self.SAMPLES = ["gA", self.H]
        r = self.select([(100, "C", "T", "", ["1", "0"])], [],
                        check=False, h37rv=self.H)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("already has a column", r.stderr)

    def test_on_by_default(self):
        with tempfile.TemporaryDirectory() as d:
            panel = vcf(os.path.join(d, "p.vcf"), ["gA"],
                        [(100, "C", "T", "", ["1"])])
            iso = vcf(os.path.join(d, "i.vcf"), ["S"], [])
            dep = write(os.path.join(d, "dep.tsv"), "c\t100\t30\n")
            out = os.path.join(d, "o.tsv")
            subprocess.run([py(), "bin/t8_select_reference.py", "--vcf", iso,
                            "--panel-snps", panel, "--depth", dep, "--out", out],
                           check=True, capture_output=True, text=True)
            got = list(csv.DictReader(open(out), delimiter="\t"))
        self.assertEqual(got[0]["reference"], self.H)


class H37RvIntervalCount(unittest.TestCase):
    """P1's tie-break needs an IS6110 interval count for every candidate;
    P0 skipped H37Rv, so a tie involving it would stop P1's summary."""

    def test_p0_counts_h37rv(self):
        s = open("bin/p0_prepare.sh").read()
        i = s.index("step_is6110_intervals() {")
        body = s[i:s.index("\n}\n", i)]
        self.assertIn("is6110-intervals", body)
        self.assertNotIn("--skip", body)


# --------------------------------------------------------------------------
class PathOntoItselfIsIdentity(unittest.TestCase):
    """odgi projecting a path onto itself answers with one of the copies
    where the path passes a node more than once."""

    def convert(self, lines):
        with tempfile.TemporaryDirectory() as d:
            write(f"{d}/frames.tsv", FRAMES)
            r = subprocess.run(
                [sys.executable, "graphframe/bin/frame_convert.py", "from-panel",
                 "--table", f"{d}/frames.tsv"],
                input="".join(l + "\n" for l in lines),
                capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        return [l.split("\t") for l in r.stdout.splitlines()], r.stderr

    def test_same_path_is_written_as_the_identity(self):
        out, err = self.convert([f"{H},90754,+\t{H},90757,+\t0\t-"])
        self.assertEqual(out[0][1], f"{H},90754,+")
        self.assertEqual(out[0][2], "0")
        self.assertEqual(out[0][3:], ["+", "+"])
        self.assertIn("1 projected onto their own path", err)

    def test_other_paths_unchanged(self):
        out, _ = self.convert([f"R1#1#c,100,+\t{H},150,+\t0\t-"])
        self.assertEqual(out[0][1], f"{H},150,+")
        self.assertEqual(out[0][3:], ["+", "-"])

    def test_node_form_unchanged(self):
        out, _ = self.convert([f"{H},100,+\t175,0,+"])
        self.assertEqual(out[0], [f"{H},100,+", "175,0,+"])

    def test_p4_reads_the_identity(self):
        p4 = load("p4_place_d21", "bin/p4_place.py")
        out, _ = self.convert([f"{H},90754,+\t{H},90757,+\t0\t-"])
        with tempfile.TemporaryDirectory() as d:
            p = write(f"{d}/h.pos", "#h\n" + "\t".join(out[0]) + "\n")
            self.assertEqual(p4.parse_pos_file(p), {90755: (90755, 0, "+")})


class Is6110SameReferenceIsIdentity(unittest.TestCase):
    """is6110_project_sites.py: an isolate matched to H37Rv projects its
    sites onto H37Rv; the fake odgi answers 1000 bases away on the other
    strand, as it can at a repeat."""

    def test_h37rv_sample_keeps_its_position(self):
        with tempfile.TemporaryDirectory() as d:
            # the leftovers fixture (sample A on R1), then A on H37Rv
            lo.OdgiInvertedStepIS6110.run_sites(self, d)
            write(f"{d}/refmap.tsv", "sample\treference\nA\tGCF_000195955\n")
            env = dict(os.environ, MTB_GRAPH_FRAMES=f"{d}/frames.tsv")
            r = subprocess.run(
                [sys.executable, "is6110/bin/is6110_project_sites.py",
                 "--reconcile", f"{d}/recon.tsv", "--refmap", f"{d}/refmap.tsv",
                 "--graph", "g.og", "--odgi", f"{d}/odgi", "--workdir", d,
                 "--out", f"{d}/out2.tsv", "--threads", "1"],
                capture_output=True, text=True, env=env)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            got = list(csv.DictReader(open(f"{d}/out2.tsv"), delimiter="\t"))[0]
        self.assertEqual((got["h37rv_pos"], got["frame_strand"]), ("5000", "+"))
        self.assertEqual(got["key"], "h37rv:5000")


class Is6110Stage2SameReference(unittest.TestCase):
    """is6110_p5_stage2.py: a key carried by a sample of reference T is
    looked up on T itself at its own position, whatever odgi or the
    projection store answers."""

    def test_store_answer_for_same_reference_is_overridden(self):
        st2 = load("is6110_p5_stage2_d21", "is6110/bin/is6110_p5_stage2.py")
        sent = []
        st2.stage1_rows = lambda a, s: [dict(key="h37rv:500:C><INS>",
                                              state="NOCALL"),
                                         dict(key="h37rv:900:C><INS>",
                                              state="NOCALL")]
        # the store holds odgi's wrong copy for the same-reference carrier
        st2.store_load = lambda a, t: {("T", 500): (569, 0)}
        st2.store_add = lambda a, t, res: None

        def fake_odgi(a, tref, tpath, fr, paths, og, queries):
            sent.extend(queries)
            return {q: (q[1] + 7, 0) for q in queries}
        st2.run_odgi = fake_odgi
        with tempfile.TemporaryDirectory() as d:
            a = types.SimpleNamespace(workdir=d)
            st2.project(a, "T", None, {"T": "T#1#c", "C": "C#1#c"}, "g.og",
                        {"B": "T"}, {"h37rv:500:": ("T", 500),
                                     "h37rv:900:": ("C", 900)})
            got = {r["short_key"]: (r["pos"], r["dist"]) for r in
                   csv.DictReader(open(f"{d}/proj/T.tsv"), delimiter="\t")}
        self.assertEqual(got["h37rv:500:"], ("500", "0"))     # was 569
        self.assertEqual(got["h37rv:900:"], ("907", "0"))     # other ref: odgi
        self.assertEqual(sent, [("C", 900)])


class P4H37RvInheritsNothing(unittest.TestCase):
    """R = H37Rv has no graph VCF column and nothing to inherit; that is
    not the 'inherited half is EMPTY' failure P4 warns of."""

    def test_no_warning_for_h37rv(self):
        s = open("bin/p4_place.py").read()
        self.assertIn("elif a.reference != a.h37rv_accession:", s)
        self.assertIn('--h37rv-accession "${H37RV_PATH%%#*}"',
                      open("bin/p4_place.sh").read())


if __name__ == "__main__":
    unittest.main(warnings="ignore")
