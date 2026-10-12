#!/usr/bin/env python3
"""Regression tests for wiring our SV caller into production (2026-10-12):
P4b reads P2's ours.events.tsv beside dysgu, applies the evaluated rules to
both callers, and P5 clusters our unsized insertions.

  - dysgu: FORMAT AF < 0.2 dropped, a record without AF kept
  - ours: BND rows are not events
  - both: tatC-helY DEL dropped in a carrier, kept in a non-carrier; the
    HaplotypeCaller hand-off drops a call at a well-supported 10-49 bp indel
  - an unsized INS of ours merges with dysgu's INS at the same place, and an
    assembled event of ours is the representative
  - `callable` marks calls touching R's uncallable mask
  - P5: an unsized INS clusters with sized ones on position, and does not
    become the cluster's size

Same pattern as tests/run_tests.py: small synthetic inputs, seconds to run.

    $MTB_PY tests/test_ours_wiring.py
"""
import csv
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import warnings

warnings.simplefilter("ignore", ResourceWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
sys.path.insert(0, os.path.join(ROOT, "bin"))

import p5_sv_matrix  # noqa: E402

S = "iso1"
DYSGU_HDR = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tiso1\n"


def dysgu_rec(pos, svtype, svlen, af, end=None):
    end = end or (pos if svtype == "INS" else pos + svlen)
    fmt, val = ("GT:AF", f"0/1:{af}") if af is not None else ("GT", "0/1")
    return (f"R\t{pos}\t.\tN\t<{svtype}>\t20\tPASS\t"
            f"SVTYPE={svtype};END={end};SVLEN={svlen};SR=3;PE=2\t{fmt}\t{val}\n")


class P4bOurs(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        j = lambda x: os.path.join(self.d, x)  # noqa: E731
        with open(j("dysgu.vcf"), "w") as f:
            f.write(DYSGU_HDR)
            f.write(dysgu_rec(10000, "DEL", 300, 0.1))     # AF below 0.2: dropped
            f.write(dysgu_rec(20000, "DEL", 300, None))    # no AF: kept
            f.write(dysgu_rec(30020, "INS", 1355, 0.6))    # merges with ours unsized INS
            f.write(dysgu_rec(50000, "DEL", 68, 0.3))      # tatC-helY DEL in a carrier
        ev = [("BND", 5000, 5000, 0, 12, "lclip;not assembled"),
              ("INS", 30000, 30000, 0, 20, "clip;one_end=0;not assembled"),
              ("INS", 30010, 30010, 1358, "", "assembled"),
              ("DEL", 40000, 40400, 400, "", "assembled"),       # in the mask
              ("INS", 60000, 60000, 0, 15, "clip;not assembled")]  # HC 25 bp indel
        with open(j("ours.tsv"), "w") as f:
            f.write("sample\ttype\tstart\tend\tsize\tsupport\tleft\tright\tdepth_ratio\tsource\n")
            for t, a, b, n, sup, src in ev:
                f.write(f"{S}\t{t}\t{a}\t{b}\t{n}\t{sup}\t\t\t\t{src}\n")
        with open(j("egt.tsv"), "w") as f:
            f.write("sample\treference\tsite_start\tsite_end\tref_has_element\tinside_depth\t"
                    "flank_depth\tratio\tgenotype\n")
            f.write(f"{S}\tR\t49990\t50058\t1\t30\t100\t0.3\tpresent\n")
        with gzip.open(j("hc.vcf.gz"), "wt") as f:
            f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tiso1\n")
            for p in range(1000, 1100):  # SNPs: median DP 100
                f.write(f"R\t{p}\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP\t1:0,100:100\n")
            f.write("R\t60010\t.\tA\t" + "A" + "C" * 25 + "\t50\tPASS\t.\tGT:AD:DP\t1:2,98:100\n")
        with open(j("mask.bed"), "w") as f:
            f.write("R\t40100\t40200\tparalog\n")
        with open(j("pos.tsv"), "w") as f:  # identity projection for every breakpoint
            for p in sorted({5000, 10000, 10300, 20000, 20300, 30000, 30010, 30020,
                             40000, 40400, 50000, 50068, 60000}):
                f.write(f"R,{p - 1},+\tH,{p - 1},+\t0\t+\n")
        r = subprocess.run([sys.executable, "bin/p4b_place_sv.py", "--sample", S, "--reference", "R",
                            "--build-id", "t", "--dysgu", j("dysgu.vcf"), "--ours", j("ours.tsv"),
                            "--element-gt", j("egt.tsv"), "--hc-vcf", j("hc.vcf.gz"),
                            "--callable-mask", j("mask.bed"), "--positions", j("pos.tsv"),
                            "--mask", j("mask.bed"), "--no-inherited", "--out", j("out.tsv")],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.log = r.stdout
        self.rows = list(csv.DictReader(open(j("out.tsv")), delimiter="\t"))

    def tearDown(self):
        if self.d.startswith(tempfile.gettempdir()):
            shutil.rmtree(self.d)

    def at(self, pos):
        return [r for r in self.rows if abs(int(r["r_pos"]) - pos) <= 50]

    def test_dysgu_af(self):
        self.assertEqual(self.at(10000), [])
        self.assertEqual(len(self.at(20000)), 1)

    def test_bnd_not_read(self):
        self.assertEqual(self.at(5000), [])

    def test_element_and_hc_rules(self):
        self.assertEqual(self.at(50000), [])
        self.assertEqual(self.at(60000), [])
        self.assertIn("element_dysgu 1", self.log)
        self.assertIn("hc_ours 1", self.log)

    def test_unsized_ins_merges_and_assembled_represents(self):
        ins = self.at(30000)
        self.assertEqual(len(ins), 1, ins)
        self.assertEqual(ins[0]["src"], "dysgu,ours")
        self.assertEqual((ins[0]["r_pos"], ins[0]["svlen"], ins[0]["qual_caller"]), ("30010", "1358", "ours"))

    def test_callable_flag(self):
        self.assertEqual(self.at(40000)[0]["callable"], "no")
        self.assertEqual(self.at(20000)[0]["callable"], "yes")


class P5Unsized(unittest.TestCase):
    def test_unsized_matches_on_position(self):
        a = dict(svtype="INS", pos=1000, svlen=0)
        b = dict(svtype="INS", pos=1100, svlen=1358)
        c = dict(svtype="INS", pos=1100, svlen=60)
        self.assertTrue(p5_sv_matrix.same_event(a, b))
        self.assertFalse(p5_sv_matrix.same_event(b, c))


if __name__ == "__main__":
    unittest.main()
