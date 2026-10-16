#!/usr/bin/env python3
"""Regression tests for fix SV-EVAL-1 (2026-10-16): an IS6110 stage-2 key
whose projection lands off a non-carrier's path is no longer ABSENT from the
graph alone.

On new1k_p200 that rule put 46 H37Rv-frame sites ABSENT in >= 90% of isolates:
at a ref_shared site the carrier's position sits beside the element its own
reference carries, so every reference without it is reached ~1,355 bp off
the path (4,230,303: dist 1,358 in all 997 non-carriers, one of which has
116x depth across the site). Now (is6110/bin/is6110_absent_check.py):

  H37Rv-frame key  the sample's P1 H37Rv alignment decides: IS6110 junction
                   reads -> ALT, spanning reads -> REF (NOCALL absent_masked
                   in the repeat mask), no reads at all -> ABS-1's region
                   test (ABSENT junction_joined when reads join R's
                   deletion), else NOCALL with a reason
  node-frame key   ABSENT node_not_on_r where R's path lacks the node,
                   NOCALL off_path_node_on_r where it carries it
  merged VCF       the stage-2 `ev` code is FORMAT/EV on IS6110 cells

The alignments are synthetic, written with pysam: a random 20 kb "H37Rv" with
one IS6110-like element at 10,001-11,355.

    $MTB_PY tests/test_audit_is6110_absent.py -v
"""
import csv
import gzip
import importlib.util
import os
import random
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
sys.path.insert(0, os.path.join(ROOT, "is6110", "bin"))

SAMTOOLS = os.environ.get("MTB_SAMTOOLS", "")
try:
    import pysam
except ImportError:                   # pragma: no cover
    pysam = None


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, rel))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


iac = load("is6110_absent_check", "is6110/bin/is6110_absent_check.py")
ic = load("inherited_check", "bin/inherited_check.py")

CTG = "NC_000962.3"
L = 20000
E0, E1 = 10000, 11355                 # element, 0-based half-open
P = 5000                              # the site: the element would sit after base P
rng = random.Random(7)
G = "".join(rng.choice("ACGT") for _ in range(L))
E = G[E0:E1]


def need_reads(tc):
    if pysam is None:
        tc.skipTest("pysam not importable; run with $MTB_PY")


def write_assets(d):
    with open(f"{d}/h37rv.fasta", "w") as fh:
        fh.write(f">{CTG}\n")
        for i in range(0, L, 80):
            fh.write(G[i:i + 80] + "\n")
    with open(f"{d}/is.gff", "w") as fh:
        fh.write(f"GCF_000195955#1#{CTG}\tRefSeq\tmobile_genetic_element\t{E0 + 1}"
                 f"\t{E1}\t.\t+\t.\tNote=IS6110-1\n")
    with open(f"{d}/mask.bed", "w") as fh:
        fh.write(f"{CTG}\t15000\t15200\n")
    return f"{d}/h37rv.fasta", f"{d}/is.gff", f"{d}/mask.bed"


def write_bam(path, reads):
    """reads: (start0, cigar tuples, query sequence[, name[, SA tag]])."""
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CTG, "LN": L}]}
    tmp = path + ".unsorted.bam"
    with pysam.AlignmentFile(tmp, "wb", header=hdr) as out:
        for i, r in enumerate(reads):
            a = pysam.AlignedSegment()
            a.query_name = r[3] if len(r) > 3 and r[3] else f"r{i}"
            a.reference_id = 0
            a.reference_start = r[0]
            a.mapping_quality = 60
            a.cigartuples = r[1]
            a.query_sequence = r[2]
            a.query_qualities = pysam.qualitystring_to_array("I" * len(r[2]))
            if len(r) > 4 and r[4]:
                a.set_tag("SA", r[4])
            a.flag = 0
            out.write(a)
    pysam.sort("-o", path, tmp)
    os.unlink(tmp)
    pysam.index(path)
    return path


def spanning(n, lo=P - 79, step=2):
    return [(lo + i * step, [(0, 100)], G[lo + i * step:lo + i * step + 100])
            for i in range(n)]


def right_junction(n):
    """60 bp of genome ending at P, then 40 bp of the element's start."""
    return [(P - 60, [(0, 60), (4, 40)], G[P - 60:P] + E[:40], f"jr{i}") for i in range(n)]


def left_junction(n):
    """40 bp of the element's end, then 60 bp of genome from P + 1."""
    return [(P, [(4, 40), (0, 60)], E[-40:] + G[P:P + 60], f"jl{i}") for i in range(n)]


def background():
    """Depth elsewhere, so the sample's median depth is not zero."""
    return [(p, [(0, 100)], G[p:p + 100]) for p in range(0, L - 100, 50)
            if not (P - 400 < p < P + 400)]


class Verdict(unittest.TestCase):
    """The counting rule, as measured (is6110_absent_check.verdict)."""

    def test_rules(self):
        self.assertEqual(iac.verdict(40, 0, 20), ("ALT", "absent_checked_alt"))
        self.assertEqual(iac.verdict(0, 96, 22), ("REF", "absent_checked_ref"))
        self.assertEqual(iac.verdict(1, 30, 22), ("REF", "absent_checked_ref"))
        self.assertEqual(iac.verdict(0, 0, 20), ("NOCALL", "absent_low_depth"))
        self.assertEqual(iac.verdict(2, 10, 20), ("NOCALL", "absent_low_depth"))
        self.assertEqual(iac.verdict(2, 30, 20), ("NOCALL", "absent_unresolved"))

    def test_junctions_against_many_spanning_reads_are_mixed(self):
        # ctrl25_p200 at 888,777: 50 junction vs 98 spanning reads, and the
        # assembly has no insertion there; the 0.2 share called it ALT
        self.assertEqual(iac.verdict(50, 98, 30), ("NOCALL", "absent_mixed"))
        self.assertEqual(iac.verdict(19, 54, 26), ("NOCALL", "absent_mixed"))


class SiteCheck(unittest.TestCase):
    """The P1 H37Rv alignment read at one site."""

    def setUp(self):
        need_reads(self)
        self.d = tempfile.mkdtemp()
        self.fa, self.gff, self.bed = write_assets(self.d)

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def chk(self, name, reads, mask=True):
        bam = write_bam(f"{self.d}/{name}.bam", reads + background())
        return iac.H37rvSiteCheck(bam, self.gff, self.fa,
                                  ic.Mask.from_bed(self.bed) if mask else None)

    def test_ref_shared_site_with_spanning_reads_is_ref(self):
        # the 4,230,303 case: off the path by ~the element's length, the
        # sample's reads straight across the site
        st, ev, j, s = self.chk("span", spanning(30)).check(P)
        self.assertEqual((st, ev, j), ("REF", "absent_checked_ref", 0))
        self.assertGreaterEqual(s, 20)

    def test_junction_reads_are_alt(self):
        st, ev, j, s = self.chk("junc", right_junction(6) + left_junction(6)).check(P)
        self.assertEqual((st, ev), ("ALT", "absent_checked_alt"))
        self.assertEqual((j, s), (12, 0))

    def test_deletion_edge_with_an_element_junction_is_alt(self):
        # the region right of P is gone in this sample, and the element sits
        # at its edge: one-sided junctions, nothing across, nothing beyond
        reads = right_junction(8) + [(p, [(0, 100)], G[p:p + 100])
                                     for p in range(P - 300, P - 100, 10)]
        st, ev, j, s = self.chk("edge", reads).check(P)
        self.assertEqual((st, ev, j, s), ("ALT", "absent_checked_alt", 8, 0))

    def test_clip_into_the_element_by_sa_tag_alone(self):
        # hard-clipped: no clipped bases to read, the SA tag says where the
        # far side continues -- the element's first base
        sa = f"{CTG},{E0 + 1},+,40M60H,60,0;"
        reads = [(P - 60, [(0, 60), (5, 40)], G[P - 60:P], f"h{i}", sa) for i in range(4)]
        st, ev, j, s = self.chk("sa", reads).check(P)
        self.assertEqual((st, j), ("ALT", 4))

    def test_clip_of_other_sequence_is_not_a_junction(self):
        other = "".join(rng.choice("ACGT") for _ in range(40))
        reads = [(P - 60, [(0, 60), (4, 40)], G[P - 60:P] + other, f"o{i}") for i in range(6)]
        self.assertEqual(self.chk("other", reads).counts(P)[0], 0)

    def test_low_depth_is_nocall(self):
        st, ev, _, _ = self.chk("low", spanning(2)).check(P)
        self.assertEqual((st, ev), ("NOCALL", "absent_low_depth"))

    def test_ref_inside_the_mask_is_nocall_masked(self):
        reads = [(p, [(0, 100)], G[p:p + 100]) for p in range(14960, 15040, 3)]
        st, ev, _, _ = self.chk("masked", reads).check(15050)
        self.assertEqual((st, ev), ("NOCALL", "absent_masked"))

    def test_missing_bam_is_fatal(self):
        with self.assertRaises(SystemExit):
            iac.H37rvSiteCheck(f"{self.d}/none.bam", self.gff, self.fa)


class RegionCheck(unittest.TestCase):
    """No reads at the site: ABS-1's test of R's deletion over it."""

    def setUp(self):
        need_reads(self)
        if not (SAMTOOLS and os.path.exists(SAMTOOLS)):
            self.skipTest("MTB_SAMTOOLS not set; source config/project_env.sh")
        self.d = tempfile.mkdtemp()
        self.fa, self.gff, self.bed = write_assets(self.d)
        # R deletes H37Rv 4,801-5,200, over the site
        with open(f"{self.d}/g.vcf", "w") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER"
                     "\tINFO\tFORMAT\tR\n")
            fh.write(f"{CTG}\t4800\t.\t{G[4799:5200]}\t{G[4799]}\t60\t.\t.\tGT\t1\n")
        self.pd = f"{self.d}/pd.tsv.gz"
        r = subprocess.run([sys.executable, "bin/panel_differences.py",
                            "--graph-vcf", f"{self.d}/g.vcf", "--out", self.pd],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def test_reads_joining_rs_deletion_make_it_absent(self):
        join = [(4740 + i, [(0, 60 - i), (2, 400), (0, 40 + i)],
                 G[4740 + i:4800] + G[5200:5240 + i]) for i in range(10)]
        bam = write_bam(f"{self.d}/del.bam", join + background())
        chk = iac.H37rvSiteCheck(bam, self.gff, self.fa, ic.Mask.from_bed(self.bed))
        chk.use_region_check(self.pd, "R", SAMTOOLS)
        st, ev, j, s = chk.check(P)
        self.assertEqual((st, ev, j, s), ("ABSENT", "junction_joined", 0, 0))

    def test_no_reads_and_no_join_stays_nocall(self):
        bam = write_bam(f"{self.d}/none.bam", background())
        chk = iac.H37rvSiteCheck(bam, self.gff, self.fa, ic.Mask.from_bed(self.bed))
        chk.use_region_check(self.pd, "R", SAMTOOLS)
        self.assertEqual(chk.check(P)[:2], ("NOCALL", "absent_low_depth"))


def tsv(path, cols, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})


class Stage2OffPath(unittest.TestCase):
    """is6110_p5_stage2.py --mode sample on off-path keys, end to end."""

    def setUp(self):
        need_reads(self)
        self.d = d = tempfile.mkdtemp()
        self.fa, self.gff, self.bed = write_assets(d)
        write_bam(f"{d}/B.h37rv.bam", spanning(30) + background())
        with open(f"{d}/refmap.tsv", "w") as fh:
            fh.write("sample\treference\nA\tC\nB\tR\n")
        tsv(f"{d}/keys.tsv", ["sample", "reference", "r_pos", "frame", "h37rv_pos",
                              "h37rv_state", "node"],
            [dict(sample="A", reference="C", r_pos=100, frame="h37rv",
                  h37rv_pos=P, h37rv_state="empty"),
             dict(sample="A", reference="C", r_pos=200, frame="h37rv",
                  h37rv_pos=E0, h37rv_state="occupied"),
             dict(sample="A", reference="C", r_pos=300, frame="node", node="77:5"),
             dict(sample="A", reference="C", r_pos=400, frame="node", node="88:0"),
             dict(sample="A", reference="C", r_pos=500, frame="node", node="99:0")])
        os.makedirs(f"{d}/s1"); os.makedirs(f"{d}/w/proj")
        self.keys = [f"h37rv:{P}:C><INS>", f"h37rv:{E0}:C><INS>",
                     "node:77:5:N><INS>", "node:88:0:N><INS>", "node:99:0:N><INS>"]
        with open(f"{d}/s1/B.tsv", "w") as fh:
            fh.write("sample\tkey\tstate\tallele\n")
            for k in self.keys:
                fh.write(f"B\t{k}\tNOCALL\t\n")
        with open(f"{d}/w/proj/R.tsv", "w") as fh:
            fh.write("short_key\tpos\tdist\n"
                     f"h37rv:{P}:\t9000\t1358\nh37rv:{E0}:\t9500\t3\n"
                     "node:77:5:\t100\t5\nnode:88:0:\t200\t7\nnode:99:0:\t300\t2\n")
        with open(f"{d}/node_paths.tsv", "w") as fh:
            fh.write("node\tn_paths\tpaths\n77\t1\tC\n88\t2\tC,R\n")

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def run_stage2(self, *extra, p1=True):
        d = self.d
        cmd = [sys.executable, "is6110/bin/is6110_p5_stage2.py", "--mode", "sample",
               "--sample", "B", "--refmap", f"{d}/refmap.tsv",
               "--cohort-keys", f"{d}/keys.tsv", "--stage1-dir", f"{d}/s1",
               "--isclean-dir", d, "--crossmap-dir", d, "--workdir", f"{d}/w",
               "--out", f"{d}/o.tsv", "--is6110-gff", self.gff,
               "--h37rv-fasta", self.fa, "--repeat-mask", self.bed,
               "--node-paths", f"{d}/node_paths.tsv",
               "--panel-diffs", f"{d}/pd.tsv.gz"]
        if p1:
            cmd += ["--p1-work", d]
        return subprocess.run(cmd + list(extra), capture_output=True, text=True)

    def rows(self):
        with open(f"{self.d}/w/samples/B.tsv", newline="") as fh:
            return {r["key"]: r for r in csv.DictReader(fh, delimiter="\t")}

    def test_off_path_keys(self):
        r = self.run_stage2()
        self.assertEqual(r.returncode, 0, r.stderr)
        got = {k: (v["state"], v["ev"]) for k, v in self.rows().items()}
        self.assertEqual(got[self.keys[0]], ("REF", "absent_checked_ref"))  # was ABSENT
        self.assertEqual(got[self.keys[1]], ("NOCALL", "absent_occupied"))
        self.assertEqual(got[self.keys[2]], ("ABSENT", "node_not_on_r"))
        self.assertEqual(got[self.keys[3]], ("NOCALL", "off_path_node_on_r"))
        self.assertEqual(got[self.keys[4]], ("NOCALL", "node_unknown"))

    def test_opt_out_keeps_the_old_absent_labelled(self):
        r = self.run_stage2("--no-absent-check", p1=False)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual((self.rows()[self.keys[0]]["state"],
                          self.rows()[self.keys[0]]["ev"]), ("ABSENT", "absent_unchecked"))

    def test_no_p1_work_is_fatal(self):
        r = self.run_stage2(p1=False)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("--p1-work", r.stderr)

    def test_missing_bam_is_fatal(self):
        os.remove(f"{self.d}/B.h37rv.bam")
        r = self.run_stage2()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("no P1 H37Rv alignment", r.stderr)


class MergedVcfEvidence(unittest.TestCase):
    """bin/merge_cohort_vcf.py writes the stage-2 code to FORMAT/EV."""

    COLS = ["sample", "key", "frame", "h37rv_pos", "h37rv_state", "node",
            "state", "evidence", "site_class"]

    def setUp(self):
        self.bgzip = os.environ.get("MTB_BGZIP", "")
        _t = os.path.join(os.path.dirname(self.bgzip), "tabix") if self.bgzip else ""
        self.tabix = os.environ.get("MTB_TABIX", _t)
        if not (self.bgzip and os.path.exists(self.bgzip) and os.path.exists(self.tabix)):
            self.skipTest("MTB_BGZIP not set; source config/project_env.sh")

    def merge(self, d, st_cols, st_rows):
        with open(f"{d}/m.tsv", "w") as fh:
            fh.write("key\tframe\tA\tB\tC\n")
        tsv(f"{d}/is.tsv", self.COLS,
            [dict(sample="A", key="h37rv:500", frame="h37rv", h37rv_pos=500,
                  h37rv_state="empty", state="ALT", evidence="two_sided",
                  site_class="ref_shared")])
        tsv(f"{d}/st2.tsv", st_cols, st_rows)
        r = subprocess.run(
            [sys.executable, "bin/merge_cohort_vcf.py", "--matrix", f"{d}/m.tsv",
             "--is6110-keys", f"{d}/is.tsv", "--is6110-states", f"{d}/st2.tsv",
             "--cohort-name", "t", "--ancestral", "", "--build-id", "testbuild",
             "--bgzip", self.bgzip, "--tabix", self.tabix, "--out", f"{d}/o.vcf.gz"],
            capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        with gzip.open(f"{d}/o.vcf.gz", "rt") as fh:
            return [l for l in fh if not l.startswith("#")][0].rstrip("\n").split("\t")

    def test_ev_written_and_checked_alt_accepted(self):
        cols = ["sample", "key", "state", "allele", "evidence", "ev"]
        rows = [dict(sample="A", key="h37rv:500:C><INS>", state="REF", ev="absent_checked_ref"),
                dict(sample="B", key="h37rv:500:C><INS>", state="ALT", ev="absent_checked_alt"),
                dict(sample="C", key="h37rv:500:C><INS>", state="NOCALL", ev="absent_masked")]
        with tempfile.TemporaryDirectory() as d:
            f = self.merge(d, cols, rows)
        self.assertEqual(f[8], "GT:ST:EV")
        # A's reported ALT is never overwritten by a computed state
        self.assertEqual(f[9:], ["1:ALT", "1:ALT:absent_checked_alt", ".:NOCALL:absent_masked"])

    def test_a_table_without_the_column_still_merges(self):
        cols = ["sample", "key", "state", "allele", "evidence"]
        rows = [dict(sample="B", key="h37rv:500:C><INS>", state="ABSENT"),
                dict(sample="C", key="h37rv:500:C><INS>", state="REF")]
        with tempfile.TemporaryDirectory() as d:
            f = self.merge(d, cols, rows)
        self.assertEqual(f[8], "GT:ST")
        self.assertEqual(f[9:], ["1:ALT", "2:ABSENT", "0:REF"])


if __name__ == "__main__":
    unittest.main(warnings="ignore")
