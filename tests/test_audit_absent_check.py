#!/usr/bin/env python3
"""Regression tests for HANDOFF fix list ABS-1 and ABS-2 (2026-10-14).

ABS-1: an ABSENT that the matched reference R's deletion explains is checked
against the sample's own P1 alignment to H37Rv, at the ends of R's deletion;
a masked site is NOCALL "masked". ABS-2: an inherited SV is read in BOTH
frames, insertions included; masked and tandem events are NOCALL with their
reason; P4b's same-event suppression and masked label; the merged VCF's
FORMAT/EV.

The fixtures are synthetic alignments written with pysam, plus the 122
calibration regions of analysis/inherited_absent/ (tests/fixtures/
absent_calibration.tsv: their recorded junction counts and the SPAdes
reassembly's verdict), and the ctpV geometry of analysis/ctpV_check.md.

    $MTB_PY tests/test_audit_absent_check.py -v
"""
import csv
import gzip
import hashlib
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


ic = load("inherited_check", "bin/inherited_check.py")


def need_reads(tc):
    if pysam is None:
        tc.skipTest("pysam not importable; run with $MTB_PY")
    if not (SAMTOOLS and os.path.exists(SAMTOOLS)):
        tc.skipTest("MTB_SAMTOOLS not set; source config/project_env.sh")


def write_fasta(path, name, seq):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i + 80] + "\n")
    return path


def write_bam(path, ctg, length, reads, seq):
    """reads: (start0, cigar tuples, mapq[, name]); sequence taken from `seq`
    along the alignment (soft clips filled with N)."""
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"},
           "SQ": [{"SN": ctg, "LN": length}]}
    tmp = path + ".unsorted.bam"
    with pysam.AlignmentFile(tmp, "wb", header=hdr) as out:
        for i, r in enumerate(reads):
            start, cig, mq = r[:3]
            a = pysam.AlignedSegment()
            a.query_name = r[3] if len(r) > 3 else f"r{i}"
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = mq
            a.cigartuples = cig
            q, p = [], start
            for op, n in cig:
                if op == 0:
                    q.append(seq[p:p + n]); p += n
                elif op == 2:
                    p += n
                elif op in (1, 4):
                    q.append("N" * n)
            a.query_sequence = "".join(q)
            a.query_qualities = pysam.qualitystring_to_array("I" * len(a.query_sequence))
            a.flag = 0
            out.write(a)
    pysam.sort("-o", path, tmp)
    os.unlink(tmp)
    pysam.index(path)
    return path


def tiled(lo, hi, step=5, rl=100, mq=60):
    """Reads of length rl, gapless, every `step` bases from lo to hi (0-based)."""
    return [(p, [(0, rl)], mq) for p in range(lo, hi - rl, step)]


# ---------------------------------------------------------------------------
class CalibrationRegions(unittest.TestCase):
    """The 122 regions SPAdes reassembly decided (analysis/inherited_absent,
    job 51392271): the junction rule applied to their recorded read counts
    must give the assembly's verdict at every one -- the agreement the
    calibration measured, pinned so a change to the rule cannot lose it."""

    def test_rule_reproduces_every_assembly_verdict(self):
        rows = list(csv.DictReader(
            (l for l in open("tests/fixtures/absent_calibration.tsv") if not l.startswith("#")),
            delimiter="\t"))
        self.assertEqual(len(rows), 122)
        bad = []
        for r in rows:
            mb = max(5.0, 0.25 * float(r["baseline_depth"]))
            v = ic.junction_verdict(int(r["bridge_start"]), int(r["bridge_end"]),
                                    int(r["joined"]), mb)
            if v != r["assembly_verdict"]:
                bad.append((r["event"], v, r["assembly_verdict"]))
        self.assertEqual(bad, [])

    def test_ctpv_regions_are_present_in_samea2297133(self):
        rows = [r for r in csv.DictReader(
            (l for l in open("tests/fixtures/absent_calibration.tsv") if not l.startswith("#")),
            delimiter="\t") if r["sample"] == "SAMEA2297133"]
        self.assertEqual({(r["start"], r["end"]) for r in rows},
                         {("1078522", "1078756"), ("1078811", "1080106")})
        self.assertTrue(all(r["assembly_verdict"] == "present_here" for r in rows))


# ---------------------------------------------------------------------------
class PanelDifferences(unittest.TestCase):
    """P0's per-panel-genome table (bin/panel_differences.py) and its reader."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.vcf = os.path.join(self.d, "g.vcf")
        ref = "A" + "CGTTACGGAT" * 2                   # 21 bp at 100
        with open(self.vcf, "w") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER"
                     "\tINFO\tFORMAT\tA\tB\tC\n")
            fh.write(f"H\t100\t.\t{ref}\tA\t60\t.\t.\tGT\t1\t.\t0\n")
            fh.write("H\t500\t.\tG\tG" + "T" * 60 + ",GA\t60\t.\t.\tGT\t0\t1\t2\n")
            fh.write("H\t700\t.\tC\tT\t60\t.\t.\tGT\t1\t1\t0\n")
        self.out = os.path.join(self.d, "pd.tsv.gz")
        r = subprocess.run([sys.executable, "bin/panel_differences.py",
                            "--graph-vcf", self.vcf, "--out", self.out],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def test_rows(self):
        got = [l.rstrip("\n").split("\t") for l in gzip.open(self.out, "rt")
               if not l.startswith("#") and not l.startswith("accession")]
        self.assertIn(["A", "DEL", "101", "120", "20", "100", "21", "1"], got)
        self.assertIn(["B", "NOGT", "100", "120", "21", "100", "21", "0"], got)
        self.assertIn(["B", "INS", "500", "500", "60", "500", "1", "1"], got)
        # a 1 bp insertion and SNPs are not listed
        self.assertFalse([g for g in got if g[0] == "C"])
        self.assertEqual(len(got), 3)

    def test_reader(self):
        a = ic.PanelDiffs(self.out, "A", self.vcf)
        self.assertEqual(a.deletion_over(110, 110)[:2], (101, 120))
        self.assertIsNone(a.deletion_over(121, 130))
        b = ic.PanelDiffs(self.out, "B", self.vcf)
        self.assertTrue(b.no_genotype_over(105, 106))
        self.assertEqual(b.ins, {(500, 60): 1})
        # H37Rv itself is the VCF's reference path: no differences, not an error
        h = ic.PanelDiffs(self.out, ic.H37RV_ACC)
        self.assertEqual(h.dels, [])

    def test_refuses_another_graph_vcf_or_genome(self):
        other = os.path.join(self.d, "other.vcf")
        with open(other, "w") as fh:
            fh.write("##fileformat=VCFv4.2\n")
        with self.assertRaises(SystemExit):
            ic.PanelDiffs(self.out, "A", other)
        with self.assertRaises(SystemExit):
            ic.PanelDiffs(self.out, "NOT_IN_PANEL")
        with self.assertRaises(SystemExit):
            ic.PanelDiffs(os.path.join(self.d, "missing.tsv.gz"), "A")


class MaskRule(unittest.TestCase):
    """ABS-1 item 2: the site, or either end of R's deletion within 150 bp.
    The end rule was dropped and restored on 2026-10-14: span alone was
    523 right / 187 wrong on ctrl25_p200, against 472 / 11 with it."""

    def test_site_and_ends(self):
        with tempfile.TemporaryDirectory() as d:
            bed = os.path.join(d, "m.bed")
            with open(bed, "w") as fh:
                fh.write("c\t999\t1100\tpe_ppe|PE1|x\n")       # 1-based 1000..1100
            m = ic.Mask.from_bed(bed)
        self.assertTrue(m.site(1050, 1050))                       # in it
        self.assertTrue(m.site(900, 1000))                        # span touches it
        self.assertFalse(m.site(2000, 2000))
        # R's deletion 1200..3000: its start is 100 bp from the mask
        self.assertTrue(m.site(2000, 2000, (1200, 3000)))
        # 151 bp away is outside a read length
        self.assertFalse(m.site(2000, 2000, (1252, 3000)))
        self.assertTrue(m.site(2000, 2000, (1250, 3000)))


# ---------------------------------------------------------------------------
class H37rvJunction(unittest.TestCase):
    """The H37Rv-frame check of R's deletion: reads across both ends, reads
    joining them, or neither (junction_check.py's rule)."""

    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d, ignore_errors=True)

    def setUp(self):
        need_reads(self)
        rng = random.Random(11)
        self.seq = "".join(rng.choice("ACGT") for _ in range(6000))

    def frame(self, name, reads, vcf=None):
        bam = write_bam(os.path.join(self.d, name + ".bam"), "c", 6000, reads, self.seq)
        return ic.H37rvFrame(bam, SAMTOOLS, vcf=vcf, h37=self.seq)

    def test_present_here(self):
        f = self.frame("present", tiled(1500, 3000))
        v, ev = f.junction(2001, 2300)
        self.assertEqual(v, "present_here", ev)
        self.assertGreaterEqual(ev["bridge_start"], 5)

    def test_deleted_here(self):
        # 10 reads joining 2000|2301 with a 300 bp CIGAR deletion, nothing across
        reads = [(1950 - i, [(0, 50 + i), (2, 300), (0, 50)], 60) for i in range(10)]
        f = self.frame("deleted", reads + tiled(1200, 1990) + tiled(2310, 3200))
        v, ev = f.junction(2001, 2300)
        self.assertEqual(v, "deleted_here", ev)
        self.assertGreaterEqual(ev["joined"], 3)

    def test_unresolved_and_low_mapq(self):
        # reads across the junction at MAPQ 10 only (a repeat): they do not
        # count; MAPQ 60 reads elsewhere give the sample a depth at all
        f = self.frame("lowmq", tiled(1500, 3000, mq=10) + tiled(4000, 5000))
        self.assertEqual(f.junction(2001, 2300)[0], "unresolved")
        # a sample with no MAPQ >= 20 depth anywhere is broken, not unresolved
        with self.assertRaises(SystemExit):
            self.frame("nodepth", tiled(1500, 3000, mq=10))

    def test_small_call(self):
        vcf = os.path.join(self.d, "p1.vcf")
        b = self.seq[2099]
        alt = "A" if b != "A" else "C"
        with open(vcf, "w") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER"
                     f"\tINFO\tFORMAT\tS\nc\t2100\t.\t{b}\t{alt}\t50\t.\t.\tGT\t1\n")
        f = self.frame("calls", tiled(1500, 3000), vcf=vcf)
        self.assertEqual(f.small_call(f"h37rv:2100:{b}>{alt}", 2100, 2100, 5)[0], "ALT")
        third = [x for x in "ACGT" if x not in (b, alt)][0]
        self.assertEqual(f.small_call(f"h37rv:2100:{b}>{third}", 2100, 2100, 5),
                         ("NOCALL", "present_other_allele"))
        self.assertEqual(f.small_call("h37rv:2200:X>Y", 2200, 2200, 5)[0], "REF")
        self.assertEqual(f.small_call("h37rv:4500:X>Y", 4500, 4500, 5),
                         ("NOCALL", "present_low_depth"))

    def test_missing_bam_is_fatal(self):
        with self.assertRaises(SystemExit):
            ic.H37rvFrame(os.path.join(self.d, "none.bam"), SAMTOOLS)


# ---------------------------------------------------------------------------
class CtpVCases(unittest.TestCase):
    """analysis/ctpV_check.md, in H37Rv coordinates. Every L1.2.1 sample was
    matched to a reference carrying the 235 bp deletion 1,078,522-1,078,756
    and the 1,296 bp deletion 1,078,811-1,080,106 of ctpV. SAMEA2297133 has an
    intact gene (72-107x, no clip or split read at either end) and was ABSENT
    at small-variant keys 1,078,644 and 1,078,715 and at the nested 872 bp
    interval; SAMEA5542103 carries the deletions (30-65 junction reads per
    breakpoint)."""

    L = 1_081_500

    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d, ignore_errors=True)

    def setUp(self):
        need_reads(self)
        rng = random.Random(5)
        self.seq = "".join(rng.choice("ACGT") for _ in range(self.L))
        d = self.d
        # R's two deletions, as the panel table states them
        vcf = os.path.join(d, "g.vcf")
        with open(vcf, "w") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\t"
                     "FILTER\tINFO\tFORMAT\tR\n")
            fh.write(f"H\t1078521\t.\t{self.seq[1078520:1078756]}\t"
                     f"{self.seq[1078520]}\t60\t.\t.\tGT\t1\n")
            fh.write(f"H\t1078810\t.\t{self.seq[1078809:1080106]}\t"
                     f"{self.seq[1078809]}\t60\t.\t.\tGT\t1\n")
        self.vcf = vcf
        pd = os.path.join(d, "pd.tsv.gz")
        subprocess.run([sys.executable, "bin/panel_differences.py", "--graph-vcf",
                        vcf, "--out", pd], check=True, capture_output=True)
        self.diffs = ic.PanelDiffs(pd, "R", vcf)
        bed = os.path.join(d, "m.bed")
        with open(bed, "w") as fh:
            fh.write("c\t100\t200\tpe_ppe|x|y\n")           # far away
        self.mask = ic.Mask.from_bed(bed)

    def intact(self):
        return ic.H37rvFrame(write_bam(os.path.join(self.d, "intact.bam"), "c",
                                       self.L, tiled(1077800, 1081300, step=2),
                                       self.seq), SAMTOOLS)

    def carrier(self):
        # reads joining each deletion's ends, none across them
        r = []
        for i in range(30):
            r.append((1078470 - i, [(0, 51 + i), (2, 235), (0, 40)], 60))
            r.append((1078760 - i, [(0, 50 + i), (2, 1296), (0, 40)], 60))
        r += tiled(1077800, 1078500, step=2) + tiled(1080110, 1081300, step=2)
        return ic.H37rvFrame(write_bam(os.path.join(self.d, "carrier.bam"), "c",
                                       self.L, r, self.seq), SAMTOOLS)

    def test_samea2297133_small_keys_present(self):
        f = self.intact()
        for p in (1078644, 1078715, 1079000):
            v, d, _ = ic.check_absent(self.diffs, self.mask, f, p, p)
            self.assertEqual(v, "present_here", p)
            self.assertEqual(f.small_call(f"h37rv:{p}:A>C", p, p, 5)[0], "REF")

    def test_samea5542103_stays_absent(self):
        f = self.carrier()
        for p in (1078644, 1079000):
            v, d, ev = ic.check_absent(self.diffs, self.mask, f, p, p)
            self.assertEqual(v, "deleted_here", (p, ev))

    def test_masked_end_is_nocall(self):
        bed = os.path.join(self.d, "m2.bed")
        with open(bed, "w") as fh:
            # 100 bp before R's 1,296 bp deletion starts
            fh.write("c\t1078600\t1078711\tpe_ppe|x|y\n")
        m = ic.Mask.from_bed(bed)
        # a site well inside the deletion is masked through the deletion's end
        v, _, _ = ic.check_absent(self.diffs, m, self.intact(), 1079500, 1079500)
        self.assertEqual(v, "masked")

    def test_no_r_deletion(self):
        v, d, _ = ic.check_absent(self.diffs, self.mask, self.intact(), 1081000, 1081000)
        self.assertEqual((v, d), ("no_r_deletion", None))


# ---------------------------------------------------------------------------
class MatchedFrame(unittest.TestCase):
    """ABS-2: R's junction read in the sample's alignment to R, and the
    two-frame rule."""

    def test_combine_frames(self):
        c = ic.combine_frames
        self.assertEqual(c("ALT", "carries"), ("ALT", "two_frames"))
        self.assertEqual(c("REF", "lacks"), ("REF", "two_frames"))
        self.assertEqual(c("NOCALL", "carries"), ("ALT", "matched_only"))
        self.assertEqual(c("REF", "unclear"), ("REF", "h37rv_only"))
        self.assertEqual(c("ALT", "lacks"), ("NOCALL", "frames_conflict"))
        self.assertEqual(c("NOCALL", "low_depth"), ("NOCALL", "neither_frame"))

    def test_tandem_flag(self):
        # a MIRU-like unit: the inserted sequence repeats its flank (MIRU40 is
        # R's extra unit at 802,427 in the RESULTS.md examples)
        rng = random.Random(3)
        unit = "".join(rng.choice("ACGT") for _ in range(77))
        left = "".join(rng.choice("ACGT") for _ in range(200)) + unit * 2
        right = unit + "".join(rng.choice("ACGT") for _ in range(200))
        h37 = left + right
        vp = len(left)
        q, a, b, chg = ic.local_haplotype(h37, vp, h37[vp - 1], h37[vp - 1] + unit)
        self.assertTrue(ic.is_tandem(q, a, b, chg))
        other = "".join(rng.choice("ACGT") for _ in range(77))
        q, a, b, chg = ic.local_haplotype(h37, vp, h37[vp - 1], h37[vp - 1] + other)
        self.assertFalse(ic.is_tandem(q, a, b, chg))

    def test_deletion_junction_carries_and_lacks(self):
        need_reads(self)
        rng = random.Random(9)
        r = "".join(rng.choice("ACGT") for _ in range(3000))
        with tempfile.TemporaryDirectory() as d:
            # carries: reads straight across R's junction at 1500
            bam = write_bam(os.path.join(d, "c.bam"), "R", 3000, tiled(1300, 1700), r)
            m = ic.matched_junction(pysam.AlignmentFile(bam), "R", 1500, 3000)
            self.assertEqual(m["call"], "carries", m)
            # lacks: reads soft-clipped at the junction (their other part is
            # the sequence R lacks)
            reads = [(1400 + i, [(0, 100 - i), (4, 50 + i)], 60) for i in range(0, 40, 4)]
            reads += [(1500, [(4, 50), (0, 100)], 60) for _ in range(6)]
            bam = write_bam(os.path.join(d, "l.bam"), "R", 3000, reads, r)
            m = ic.matched_junction(pysam.AlignmentFile(bam), "R", 1500, 3000)
            self.assertEqual(m["call"], "lacks", m)

    def test_insertion_carried_and_lacking(self):
        need_reads(self)
        rng = random.Random(4)
        r = "".join(rng.choice("ACGT") for _ in range(4000))
        with tempfile.TemporaryDirectory() as d:
            # R's Y = 1800..2000; carried: even depth, clean junctions
            bam = write_bam(os.path.join(d, "c.bam"), "R", 4000, tiled(1300, 2600), r)
            m = ic.matched_insertion(pysam.AlignmentFile(bam), "R", 1800, 2000, 4000)
            self.assertEqual(m["call"], "carries", m)
            # lacking: no reads across Y, reads clipped at both of its ends
            reads = tiled(1300, 1795) + tiled(2005, 2600)
            reads += [(1700, [(0, 100), (4, 40)], 60) for _ in range(4)]
            reads += [(2000, [(4, 40), (0, 100)], 60) for _ in range(4)]
            bam = write_bam(os.path.join(d, "l.bam"), "R", 4000, reads, r)
            m = ic.matched_insertion(pysam.AlignmentFile(bam), "R", 1800, 2000, 4000)
            self.assertEqual(m["call"], "lacks", m)
            # a GC-rich PE_PGRS insert the isolate CARRIES can show near-zero
            # depth and no clips (150,890 in SAMD00089643): not "lacks"
            reads = tiled(1300, 1795) + tiled(2005, 2600)
            bam = write_bam(os.path.join(d, "g.bam"), "R", 4000, reads, r)
            m = ic.matched_insertion(pysam.AlignmentFile(bam), "R", 1800, 2000, 4000)
            self.assertNotEqual(m["call"], "lacks", m)
            self.assertNotEqual(m["call"], "carries", m)


class InheritedInsertionStates(unittest.TestCase):
    """ABS-2 items 2 and 3: an inherited insertion's state, and the reasons."""

    def test_states(self):
        m = load("p5_inherited_sv", "bin/p5_inherited_sv.py")
        self.assertEqual(m.ins_state(False, False, "carries"), ("ALT", "matched_carries"))
        # "lacks" was wrong 11 of 12 on ctrl25_p200: NOCALL (user, 2026-10-14)
        self.assertEqual(m.ins_state(False, False, "lacks"),
                         ("NOCALL", "ins_lacks_unvalidated"))
        self.assertEqual(m.ins_state(True, False, "carries"), ("NOCALL", "masked"))
        self.assertEqual(m.ins_state(False, True, "carries"), ("NOCALL", "tandem"))
        for c in ("repeat", "low_depth", "unclear", "unplaced"):
            self.assertEqual(m.ins_state(False, False, c), ("NOCALL", c))


# ---------------------------------------------------------------------------
class SvGenotypeBothFrames(unittest.TestCase):
    """p5_sv_genotype.py interval mode with p5_inherited_sv.py's table: the
    ABS-2 rule for an inherited deletion, and ABS-1 for an interval absent
    from R (the ctpV 872 bp case)."""

    IV = ["svi:DEL:1001:100", "svi:DEL:2001:100", "svi:DEL:3001:100",
          "svi:DEL:4001:100", "svi:DEL:5001:100", "svi:DEL:7001:100",
          "svi:DEL:7301:100", "svi:DEL:7601:100"]

    def genotype(self, d, chk_rows, extra=()):
        def tsv(path, header, recs):
            with open(path, "w") as fh:
                fh.write("\t".join(header) + "\n")
                for r in recs:
                    fh.write("\t".join(str(r.get(c, "")) for c in header) + "\n")
        ivs = []
        for iv in self.IV:
            s = int(iv.split(":")[2])
            # 1-5: R carries them; 6-8 lie inside R's 7000-7999 deletion
            ivs.append(dict(interval=iv, source="graph", svtype="DEL", start=s,
                            end=s + 99, svlen=100, support_tier="A_multi_assembly",
                            is6110_prox=0, n_ref_carriers=1,
                            ref_carriers="R" if s < 7000 else "X"))
        tsv(f"{d}/iv.tsv", ["interval", "source", "svtype", "start", "end", "svlen",
                            "support_tier", "is6110_prox", "n_ref_carriers",
                            "ref_carriers"], ivs)
        with open(f"{d}/proj.pos", "w") as fh:
            for p in range(1, 9001):
                if 6800 <= p <= 8100:
                    continue                  # R lacks 6800-8100: nothing projects
                fh.write(f"H#1#c,{p - 1},+\tR#1#c,{p - 1},+\t0\t+\t+\n")
        with open(f"{d}/g.vcf", "w") as fh:
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
            fh.write("c\t1\t.\tA\t<NON_REF>\t.\t.\tEND=9000\tGT:DP\t0:30\n")
        h = {"svi:DEL:1001:100": "REF", "svi:DEL:2001:100": "ALT",
             "svi:DEL:3001:100": "NOCALL", "svi:DEL:4001:100": "ALT",
             "svi:DEL:5001:100": "ALT", "svi:DEL:7001:100": "REF",
             "svi:DEL:7301:100": "REF", "svi:DEL:7601:100": "REF"}
        tsv(f"{d}/tf.tsv", ["sample", "interval", "h_state", "h_clip_sides"],
            [dict(sample="S", interval=k, h_state=v, h_clip_sides=0) for k, v in h.items()])
        tsv(f"{d}/chk.tsv", ["sample", "kind", "id", "masked", "tandem", "call",
                             "h_junction"], chk_rows)
        r = subprocess.run(
            [sys.executable, "bin/p5_sv_genotype.py", "--intervals", f"{d}/iv.tsv",
             "--reference", "R", "--sample", "S", "--gvcf", f"{d}/g.vcf",
             "--projected", f"{d}/proj.pos", "--twoframe", f"{d}/tf.tsv",
             "--mapq-scope", "never", *extra, "--out", f"{d}/out.tsv"],
            capture_output=True, text=True)
        return r

    def chk(self):
        base = dict(sample="S", kind="DEL", masked="no", tandem="no")
        m = {"svi:DEL:1001:100": dict(call="carries"),           # H37Rv REF
             "svi:DEL:2001:100": dict(call="carries"),           # H37Rv ALT
             "svi:DEL:3001:100": dict(call="carries"),           # H37Rv NOCALL
             "svi:DEL:4001:100": dict(call="carries", masked="yes"),
             "svi:DEL:5001:100": dict(call="carries", tandem="yes"),
             "svi:DEL:7001:100": dict(h_junction="present_here"),
             "svi:DEL:7301:100": dict(h_junction="deleted_here"),
             "svi:DEL:7601:100": dict(h_junction="unresolved")}
        return [dict(base, id=k, **v) for k, v in m.items()]

    def test_rules(self):
        with tempfile.TemporaryDirectory() as d:
            r = self.genotype(d, self.chk(), ["--inherited-check", f"{d}/chk.tsv"])
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            got = {q["interval"]: (q["state"], q["evidence"])
                   for q in csv.DictReader(open(f"{d}/out.tsv"), delimiter="\t")}
        self.assertEqual(got["svi:DEL:1001:100"],
                         ("NOCALL", "inherited_reference_frames_conflict"))
        self.assertEqual(got["svi:DEL:2001:100"], ("ALT", "inherited_reference_two_frames"))
        self.assertEqual(got["svi:DEL:3001:100"], ("ALT", "inherited_reference_matched_only"))
        # masked inherited deletions are genotyped by both frames (2026-10-14)
        self.assertEqual(got["svi:DEL:4001:100"],
                         ("ALT", "inherited_reference_two_frames_inmask"))
        self.assertEqual(got["svi:DEL:5001:100"], ("NOCALL", "inherited_reference_tandem"))
        # ABS-1 for intervals absent from R (flanks inside R's gap too)
        self.assertEqual(got["svi:DEL:7001:100"], ("REF", "absent_checked_present_ref"))
        self.assertEqual(got["svi:DEL:7301:100"], ("ABSENT", "absent_junction_joined"))
        self.assertEqual(got["svi:DEL:7601:100"], ("NOCALL", "absent_unresolved"))

    def test_table_required_and_complete(self):
        with tempfile.TemporaryDirectory() as d:
            r = self.genotype(d, self.chk())
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("--inherited-check", r.stderr)
            r = self.genotype(d, self.chk()[:-1], ["--inherited-check", f"{d}/chk.tsv"])
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("another catalogue", r.stderr)
            r = self.genotype(d, self.chk(), ["--no-inherited-check"])
            self.assertEqual(r.returncode, 0, r.stderr)
            got = {q["interval"]: q["state"]
                   for q in csv.DictReader(open(f"{d}/out.tsv"), delimiter="\t")}
            # as before: the H37Rv frame alone, and ABSENT unread
            self.assertEqual(got["svi:DEL:1001:100"], "REF")
            self.assertEqual(got["svi:DEL:7001:100"], "ABSENT")


# ---------------------------------------------------------------------------
class P4bInheritedRows(unittest.TestCase):
    """ABS-2 item 4: a called SV suppresses an inherited one only when it is
    the same type within 50% in size; inherited rows carry a masked label."""

    def run_p4b(self, d, calls):
        def w(name, text):
            p = os.path.join(d, name)
            with open(p, "w") as fh:
                fh.write(text)
            return p
        hdr = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        dys = w("d.vcf", hdr + "".join(
            f"R\t{p}\t.\tN\t<{t}>\t60\tPASS\tSVTYPE={t};END={p + (ln if t == 'DEL' else 1)};"
            f"SVLEN={ln}\tGT:AF\t1:0.5\n" for p, t, ln in calls))
        ours = w("o.tsv", "sample\ttype\tstart\tend\tsize\tsupport\tsource\n")
        egt = w("e.tsv", "sample\tstate\n")
        subprocess.run(["bash", "-c", f"printf '{hdr}' | gzip > {d}/hc.vcf.gz"], check=True)
        cm = w("cm.bed", "")
        pos = w("p.pos", "".join(
            f"R#1#c,{q - 1},+\tH#1#c,{q - 1},+\t0\t+\t+\n"
            for p, t, ln in calls for q in (p, p + (ln if t == "DEL" else 1))))
        mask = w("m.bed", "c\t4999\t5100\tpe_ppe|x|y\n")
        ref = "A" + "C" * 400
        g = w("g.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER"
              "\tINFO\tFORMAT\tR\n"
              f"H\t1000\t.\t{ref}\tA\t60\t.\tLV=0\tGT\t1\n"
              f"H\t5050\t.\tA\t{'A' + 'G' * 80}\t60\t.\tLV=0\tGT\t1\n")
        r = subprocess.run([sys.executable, "bin/p4b_place_sv.py", "--sample", "S",
                            "--reference", "R", "--build-id", "b", "--dysgu", dys,
                            "--ours", ours, "--element-gt", egt,
                            "--hc-vcf", f"{d}/hc.vcf.gz", "--callable-mask", cm,
                            "--positions", pos, "--mask", mask, "--graph-vcf", g,
                            "--out", f"{d}/out.tsv"], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        return [x for x in csv.DictReader(open(f"{d}/out.tsv"), delimiter="\t")
                if x["component"] == "inherited"]

    def test_other_type_nearby_does_not_suppress(self):
        with tempfile.TemporaryDirectory() as d:
            inh = self.run_p4b(d, [(1100, "INS", 60)])     # an INS beside R's 400 bp DEL
            self.assertIn("DEL", [x["svtype"] for x in inh])

    def test_same_event_suppresses(self):
        with tempfile.TemporaryDirectory() as d:
            inh = self.run_p4b(d, [(1050, "DEL", 380)])
            self.assertNotIn("DEL", [x["svtype"] for x in inh])

    def test_masked_label(self):
        with tempfile.TemporaryDirectory() as d:
            inh = {x["svtype"]: x for x in self.run_p4b(d, [])}
            self.assertEqual(inh["DEL"]["callable"], "yes")
            self.assertEqual(inh["INS"]["callable"], "no")      # in the mask


# ---------------------------------------------------------------------------
class P5StatesAbsentCheck(unittest.TestCase):
    """p5_states.py end to end: an ABSENT R's deletion explains is checked
    against the P1 alignment, and every ABSENT/NOCALL carries its reason."""

    def setUp(self):
        need_reads(self)

    def test_absent_replaced_by_the_sample_call_and_reasons_written(self):
        rng = random.Random(7)
        h = "".join(rng.choice("ACGT") for _ in range(4_000_100))
        r = h[:2000] + h[2300:]                       # R lacks 2001..2300
        d = tempfile.mkdtemp()
        try:
            os.makedirs(f"{d}/refs")
            write_fasta(f"{d}/refs/GCF_000195955.fasta", "GCF_000195955", h)
            write_fasta(f"{d}/refs/RREF.fasta", "RREF", r)
            other = "A" if h[2099] != "A" else "C"
            keys = [dict(key=f"h37rv:2100:{h[2099]}>{other}", frame="h37rv",
                         h37rv_pos=2100, canonical_ref=h[2099], canonical_alt=other,
                         kind="SNP", region="core")]
            with open(f"{d}/keys.tsv", "w") as fh:
                cols = ["key", "frame", "h37rv_pos", "node", "node_offset",
                        "canonical_ref", "canonical_alt", "region", "kind", "acc_locus"]
                fh.write("\t".join(cols) + "\n")
                for k in keys:
                    fh.write("\t".join(str(k.get(c, "")) for c in cols) + "\n")
            with open(f"{d}/placed.tsv", "w") as fh:
                fh.write("sample\tframe\th37rv_pos\tref\talt\tr_pos\tnode\tnode_offset\n")
            with open(f"{d}/proj.pos", "w") as fh:
                # odgi lands a deleted base beside the gap, at the deleted length
                fh.write(f"H#1#c,2099,+\tRREF#1#c,2000,+\t300\t+\t+\n")
            with open(f"{d}/s.g.vcf", "w") as fh:
                fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\t"
                         "FILTER\tINFO\tFORMAT\tS\n")
                fh.write(f"c\t1\t.\t{r[0]}\t<NON_REF>\t.\t.\tEND={len(r)}\tGT:DP\t0:30\n")
            with open(f"{d}/frames.tsv", "w") as fh:
                fh.write("accession\tstrand\toffset\tpanel_len\nRREF\t+\t0\t4000000\n")
            g = f"{d}/g.vcf"
            with open(g, "w") as fh:
                fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\t"
                         "FILTER\tINFO\tFORMAT\tRREF\n")
                fh.write(f"H\t2000\t.\t{h[1999:2300]}\t{h[1999]}\t60\t.\t.\tGT\t1\n")
            subprocess.run([sys.executable, "bin/panel_differences.py", "--graph-vcf",
                            g, "--out", f"{d}/pd.tsv.gz"], check=True, capture_output=True)
            with open(f"{d}/m.bed", "w") as fh:
                fh.write("c\t100\t200\tpe_ppe|x|y\n")
            # the sample HAS the sequence: gapless reads across both ends
            bam = write_bam(f"{d}/p1.bam", "c", len(h), tiled(1500, 2900), h)
            with open(f"{d}/p1.vcf", "w") as fh:
                fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\t"
                         "FILTER\tINFO\tFORMAT\tS\n")
            base = [sys.executable, "bin/p5_states.py", "--sample", "S",
                    "--reference", "RREF", "--keys", f"{d}/keys.tsv",
                    "--placed", f"{d}/placed.tsv", "--gvcf", f"{d}/s.g.vcf",
                    "--projected", f"{d}/proj.pos",
                    "--h37rv", f"{d}/refs/GCF_000195955.fasta",
                    "--node-paths", f"{d}/np.tsv", "--node-positions", f"{d}/npos.tsv"]
            with open(f"{d}/np.tsv", "w") as fh:
                fh.write("node\tpaths\n")
            with open(f"{d}/npos.tsv", "w") as fh:
                fh.write("node\taccession\tstart\tstrand\tn_occurrences\tlength\n")
            env = dict(os.environ, MTB_GRAPH_FRAMES=f"{d}/frames.tsv")

            def cells(path):
                return [l.rstrip("\n").split("\t") for l in open(path)
                        if l[:1].isdigit()]

            # without the inputs: FATAL, not a silently unchecked ABSENT
            p = subprocess.run(base + ["--out", f"{d}/x.tsv"], capture_output=True,
                               text=True, env=env)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("--h37rv-bam", p.stderr)
            # opted out: ABSENT, labelled as resting on R alone
            p = subprocess.run(base + ["--no-absent-check", "--out", f"{d}/old.tsv"],
                               capture_output=True, text=True, env=env)
            self.assertEqual(p.returncode, 0, p.stderr + p.stdout)
            self.assertEqual(cells(f"{d}/old.tsv"), [["0", "ABSENT", "", "unchecked"]])
            # checked: the reads show the sequence, so the sample's own call
            p = subprocess.run(base + ["--h37rv-bam", bam, "--h37rv-vcf", f"{d}/p1.vcf",
                                       "--panel-diffs", f"{d}/pd.tsv.gz",
                                       "--graph-vcf", g, "--repeat-mask", f"{d}/m.bed",
                                       "--samtools", SAMTOOLS,
                                       "--out", f"{d}/new.tsv"],
                               capture_output=True, text=True, env=env)
            self.assertEqual(p.returncode, 0, p.stderr + p.stdout)
            self.assertEqual(cells(f"{d}/new.tsv"),
                             [["0", "REF", h[2099], "present_here_ref"]])
            self.assertIn("#absent_check\tv1", open(f"{d}/new.tsv").read())
            # and the evidence reaches the states array
            io = load("p5_states_io", "bin/p5_states_io.py")
            os.makedirs(f"{d}/st")
            shutil.copy(f"{d}/new.tsv", f"{d}/st/S.states.tsv")
            kk = list(csv.DictReader(open(f"{d}/keys.tsv"), delimiter="\t"))
            M, EV, codes = io.load_states(f"{d}/st", kk, ["S"], quiet=True,
                                          with_evidence=True)
            self.assertEqual(codes.names[EV[0, 0]], "present_here_ref")
            self.assertEqual(io.NAME[M[0, 0]], "REF")
        finally:
            shutil.rmtree(d, ignore_errors=True)


# ---------------------------------------------------------------------------
class MergedVcfEvidence(unittest.TestCase):
    """ABS-1 item 4 and ABS-2 item 2 in the merged VCF: FORMAT/EV, and an
    inherited insertion written as measured."""

    def test_cells_and_format(self):
        m = load("merge_cohort_vcf_ev", "bin/merge_cohort_vcf.py")
        self.assertEqual(m.vcf_cell("ABSENT", "junction_joined"), "2:ABSENT:junction_joined")
        self.assertEqual(m.vcf_cell("REF"), "0:REF")
        self.assertEqual(m.vcf_cell("NOCALL", "a b;c=d:e,f"), ".:NOCALL:a_b_c_d_e_f")
        line, _ = m.record_line(("c", 5, "k", "A", ["C"], ["CLASS=small"],
                                 "0:REF\t2:ABSENT:masked"))
        self.assertEqual(line.split("\t")[8], "GT:ST:EV")
        line, _ = m.record_line(("c", 5, "k", "A", ["C"], ["CLASS=small"], "0:REF\t1:ALT"))
        self.assertEqual(line.split("\t")[8], "GT:ST")

    def test_inherited_insertion_written_as_measured(self):
        with tempfile.TemporaryDirectory() as d:
            def w(name, text):
                with open(f"{d}/{name}", "w") as fh:
                    fh.write(text)
                return f"{d}/{name}"
            w("m.tsv", "key\tframe\tregion\tkind\th37rv_pos\tnode\tnode_offset\t"
                       "canonical_ref\tA\tB\tC\nh37rv:10:A>C\th37rv\tcore\tSNP\t10\t\t\tA\tALT\tREF\tREF\n")
            svm = w("sv.tsv", "key\tsvtype\th37rv_pos\tsvlen\tsize_band\tn_alt\t"
                              "n_both_callers\tsrc\tmax_qual\tqual_band\tstage11_ppv\t"
                              "qual_from\tcomponent\tA\tB\tC\n"
                              "sv:INS:500:80\tINS\t500\t80\t50-499\t3\t0\tgraph_vcf\t\t"
                              "composed\t0.631\tgraph_vcf\tmixed\tALT\tALT\tALT\n")
            mem = w("mem.tsv", "key\tsample\tmember\tcomponent\n"
                               "sv:INS:500:80\tA\tsv:INS:500:500\tinherited\n"
                               "sv:INS:500:80\tB\tsv:INS:500:500\tinherited\n"
                               "sv:INS:500:80\tC\tsv:INS:502:502\tcalled\n")
            inh = w("inh.tsv", "sample\tkind\tid\tstate\tevidence\n"
                               "A\tINS\tsv:INS:500:500\tNOCALL\ttandem\n"
                               "B\tINS\tsv:INS:500:500\tALT\tmatched_carries\n")
            r = subprocess.run([sys.executable, "bin/merge_cohort_vcf.py",
                                "--matrix", f"{d}/m.tsv", "--sv-matrix", svm,
                                "--inherited-sv-states", inh, "--sv-members", mem,
                                "--cohort-name", "t", "--build-id", "b",
                                "--no-harmonize", "--out", f"{d}/o.vcf"],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            rec = [l.rstrip("\n").split("\t") for l in open(f"{d}/o.vcf")
                   if l.startswith("NC_000962.3\t500\t")][0]
            self.assertEqual(rec[8], "GT:ST:EV")
            # A: inheritance alone, measured tandem; B: carried; C: a caller's call
            self.assertEqual(rec[9:], [".:NOCALL:tandem", "1:ALT:matched_carries", "1:ALT"])
            # the two inputs go together
            r = subprocess.run([sys.executable, "bin/merge_cohort_vcf.py",
                                "--matrix", f"{d}/m.tsv", "--sv-matrix", svm,
                                "--inherited-sv-states", inh,
                                "--cohort-name", "t", "--build-id", "b",
                                "--out", f"{d}/o2.vcf"], capture_output=True, text=True)
            self.assertNotEqual(r.returncode, 0)


class ReferenceOnlySummary(unittest.TestCase):
    """ABS-1 item 4: per sample, calls resting only on the reference, and the
    deeper-sublineage flag."""

    def test_counts_and_flag(self):
        with tempfile.TemporaryDirectory() as d:
            with open(f"{d}/v.vcf", "w") as fh:
                fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\t"
                         "FILTER\tINFO\tFORMAT\tA\tB\n")
                fh.write("c\t1\tk1\tA\tC\t.\tPASS\t.\tGT:ST:EV\t2:ABSENT:node_not_on_r\t"
                         "2:ABSENT:junction_joined\n")
                fh.write("c\t2\tk2\tA\tC\t.\tPASS\t.\tGT:ST:EV\t.:NOCALL:masked\t0:REF\n")
                fh.write("c\t3\tk3\tA\tC\t.\tPASS\t.\tGT:ST\t1:ALT\t0:REF\n")
            with open(f"{d}/refmap.tsv", "w") as fh:
                fh.write("sample\tsample_sublineage\treference\tref_sublineage\n"
                         "A\tlineage1.2.1\tR1\tlineage1.2.1.2.1\n"
                         "B\tlineage4\tR2\tlineage4\n")
            r = subprocess.run([sys.executable, "bin/p5_reference_only.py", "--vcf",
                                f"{d}/v.vcf", "--refmap", f"{d}/refmap.tsv",
                                "--out", f"{d}/o.tsv"], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            got = {x["sample"]: x for x in csv.DictReader(open(f"{d}/o.tsv"), delimiter="\t")}
            self.assertEqual(got["A"]["ref_deeper"], "yes")
            self.assertEqual(got["B"]["ref_deeper"], "no")
            self.assertEqual(got["A"]["reference_only"], "1")
            self.assertEqual(got["A"]["masked"], "1")
            self.assertEqual(got["B"]["absent_confirmed"], "1")


class StaleCommentFixed(unittest.TestCase):
    """ABS-2: p5_sv_matrix.py said inherited rows were excluded by the default
    --keep-filter PASS; the default keeps them."""

    def test_comment(self):
        src = open("bin/p5_sv_matrix.py").read()
        self.assertNotIn("excluded here by the\n    # default --keep-filter PASS", src)
        self.assertNotIn("nothing in this script reads `component` yet", src)
        self.assertIn('ap.add_argument("--keep-filter", default="PASS,INHERITED"', src)


if __name__ == "__main__":
    unittest.main()
