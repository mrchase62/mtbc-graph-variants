#!/usr/bin/env python3
"""Regression tests for the p4_p5 audit findings fixed on 2026-10-05
(analysis/audit/p4_p5.md, and P3IS-2 of p3_is6110_accessory.md).

Same pattern as tests/run_tests.py: unittest, small synthetic inputs, the
standard library plus the pipeline's own scripts, seconds to run.

    $MTB_PY tests/test_audit_p4p5.py -v
"""
import csv
import gzip
import importlib.util
import os
import random
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
    return list(csv.DictReader(open(path, newline=""), delimiter="\t"))


def tsv(path, header, recs):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in recs:
            fh.write("\t".join(str(r.get(c, "")) for c in header) + "\n")
    return path


def sparse_states(path, keys):
    """{key: state} from a sparse states file, REF filled in."""
    got = {}
    for line in open(path):
        if line.startswith("#") or line.startswith("idx"):
            continue
        i, st, _ = line.rstrip("\n").split("\t")
        got[keys[int(i)]["key"]] = st
    return {k["key"]: got.get(k["key"], "REF") for k in keys}


# --------------------------------------------------------------------------
class ProjectionAndReferenceAllele(unittest.TestCase):
    """P4P5-1: a non-zero odgi distance was written ABSENT whatever it meant.
    P4P5-2: a reference block was REF even where R carries the key's ALT.

    Synthetic H37Rv (4 Mb, the size check) and an R that differs from it by a
    SNP at 1000, a 30 bp deletion of 2001-2030 and a 12 bp insertion after
    3000. The projection rows are what odgi gives for each case (measured on
    GCF_014900175): the SNP base lands on R's own base at dist 1, a deleted
    base on the base after the gap at the deleted node's length."""

    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp()
        d = cls.d
        rng = random.Random(7)
        h = "".join(rng.choice("ACGT") for _ in range(4_000_100))
        snp_alt = {"A": "C", "C": "G", "G": "T", "T": "A"}[h[999]]
        ins = "GATTACAGATTC"
        if ins[-1] == h[2999]:
            ins = ins[:-1] + ("A" if h[2999] != "A" else "C")
        comp = str.maketrans("ACGT", "TGCA")
        # ... and 6001-6100 inverted in R, which odgi reports one base off
        r = (h[:999] + snp_alt + h[1000:2000] + h[2030:3000] + ins
             + h[3000:6000] + h[6000:6100].translate(comp)[::-1] + h[6100:])
        cls.h, cls.r, cls.snp_alt, cls.ins = h, r, snp_alt, ins

        def t_of(p):              # R position of H37Rv position p
            return p if p <= 2000 else (p - 30 if p <= 3000 else p - 30 + len(ins))
        os.makedirs(f"{d}/refs")
        for name, s in (("GCF_000195955", h), ("RREF", r)):
            with open(f"{d}/refs/{name}.fasta", "w") as fh:
                fh.write(f">{name}\n")
                for i in range(0, len(s), 80):
                    fh.write(s[i:i + 80] + "\n")
        third = [b for b in "ACGT" if b not in (h[999], snp_alt)][0]
        other = [b for b in "ACGT" if b != h[499]][0]
        cls.keys = [
            # P4P5-1: R carries this SNP; the sample matches R
            dict(key=f"h37rv:1000:{h[999]}>{snp_alt}", h37rv_pos=1000,
                 canonical_ref=h[999], canonical_alt=snp_alt, kind="SNP"),
            # P4P5-1: same site, an allele neither R nor H37Rv has
            dict(key=f"h37rv:1000:{h[999]}>{third}", h37rv_pos=1000,
                 canonical_ref=h[999], canonical_alt=third, kind="SNP"),
            # inside R's deletion: ABSENT is right and must stay
            dict(key=f"h37rv:2015:{h[2014]}>A", h37rv_pos=2015,
                 canonical_ref=h[2014], canonical_alt="A", kind="SNP"),
            # P4P5-2: R carries this insertion; the target is the anchor
            dict(key=f"h37rv:3000:{h[2999]}>{h[2999] + ins}", h37rv_pos=3000,
                 canonical_ref=h[2999], canonical_alt=h[2999] + ins,
                 kind="INDEL"),
            # an ordinary site R shares with H37Rv
            dict(key=f"h37rv:500:{h[499]}>{other}", h37rv_pos=500,
                 canonical_ref=h[499], canonical_alt=other, kind="SNP"),
            # off R's path with an unrelated target: nothing is shown
            dict(key=f"h37rv:5000:{h[4999]}>{other}", h37rv_pos=5000,
                 canonical_ref=h[4999], canonical_alt=other, kind="SNP"),
            # inside R's inversion, R's base is H37Rv's complemented
            dict(key=f"h37rv:6050:{h[6049]}>{[b for b in 'ACGT' if b != h[6049]][0]}",
                 h37rv_pos=6050, canonical_ref=h[6049],
                 canonical_alt=[b for b in "ACGT" if b != h[6049]][0],
                 kind="SNP"),
        ]
        for k in cls.keys:
            k.update(frame="h37rv", node="", node_offset="", region="core",
                     acc_locus="")
        tsv(f"{d}/keys.tsv", ["key", "frame", "h37rv_pos", "node", "node_offset",
                              "canonical_ref", "canonical_alt", "region", "kind",
                              "acc_locus"], cls.keys)
        proj = {1000: (1000, 1), 2015: (t_of(2031), 15), 3000: (t_of(3000), 0),
                500: (500, 0), 5000: (3_500_000, 80)}
        with open(f"{d}/proj.pos", "w") as fh:
            for p, (t, dist) in sorted(proj.items()):
                fh.write(f"H#1#c,{p - 1},+\tRREF#1#c,{t - 1},+\t{dist}\t+\t+\n")
            # the homolog of 6050 is R's t_of(6051); odgi says one past it,
            # with its own strand flag `-`
            fh.write(f"H#1#c,6049,+\tRREF#1#c,{t_of(6051)},+\t0\t+\t-\n")
        tsv(f"{d}/placed.tsv", ["sample", "frame", "h37rv_pos", "ref", "alt",
                                "r_pos", "node", "node_offset"],
            [dict(sample="S", frame="h37rv", h37rv_pos=700, ref=h[699],
                  alt=other if other != h[699] else "T")])
        with open(f"{d}/s.g.vcf", "w") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\t"
                     "FILTER\tINFO\tFORMAT\tS\n")
            fh.write(f"c\t1\t.\t{r[0]}\t<NON_REF>\t.\t.\tEND={len(r)}\tGT:DP\t0:30\n")
        write(f"{d}/frames.tsv", """
            accession\tstrand\toffset\tpanel_len
            RREF\t+\t0\t4000082
            """)
        env = dict(os.environ, MTB_GRAPH_FRAMES=f"{d}/frames.tsv")
        # --ref-fasta is not passed: it defaults to <dir of --h37rv>/<ref>.fasta
        cls.proc = subprocess.run(
            [sys.executable, "bin/p5_states.py", "--sample", "S",
             "--reference", "RREF", "--keys", f"{d}/keys.tsv",
             "--placed", f"{d}/placed.tsv", "--gvcf", f"{d}/s.g.vcf",
             "--projected", f"{d}/proj.pos",
             "--h37rv", f"{d}/refs/GCF_000195955.fasta",
             # the node tables are required now (audit cleanup); this fixture
             # has no node-frame keys, so header-only tables
             "--node-paths", write(f"{d}/node_paths.tsv", "node\tpaths\n"),
             "--node-positions", write(f"{d}/node_positions.tsv",
                                       "node\taccession\tstart\tstrand\t"
                                       "n_occurrences\tlength\n"),
             "--out", f"{d}/s.states.tsv"],
            capture_output=True, text=True, env=env)
        cls.st = (sparse_states(f"{d}/s.states.tsv", cls.keys)
                  if cls.proc.returncode == 0 else {})

    @classmethod
    def tearDownClass(cls):
        import shutil
        shutil.rmtree(cls.d, ignore_errors=True)

    def state(self, i):
        self.assertEqual(self.proc.returncode, 0, self.proc.stderr + self.proc.stdout)
        return self.st[self.keys[i]["key"]]

    def test_r_snp_at_distance_one_is_genotyped_not_absent(self):
        self.assertEqual(self.state(0), "ALT")       # was ABSENT

    def test_third_allele_at_r_snp_is_nocall_not_absent(self):
        self.assertEqual(self.state(1), "NOCALL")    # was ABSENT

    def test_base_r_deletes_stays_absent(self):
        self.assertEqual(self.state(2), "ABSENT")

    def test_unrelated_target_is_nocall_not_absent(self):
        self.assertEqual(self.state(5), "NOCALL")    # was ABSENT

    def test_insertion_r_carries_is_alt_not_ref(self):
        self.assertEqual(self.state(3), "ALT")       # P4P5-2: was REF

    def test_shared_site_stays_ref(self):
        self.assertEqual(self.state(4), "REF")

    def test_locally_inverted_target_read_at_its_homolog(self):
        self.assertEqual(self.state(6), "REF")

    def test_helpers_directly(self):
        p5 = load("p5_states", "bin/p5_states.py")
        h, r = self.h, self.r
        self.assertTrue(p5.homologous(h, r, 1000, 1000, "+"))
        self.assertFalse(p5.homologous(h, r, 2015, 2001, "+"))
        self.assertTrue(p5.deleted_in_ref(h, r, 2015, 2001, "+"))
        self.assertFalse(p5.deleted_in_ref(h, r, 1000, 1000, "+"))
        self.assertEqual(p5.r_allele(h, r, 3000, 2970, "+", h[2999],
                                     h[2999] + self.ins), "ALT")
        self.assertEqual(p5.r_allele(h, r, 500, 500, "+", h[499], "A" if
                                     h[499] != "A" else "C"), "REF")


# --------------------------------------------------------------------------
class TwoFrameConfirmation(unittest.TestCase):
    """P4P5-3: one clip cluster confirmed an inherited deletion although the
    H37Rv frame showed full depth. P4P5-14: depth inside a deletion the
    sample's own gVCF calls was counted as coverage."""

    def genotype(self):
        with tempfile.TemporaryDirectory() as d:
            ivs = [dict(interval="svi:DEL:1001:100", start=1001, end=1100),
                   dict(interval="svi:DEL:3001:100", start=3001, end=3100),
                   dict(interval="svi:DEL:5001:100", start=5001, end=5100),
                   dict(interval="svi:DEL:7001:100", start=7001, end=7100,
                        carriers="")]
            for iv in ivs:
                iv.update(source="graph", svtype="DEL", svlen=100,
                          support_tier="A_multi_assembly", is6110_prox=0)
                c = iv.pop("carriers", "R")
                iv.update(n_ref_carriers=int(bool(c)), ref_carriers=c)
            tsv(f"{d}/iv.tsv", ["interval", "source", "svtype", "start", "end",
                                "svlen", "support_tier", "is6110_prox",
                                "n_ref_carriers", "ref_carriers"], ivs)
            with open(f"{d}/proj.pos", "w") as fh:
                for p in range(1, 9001):
                    fh.write(f"H#1#c,{p - 1},+\tR#1#c,{p - 1},+\t0\t+\t+\n")
            with open(f"{d}/g.vcf", "w") as fh:
                fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
                fh.write("c\t1\t.\tA\t<NON_REF>\t.\t.\tEND=7000\tGT:DP\t0:30\n")
                # the sample's own gVCF calls 7001-7100 deleted, then GATK's
                # next reference block starts inside the deletion
                fh.write("c\t7000\t.\tA" + "C" * 100 + "\tA,<NON_REF>\t500\t.\t.\t"
                         "GT:DP\t1:20\n")
                fh.write("c\t7001\t.\tC\t<NON_REF>\t.\t.\tEND=9000\tGT:DP\t0:20\n")
            tsv(f"{d}/tf.tsv", ["sample", "interval", "h_state", "h_clip_sides"],
                [dict(sample="S", interval="svi:DEL:1001:100", h_state="REF",
                      h_clip_sides=1),
                 dict(sample="S", interval="svi:DEL:3001:100", h_state="REF",
                      h_clip_sides=0),
                 dict(sample="S", interval="svi:DEL:5001:100", h_state="ALT",
                      h_clip_sides=2)])
            r = subprocess.run(
                [sys.executable, "bin/p5_sv_genotype.py", "--intervals",
                 f"{d}/iv.tsv", "--reference", "R", "--sample", "S",
                 "--gvcf", f"{d}/g.vcf", "--projected", f"{d}/proj.pos",
                 "--twoframe", f"{d}/tf.tsv", "--mapq-scope", "never",
                 "--out", f"{d}/out.tsv"], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            return {q["interval"]: q["state"] for q in rows(f"{d}/out.tsv")}

    def test_clip_cluster_does_not_override_full_depth(self):
        st = self.genotype()
        # full H37Rv-frame depth plus one clipped end: conflicting, not ALT
        self.assertEqual(st["svi:DEL:1001:100"], "NOCALL")
        self.assertEqual(st["svi:DEL:3001:100"], "REF")
        self.assertEqual(st["svi:DEL:5001:100"], "ALT")

    def test_called_deletion_is_not_depth(self):
        st = self.genotype()
        # P4P5-14: blocks inside the sample's own called deletion are not depth
        self.assertEqual(st["svi:DEL:7001:100"], "ALT")       # was REF


# --------------------------------------------------------------------------
class CallerGenotype(unittest.TestCase):
    """P4P5-4: delly PASS records genotyped 0/0 became ALT carriers."""

    def test_hom_ref_records_are_not_read(self):
        p4b = load("p4b_place_sv", "bin/p4b_place_sv.py")
        with tempfile.TemporaryDirectory() as d:
            v = write(f"{d}/delly.vcf", """
                ##fileformat=VCFv4.2
                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
                c\t1000\ta\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=1500\tGT:GQ\t0/0:30
                c\t2000\tb\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=2500\tGT:GQ\t0/1:30
                c\t3000\tc\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=3500\tGT:GQ\t1/1:30
                c\t4000\td\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=4500\tGT:GQ\t./.:0
                c\t5000\te\tN\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=5500\tGT\t1
                """)
            got = sorted(r["pos"] for r in p4b.parse_sv_vcf(v, "delly"))
        self.assertEqual(got, [2000, 3000, 5000])


class QualBandCaller(unittest.TestCase):
    """P4P5-11: a two-caller row was banded on delly's scale when the QUAL
    kept was dysgu's."""

    def test_band_on_the_caller_whose_qual_it_is(self):
        with tempfile.TemporaryDirectory() as d:
            tsv(f"{d}/refmap.tsv", ["sample", "reference"],
                [dict(sample="S", reference="R")])
            tsv(f"{d}/S.sv_placed.tsv",
                ["sample", "svtype", "svlen", "frame", "key", "h37rv_pos",
                 "h37rv_end", "src", "n_callers", "sr", "pe", "qual", "filter",
                 "component", "qual_caller"],
                [dict(sample="S", svtype="DEL", svlen=500, frame="h37rv",
                      key="sv:DEL:1000:1500", h37rv_pos=1000, h37rv_end=1500,
                      src="delly,dysgu", n_callers=2, sr=3, pe=3, qual=25,
                      filter="PASS", component="called", qual_caller="dysgu")])
            write(f"{d}/cal.tsv", """
                caller\tbin_lo\tppv
                delly\t100\t0.9
                delly\t0\t0.03
                dysgu\t20\t0.9
                dysgu\t0\t0.5
                """)
            r = subprocess.run([sys.executable, "bin/p5_sv_matrix.py",
                                "--refmap", f"{d}/refmap.tsv", "--dir", d,
                                "--calibration", f"{d}/cal.tsv",
                                "--out", f"{d}/m.tsv"],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            got = rows(f"{d}/m.tsv")
        self.assertEqual([q["qual_band"] for q in got], ["high"])   # was medium


# --------------------------------------------------------------------------
class IntervalCatalogue(unittest.TestCase):
    """P4P5-5: clusters widened to the union of their members and compared
    the next deletion with the widened span."""

    def test_no_widening(self):
        with tempfile.TemporaryDirectory() as d:
            def rec(pos, n, gts):
                return (f"c\t{pos}\t.\t{'A' * (n + 1)}\tA\t.\t.\tLV=0\tGT\t"
                        + "\t".join(gts) + "\n")
            with open(f"{d}/g.vcf", "w") as fh:
                fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL"
                         "\tFILTER\tINFO\tFORMAT\tR1\tR2\tR3\n")
                fh.write(rec(999, 100, ["1", "0", "0"]))     # 1000-1099
                fh.write(rec(1009, 120, ["0", "1", "0"]))    # same event as it
                fh.write(rec(1099, 180, ["0", "0", "1"]))    # 1.8x the first
                # within the position tolerance but sharing no base
                fh.write(rec(1999, 58, ["1", "0", "0"]))
                fh.write(rec(2149, 58, ["0", "1", "0"]))
            r = subprocess.run([sys.executable, "bin/sv_intervals.py",
                                "--graph-vcf", f"{d}/g.vcf", "--is6110-gff", "",
                                "--out", f"{d}/iv.tsv"],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            got = {q["interval"]: q["ref_carriers"] for q in rows(f"{d}/iv.tsv")}
        self.assertEqual(got, {"svi:DEL:1000:100": "R1,R2",
                               "svi:DEL:1100:180": "R3",
                               "svi:DEL:2000:58": "R1",
                               "svi:DEL:2150:58": "R2"})


# --------------------------------------------------------------------------
class InheritedOverlap(unittest.TestCase):
    """P4P5-8: an inherited record was dropped only when a called record sat
    at the very same position, not when one overlapped it."""

    def test_overlapping_inherited_dropped(self):
        with tempfile.TemporaryDirectory() as d:
            hdr = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\t"
                   "FILTER\tINFO\tFORMAT\tS\n")
            write(f"{d}/direct.vcf", hdr)
            write(f"{d}/matched.vcf", hdr + "c\t150\t.\tTGGG\tT\t50\t.\t.\tGT\t1\n")
            write(f"{d}/hpos.tsv", "#src\ttgt\tdist\n"
                  "R#1#c,149,+\tH#1#c,149,+\t0\t+\t+\n")
            write(f"{d}/npos.tsv", "#src\tnode\n" "R#1#c,149,+\t77,3,+\n")
            write(f"{d}/mask.bed", "c\t99\t300\tpe_ppe|PPE1|named\n")
            write(f"{d}/loci.tsv", "pos\tlocus_id\n")
            write(f"{d}/graph.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF"
                  "\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tRREF\n"
                  "c\t120\t.\tT\tTTT\t.\t.\tLV=0\tGT\t1\n"
                  "c\t151\t.\tTTG\tGGG\t.\t.\tLV=0\tGT\t1\n"
                  "c\t180\t.\tA\tG\t.\t.\tLV=0\tGT\t1\n")
            r = subprocess.run(
                [sys.executable, "bin/p4_place.py", "--sample", "S",
                 "--reference", "RREF", "--build-id", "b",
                 "--direct", f"{d}/direct.vcf", "--matched", f"{d}/matched.vcf",
                 "--h37rv-pos", f"{d}/hpos.tsv", "--node-pos", f"{d}/npos.tsv",
                 "--mask", f"{d}/mask.bed", "--loci", f"{d}/loci.tsv",
                 "--graph-vcf", f"{d}/graph.vcf", "--out", f"{d}/placed.tsv"],
                capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            inh = sorted(int(q["h37rv_pos"]) for q in rows(f"{d}/placed.tsv")
                         if q["component"] == "inherited")
        self.assertEqual(inh, [120, 180])        # 151 lies inside 150-153


# --------------------------------------------------------------------------
class MergeInsertionsAndPresence(unittest.TestCase):
    """P3IS-2: loci whose sequence H37Rv already has were written absent in
    every sample. P4P5-6: a caller insertion was written beside the accessory
    or IS6110 record for the same event."""

    def test_merge(self):
        with tempfile.TemporaryDirectory() as d:
            samples = ["S1", "S2", "S3"]
            keys = [dict(key="h37rv:100:A>G", frame="h37rv", h37rv_pos=100,
                         node="", node_offset="", canonical_ref="A",
                         canonical_alt="G", region="core", kind="SNP",
                         acc_locus="")]
            tsv(f"{d}/keys.tsv", list(keys[0]), keys)
            import hashlib
            sha = hashlib.sha1(keys[0]["key"].encode()).hexdigest()
            for sm in samples:
                with open(f"{d}/{sm}.states.tsv", "w") as fh:
                    fh.write(f"#format\tsparse-v1\n#sample\t{sm}\n#default\tREF\n"
                             f"#n_keys\t1\n#keys_sha1\t{sha}\nidx\tstate\tallele\n"
                             + ("0\tALT\tG\n" if sm == "S1" else ""))
            tsv(f"{d}/refmap.tsv", ["sample", "reference"],
                [dict(sample=s, reference="R") for s in samples])
            os.makedirs(f"{d}/acc")
            cat = [dict(locus_id="ACC_0010000", pos=10000, klass="reference_gap",
                        novelty="copy_number", h37rv_cov="1.0", graph_len=1358),
                   dict(locus_id="ACC_0020000", pos=20000, klass="polymorphic",
                        novelty="novel", h37rv_cov="0.0", graph_len=2000)]
            tsv(f"{d}/cat.tsv", ["locus_id", "pos", "klass", "novelty",
                                 "h37rv_cov", "graph_len"], cat)
            for sm in samples:
                tsv(f"{d}/acc/{sm}.presence.tsv", ["sample", "locus", "state"],
                    [dict(sample=sm, locus="ACC_0010000", state="ABSENT"),
                     dict(sample=sm, locus="ACC_0020000",
                          state="PRESENT" if sm == "S1" else "ABSENT")])
            sv = [dict(key="sv:INS:20000:20000", svtype="INS", h37rv_pos=20000,
                       svlen=1999, S1="ALT"),
                  dict(key="sv:INS:30005:30005", svtype="INS", h37rv_pos=30005,
                       svlen=1355, S2="ALT"),
                  dict(key="sv:INS:40000:40000", svtype="INS", h37rv_pos=40000,
                       svlen=500, S3="ALT")]
            for r in sv:
                r.update(component="called", qual_band="high")
                for s in samples:
                    r.setdefault(s, "NOCALL")
            tsv(f"{d}/sv.tsv", ["key", "svtype", "h37rv_pos", "svlen",
                                "component", "qual_band"] + samples, sv)
            tsv(f"{d}/is.tsv", ["key", "sample", "state", "frame", "h37rv_pos",
                                "node", "evidence", "site_class"],
                [dict(key="h37rv:30000", sample=s, frame="h37rv",
                      h37rv_pos=30000, state="ALT" if s == "S2" else "REF",
                      evidence="two_sided", site_class="ref_lacking")
                 for s in samples])
            r = subprocess.run(
                [sys.executable, "bin/merge_cohort_vcf.py", "--states-dir", d,
                 "--refmap", f"{d}/refmap.tsv", "--keys", f"{d}/keys.tsv",
                 "--cohort-name", "t", "--build-id", "b", "--ancestral", "",
                 "--sv-matrix", f"{d}/sv.tsv",
                 "--accessory-presence", f"{d}/acc",
                 "--accessory-catalogue", f"{d}/cat.tsv",
                 "--is6110-keys", f"{d}/is.tsv", "--out", f"{d}/o.vcf"],
                capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            ids = {l.split("\t")[2] for l in open(f"{d}/o.vcf")
                   if not l.startswith("#")}
        self.assertNotIn("acc:ACC_0010000", ids)    # P3IS-2: unmeasurable
        self.assertIn("acc:ACC_0020000", ids)
        self.assertNotIn("sv:INS:20000:20000", ids)  # same event as the locus
        self.assertNotIn("sv:INS:30005:30005", ids)  # same event as IS6110
        self.assertIn("sv:INS:40000:40000", ids)
        self.assertIn("h37rv:30000", ids)

    def test_presence_rule(self):
        lp = load("locus_presence", "accessory/bin/locus_presence.py")
        self.assertTrue(lp.read_route_blind(dict(novelty="copy_number",
                                                 h37rv_cov="1.0")))
        self.assertTrue(lp.read_route_blind(dict(novelty="", h37rv_cov="0.95")))
        self.assertFalse(lp.read_route_blind(dict(novelty="novel",
                                                  h37rv_cov="0.1")))


if __name__ == "__main__":
    unittest.main(warnings="ignore")
