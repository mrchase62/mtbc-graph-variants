#!/usr/bin/env python3
"""Regression tests for HANDOFF SV-EVAL-2 (2026-10-16): level-1 accessory
presence.

On new1k_p200 the merged VCF held 99 of the build catalogue's 564 loci and
43% of their cells were NOCALL with no reason. Two causes:

  1. a present MOSAIC locus -- part of it sequence H37Rv also has -- was
     measured by breadth over its whole length, but reads from the shared
     part align to H37Rv and never reach the unmapped/clipped pool, so it
     stayed UNCERTAIN. It is now measured over its H37Rv-unique positions.
  2. a locus whose sequence H37Rv carries (read_route_blind) was UNMEASURABLE
     and the merge skipped it. It is now read at its H37Rv anchor: insertion
     junction reads against reads running across the anchor.

The merge writes every catalogue locus, with FORMAT/EV on each NOCALL, and
still reads tables written before (no `route` column).

The fixtures are a synthetic 24 kb "H37Rv" and reads simulated from samples
that do or do not carry the loci, aligned with the pipeline's bwa.

    $MTB_PY tests/test_audit_accessory_presence.py -v
"""
import csv
import importlib.util
import os
import random
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
BWA = os.environ.get("MTB_BWA", "")
try:
    import pysam
except ImportError:                   # pragma: no cover
    pysam = None


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, rel))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


lp = load("locus_presence", "accessory/bin/locus_presence.py")
CTG = "NC_000962.3"


def rnd(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def tsv(path, header, recs):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in recs:
            fh.write("\t".join(str(r.get(c, "")) for c in header) + "\n")
    return path


def fasta(path, seqs):
    with open(path, "w") as fh:
        for n, s in seqs.items():
            fh.write(f">{n}\n")
            for i in range(0, len(s), 80):
                fh.write(s[i:i + 80] + "\n")
    return path


class Fixture:
    """H37Rv G; MOSAIC = 600 bp of G (5,001-5,600) + 600 novel bases;
    NOVEL = 800 novel bases; BLIND = an exact copy of G 15,001-16,360 (an
    'IS6110' H37Rv already carries). The carrier inserts MOSAIC after G base
    9,000 and BLIND after G base 3,000; the non-carrier is G itself."""

    def __init__(self, d):
        rng = random.Random(7)
        self.d = d
        G = rnd(24000, rng)
        self.G = G
        self.loci = {"ACC_0009000": G[5000:5600] + rnd(600, rng),
                     "ACC_0012000": rnd(800, rng),
                     "ACC_0003000": G[15000:16360]}
        self.pos = {"ACC_0009000": 9000, "ACC_0012000": 12000, "ACC_0003000": 3000}
        self.h37 = fasta(f"{d}/h37rv.fa", {CTG: G})
        self.cat_fa = fasta(f"{d}/cat.fasta", self.loci)
        self.cat = tsv(f"{d}/cat.tsv",
                       ["locus_id", "pos", "klass", "graph_len", "novelty",
                        "h37rv_cov95", "panel_carriers"],
                       [dict(locus_id=k, pos=self.pos[k], klass="reference_gap",
                             graph_len=len(s),
                             novelty="copy_number" if k == "ACC_0003000" else "mosaic",
                             h37rv_cov95="1.0" if k == "ACC_0003000" else "0.5",
                             panel_carriers="R1") for k, s in self.loci.items()])
        carrier = (G[:3000] + self.loci["ACC_0003000"] + G[3000:9000]
                   + self.loci["ACC_0009000"] + G[9000:])
        self.bam = {"C": self.reads(carrier, "C"), "N": self.reads(G, "N")}

    def reads(self, genome, tag, rl=150, step=4):
        fq = f"{self.d}/{tag}.fq"
        with open(fq, "w") as fh:
            for i in range(0, len(genome) - rl, step):
                s = genome[i:i + rl]
                fh.write(f"@{tag}{i}\n{s}\n+\n{'I' * rl}\n")
            # every real isolate has some reads H37Rv cannot place; without
            # them the pool is empty and the step (rightly) refuses to run
            jr = random.Random(5)
            for i in range(20):
                fh.write(f"@{tag}junk{i}\n{rnd(rl, jr)}\n+\n{'I' * rl}\n")
        if not os.path.exists(self.h37 + ".bwt"):
            subprocess.run([BWA, "index", self.h37], check=True, capture_output=True)
        out = f"{self.d}/{tag}.bam"
        sam = subprocess.run([BWA, "mem", self.h37, fq], check=True,
                             capture_output=True, text=True).stdout
        subprocess.run([SAMTOOLS, "sort", "-o", out, "-"], input=sam, text=True,
                       check=True, capture_output=True)
        subprocess.run([SAMTOOLS, "index", out], check=True)
        return out

    def call(self, tag, extra=()):
        out = f"{self.d}/{tag}.presence.tsv"
        r = subprocess.run([sys.executable, "accessory/bin/locus_presence.py",
                            "--sample", tag, "--reference", "R1",
                            "--catalogue", self.cat, "--fasta", self.cat_fa,
                            "--cram", self.bam[tag], "--h37rv", self.h37,
                            "--index-dir", f"{self.d}/idx",
                            "--samtools", SAMTOOLS, "--bwa", BWA, "--out", out,
                            *extra], capture_output=True, text=True)
        return r, ({x["locus"]: x for x in csv.DictReader(open(out), delimiter="\t")}
                   if r.returncode == 0 else {})


def need_tools(tc):
    if pysam is None:
        tc.skipTest("pysam not importable; run with $MTB_PY")
    if not (SAMTOOLS and os.path.exists(SAMTOOLS) and BWA and os.path.exists(BWA)):
        tc.skipTest("MTB_SAMTOOLS / MTB_BWA not set; source config/project_env.sh")


class PresenceCall(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.fx = None
        if pysam is not None and SAMTOOLS and BWA and os.path.exists(BWA):
            cls.fx = Fixture(cls.tmp.name)
            cls.rc, cls.c = cls.fx.call("C")
            cls.rn, cls.n = cls.fx.call("N")

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def setUp(self):
        need_tools(self)
        self.assertEqual(self.rc.returncode, 0, self.rc.stderr)
        self.assertEqual(self.rn.returncode, 0, self.rn.stderr)

    def test_unique_positions_exclude_the_h37rv_half(self):
        idx = [os.path.join(r, "unique.tsv") for r, _, f in os.walk(f"{self.fx.d}/idx")
               if "unique.tsv" in f]
        self.assertEqual(len(idx), 1)
        u = {r["locus"]: int(r["n_unique"]) for r in csv.DictReader(open(idx[0]), delimiter="\t")}
        self.assertLess(u["ACC_0009000"], 700)      # the novel half and its edge
        self.assertGreater(u["ACC_0009000"], 500)
        self.assertEqual(u["ACC_0012000"], 800)
        self.assertLess(u["ACC_0003000"], lp.MIN_UNIQUE)

    def test_mosaic_present_by_unique_breadth(self):
        r = self.c["ACC_0009000"]
        # whole-locus breadth is about half: the old rule's UNCERTAIN
        self.assertLess(float(r["cov_unfiltered"]), 0.8)
        self.assertEqual((r["route"], r["state"], r["evidence_code"]),
                         ("pool", "PRESENT", "unique_breadth"))
        self.assertEqual(self.n["ACC_0009000"]["state"], "ABSENT")

    def test_novel_absent(self):
        for t in (self.c, self.n):
            self.assertEqual((t["ACC_0012000"]["route"], t["ACC_0012000"]["state"]),
                             ("pool", "ABSENT"))

    def test_blind_locus_read_at_its_anchor(self):
        c, n = self.c["ACC_0003000"], self.n["ACC_0003000"]
        self.assertEqual((c["route"], c["state"], c["evidence_code"]),
                         ("junction", "PRESENT", "junction_reads"))
        self.assertEqual((n["route"], n["state"], n["evidence_code"]),
                         ("junction", "ABSENT", "spanning_reads"))
        self.assertEqual(c["presence_version"], lp.PRESENCE_VERSION)

    def test_the_p1_alignment_is_read_when_given(self):
        r, t = self.fx.call("C", ["--h37rv-bam", self.fx.bam["N"]])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(t["ACC_0003000"]["state"], "ABSENT")   # N's reads there
        self.assertEqual(t["ACC_0009000"]["state"], "PRESENT")  # pool from --cram

    def test_missing_sequence_is_fatal(self):
        fasta(self.fx.cat_fa + ".short.fa", {"ACC_0009000": self.fx.loci["ACC_0009000"]})
        r = subprocess.run([sys.executable, "accessory/bin/locus_presence.py",
                            "--sample", "C", "--catalogue", self.fx.cat,
                            "--fasta", self.fx.cat_fa + ".short.fa",
                            "--cram", self.fx.bam["C"], "--h37rv", self.fx.h37,
                            "--samtools", SAMTOOLS, "--bwa", BWA,
                            "--out", f"{self.fx.d}/x.tsv"], capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("no sequence", r.stderr)


class JunctionRule(unittest.TestCase):
    """junction_call on hand-made reads: clips must carry the locus's bases."""

    def test_clips_of_other_sequence_do_not_count(self):
        need_tools(self)
        rng = random.Random(3)
        G = rnd(4000, rng)
        S = rnd(500, rng)
        other = rnd(500, rng)
        with tempfile.TemporaryDirectory() as d:
            hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CTG, "LN": len(G)}]}
            path = f"{d}/a.bam"
            with pysam.AlignmentFile(path + ".u.bam", "wb", header=hdr) as out:
                for i in range(10):     # right-clipped at 2000, clip = `other`
                    a = pysam.AlignedSegment()
                    a.query_name = f"r{i}"; a.reference_id = 0
                    a.reference_start = 2000 - 100 + i
                    a.mapping_quality = 60
                    a.cigartuples = [(0, 100 - i), (4, 50 + i)]
                    a.query_sequence = G[a.reference_start:2000] + other[:50 + i]
                    a.query_qualities = pysam.qualitystring_to_array("I" * 150)
                    out.write(a)
            pysam.sort("-o", path, path + ".u.bam"); pysam.index(path)
            bam = pysam.AlignmentFile(path)
            st, code, nj, nr = lp.junction_call(bam, CTG, 2000, S, G)
            self.assertEqual(nj, 0)
            st2, _, nj2, _ = lp.junction_call(bam, CTG, 2000, other, G)
        self.assertEqual((st, code), ("UNCERTAIN", "junction_other_sequence"))
        self.assertEqual((st2, nj2), ("PRESENT", 10))


class MergeWritesEveryLocus(unittest.TestCase):
    """bin/merge_cohort_vcf.py: every catalogue locus is a record; NOCALL
    cells carry EV; an older table (no `route`) still merges, its blind loci
    NOCALL read_route_blind; IS6110 copies are flagged ACCIS6110."""

    def test_merge(self):
        rng = random.Random(11)
        with tempfile.TemporaryDirectory() as d:
            samples = ["S1", "S2", "S3"]
            keys = [dict(key="h37rv:100:A>G", frame="h37rv", h37rv_pos=100,
                         node="", node_offset="", canonical_ref="A",
                         canonical_alt="G", region="core", kind="SNP", acc_locus="")]
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
            G = rnd(60000, rng)
            fasta(f"{d}/h37.fa", {CTG: G})
            with open(f"{d}/is.gff", "w") as fh:
                fh.write(f"{CTG}\tRefSeq\tmobile_genetic_element\t40001\t41355\t.\t+\t.\t"
                         f"Note=IS6110-1\n")
            loci = {"ACC_0010000": G[40000:41355],          # an IS6110 copy, blind
                    "ACC_0020000": rnd(2000, rng),           # measurable
                    "ACC_0030000": rnd(900, rng)}            # no table has a row
            fasta(f"{d}/cat.fasta", loci)
            tsv(f"{d}/cat.tsv", ["locus_id", "pos", "klass", "novelty",
                                 "h37rv_cov95", "graph_len"],
                [dict(locus_id="ACC_0010000", pos=10000, klass="reference_gap",
                      novelty="copy_number", h37rv_cov95="1.0", graph_len=1355),
                 dict(locus_id="ACC_0020000", pos=20000, klass="polymorphic",
                      novelty="novel", h37rv_cov95="0.0", graph_len=2000),
                 dict(locus_id="ACC_0030000", pos=30000, klass="polymorphic",
                      novelty="novel", h37rv_cov95="0.0", graph_len=900)])
            os.makedirs(f"{d}/acc")
            v2 = ["sample", "locus", "state", "evidence_code", "route"]
            # S1, S2: version 2; S3: a table from before
            tsv(f"{d}/acc/S1.presence.tsv", v2,
                [dict(sample="S1", locus="ACC_0010000", state="PRESENT",
                      evidence_code="junction_reads", route="junction"),
                 dict(sample="S1", locus="ACC_0020000", state="PRESENT",
                      evidence_code="unique_breadth", route="pool")])
            tsv(f"{d}/acc/S2.presence.tsv", v2,
                [dict(sample="S2", locus="ACC_0010000", state="UNCERTAIN",
                      evidence_code="junction_low_depth", route="junction"),
                 dict(sample="S2", locus="ACC_0020000", state="UNCERTAIN",
                      evidence_code="partial_unique_coverage", route="pool")])
            tsv(f"{d}/acc/S3.presence.tsv", ["sample", "locus", "state"],
                [dict(sample="S3", locus="ACC_0010000", state="ABSENT"),
                 dict(sample="S3", locus="ACC_0020000", state="UNCERTAIN")])
            r = subprocess.run(
                [sys.executable, "bin/merge_cohort_vcf.py", "--states-dir", d,
                 "--refmap", f"{d}/refmap.tsv", "--keys", f"{d}/keys.tsv",
                 "--cohort-name", "t", "--build-id", "b", "--ancestral", "",
                 "--h37rv-fasta", f"{d}/h37.fa", "--h37rv-length", str(len(G)),
                 "--is6110-gff", f"{d}/is.gff",
                 "--accessory-presence", f"{d}/acc",
                 "--accessory-catalogue", f"{d}/cat.tsv", "--out", f"{d}/o.vcf"],
                capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            recs = {f[2]: f for f in (l.rstrip("\n").split("\t") for l in open(f"{d}/o.vcf")
                                      if not l.startswith("#"))}
        self.assertEqual({k for k in recs if k.startswith("acc:")},
                         {"acc:ACC_0010000", "acc:ACC_0020000", "acc:ACC_0030000"})
        self.assertEqual(recs["acc:ACC_0010000"][9:],
                         ["1:ALT", ".:NOCALL:junction_low_depth", ".:NOCALL:read_route_blind"])
        self.assertEqual(recs["acc:ACC_0020000"][9:],
                         ["1:ALT", ".:NOCALL:partial_unique_coverage", ".:NOCALL:uncertain"])
        self.assertEqual(set(recs["acc:ACC_0030000"][9:]), {".:NOCALL:no_presence_row"})
        self.assertIn("ACCIS6110", recs["acc:ACC_0010000"][7].split(";"))
        self.assertNotIn("ACCIS6110", recs["acc:ACC_0020000"][7].split(";"))
        self.assertEqual(recs["acc:ACC_0010000"][8], "GT:ST:EV")


if __name__ == "__main__":
    unittest.main(warnings="ignore")
