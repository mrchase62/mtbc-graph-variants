#!/usr/bin/env python3
"""Regression tests for the second review's findings (analysis/audit/review2/).

    $MTB_PY tests/test_audit_review2.py      # or via tests/run_tests.py
"""
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
PANSN = "GCF_000195955#1#NC_000962.3"


def write(path, text):
    with open(path, "w") as fh:
        fh.write(textwrap.dedent(text).lstrip("\n"))
    return path


def run(args, **kw):
    return subprocess.run(args, capture_output=True, text=True, **kw)


def bcftools():
    b = os.environ.get("MTB_BCFTOOLS", "")
    return b if b and os.access(b, os.X_OK) else ""


def vcf_gz(d, name, contig, samples, body):
    txt = ("##fileformat=VCFv4.2\n"
           f"##contig=<ID={contig},length=10000>\n"
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="x">\n'
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
           + "\t".join(samples) + "\n")
    for line in textwrap.dedent(body).strip().splitlines():
        txt += contig + "\t" + "\t".join(line.split()) + "\n"
    plain = write(os.path.join(d, name), txt)
    gz = plain + ".gz"
    subprocess.run([bcftools(), "view", "-Oz", "-o", gz, plain], check=True)
    subprocess.run([bcftools(), "index", "-t", gz], check=True)
    return gz


def read_fasta(p):
    out, name = {}, None
    for line in open(p):
        if line.startswith(">"):
            name = line[1:].strip()
            out[name] = ""
        else:
            out[name] += line.strip()
    return out


class CombinedAlignment(unittest.TestCase):
    """R2-TREES-1 / TP-5: the chain's cohort + panel alignment, with the fixed
    outgroup rule applied to every panel genome."""

    def setUp(self):
        if not bcftools():
            self.skipTest("MTB_BCFTOOLS not set; source config/project_env.sh")

    def build(self, d):
        # cohort: isolates S1 S2 and the reference row, two columns
        write(f"{d}/c.fasta", ">S1\nAT\n>S2\nGT\n>H37Rv\nGC\n")
        write(f"{d}/c.sites.tsv", """
            column\tchrom\tpos\tref\talt
            0\tNC_000962.3\t100\tG\tA
            1\tNC_000962.3\t150\tC\tT
            """)
        # the cohort's merged VCF: records at 100, 150 and an indel at 300
        cvcf = vcf_gz(d, "m.vcf", "NC_000962.3", ("S1", "S2"), """
            100 . G A . . . GT 1 0
            150 . C T . . . GT 0 0
            300 . A AT . . . GT 1 0
            """)
        # the panel: O (outgroup), P1, P2, and the panel's own H37Rv copy
        pvcf = vcf_gz(d, "p.vcf", PANSN, ("O", "P1", "P2", "GCF_000195955"), """
            100 . G A . . . GT 1 . 0 0
            148 . ACGT A . . . GT 0 0 1 0
            150 . C T . . . GT 0 1 . 0
            200 . C T . . . GT 0 1 0 0
            300 . A G . . . GT 0 0 1 0
            """)
        r = run([sys.executable, "assoc/bin/combined_alignment.py",
                 "--cohort-alignment", f"{d}/c.fasta",
                 "--cohort-sites", f"{d}/c.sites.tsv", "--vcf", cvcf,
                 "--panel-vcf", pvcf, "--outgroup", "O",
                 "--max-missing", "1.0",
                 "--out", f"{d}/x.fasta", "--sites-out", f"{d}/x.sites.tsv"],
                env=dict(os.environ, MTB_BCFTOOLS=bcftools()))
        return r

    def test_rows_columns_and_rules(self):
        with tempfile.TemporaryDirectory() as d:
            r = self.build(d)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            seq = read_fasta(f"{d}/x.fasta")
            sites = [l.split("\t") for l in open(f"{d}/x.sites.tsv")][1:]
        # the panel's H37Rv copy is left out; the reference row stands for it
        self.assertEqual(list(seq), ["S1", "S2", "H37Rv", "O", "P1", "P2"])
        # columns: the two cohort columns, then the panel-only SNP at 200.
        # 300 is left out: a cohort record covers it, so the isolates' state
        # there is not known
        self.assertEqual([(s[2], s[5].strip()) for s in sites],
                         [("100", "cohort"), ("150", "cohort"), ("200", "panel")])
        # 100: O ALT; P1 missing -> N (not REF); P2 REF
        # 150: P2's 148 ACGT>A deletion covers it -> N; P1 ALT
        # 200: isolates REF (no cohort record), P1 the only carrier
        self.assertEqual(seq["S1"], "ATC")
        self.assertEqual(seq["H37Rv"], "GCC")
        self.assertEqual(seq["O"], "ACC")
        self.assertEqual(seq["P1"], "NTT")
        self.assertEqual(seq["P2"], "GNC")

    def test_outgroup_must_be_a_panel_genome(self):
        with tempfile.TemporaryDirectory() as d:
            self.build(d)
            r = run([sys.executable, "assoc/bin/combined_alignment.py",
                     "--cohort-alignment", f"{d}/c.fasta",
                     "--cohort-sites", f"{d}/c.sites.tsv",
                     "--vcf", f"{d}/m.vcf.gz", "--panel-vcf", f"{d}/p.vcf.gz",
                     "--outgroup", "NOPE",
                     "--out", f"{d}/y.fasta", "--sites-out", f"{d}/y.tsv"],
                    env=dict(os.environ, MTB_BCFTOOLS=bcftools()))
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not a sample", r.stderr)


class PruneForCohort(unittest.TestCase):
    """R2-TREES-1: the combined tree is cut to the cohort, the reference tip
    and the outgroup, and the root stays on the outgroup."""

    def run_prune(self, d, tree, samples=("S1", "S2")):
        write(f"{d}/t.nwk", tree + "\n")
        vcf = write(f"{d}/m.vcf", "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                    + "\t".join(samples) + "\n")
        return run([sys.executable, "assoc/bin/prune_for_cohort.py",
                    "--tree", f"{d}/t.nwk", "--vcf", vcf, "--outgroup", "O",
                    "--out", f"{d}/o.nwk"])

    def test_keeps_cohort_reference_and_outgroup(self):
        try:
            import dendropy, pysam  # noqa: F401
        except ImportError:
            self.skipTest("needs dendropy and pysam ($MTB_PY)")
        with tempfile.TemporaryDirectory() as d:
            r = self.run_prune(
                d, "(O:1,((S1:1,P1:1):1,((S2:1,H37Rv:1):1,P2:1):1):1);")
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            t = dendropy.Tree.get(path=f"{d}/o.nwk", schema="newick",
                                  rooting="force-rooted")
        self.assertEqual(sorted(n.taxon.label for n in t.leaf_node_iter()),
                         ["H37Rv", "O", "S1", "S2"])
        kids = t.seed_node.child_nodes()
        self.assertEqual(len(kids), 2)
        self.assertIn("O", [k.taxon.label for k in kids if k.is_leaf()])

    def test_sample_missing_from_tree_refused(self):
        try:
            import dendropy, pysam  # noqa: F401
        except ImportError:
            self.skipTest("needs dendropy and pysam ($MTB_PY)")
        with tempfile.TemporaryDirectory() as d:
            r = self.run_prune(d, "(O:1,(S1:1,H37Rv:1):1);", samples=("S1", "S9"))
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not tips", r.stderr)


class ChainUsesTheCombinedTree(unittest.TestCase):
    """R2-TREES-1: the chain no longer builds a cohort-only tree."""

    def test_chain_steps(self):
        txt = open("assoc/bin/cohort_assoc_tail.sh").read()
        self.assertIn("combined_alignment.py", txt)
        self.assertIn("prune_for_cohort.py", txt)
        self.assertIn('build_snp_tree.sh" "data/trees/${C}.combined.fasta"', txt)
        self.assertIn('--outgroup-fasta "data/trees/${C}.combined.fasta"', txt)
        self.assertNotIn("og.fasta", txt)


def load(name, rel):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, rel))
    mod = importlib.util.module_from_spec(spec)
    sys.path.insert(0, os.path.join(ROOT, "bin"))
    spec.loader.exec_module(mod)
    return mod


class DeletedInRefPresenceGuard(unittest.TestCase):
    """R2-GENO-1: ABSENT only where R really lacks the position. A base whose
    H37Rv context occurs in R near the target is present, whatever the anchors
    near the target say."""

    @classmethod
    def setUpClass(cls):
        import random
        rnd = random.Random(7)
        cls.h37 = "".join(rnd.choice("ACGT") for _ in range(6000))
        cls.p5 = load("p5_states_r2", "bin/p5_states.py")

    def test_real_deletion_is_absent(self):
        # R lacks H37Rv 1001..1100 (1-based); odgi lands p near the junction
        r = self.h37[:1000] + self.h37[1100:]
        self.assertTrue(self.p5.deleted_in_ref(self.h37, r, 1050, 1000, "+",
                                               dist=50))

    def test_context_present_near_target_is_not_absent(self):
        # the same deletion, but R also carries p's H37Rv context (12 bases
        # each side) 800 bases away -- the base is in R, in another copy
        ctx = self.h37[1049 - 12:1049 + 13]
        r = self.h37[:1000] + self.h37[1100:1800] + ctx + self.h37[1800:]
        self.assertTrue(self.p5._present_in(r[0:3000], self.h37, 1050))
        self.assertFalse(self.p5.deleted_in_ref(self.h37, r, 1050, 1000, "+",
                                                dist=50))

    def test_context_with_a_substitution_at_p_still_present(self):
        ctx = self.h37[1049 - 12:1049] + ("A" if self.h37[1049] != "A" else "C") \
            + self.h37[1050:1050 + 12]
        r = self.h37[:1000] + self.h37[1100:1800] + ctx + self.h37[1800:]
        self.assertFalse(self.p5.deleted_in_ref(self.h37, r, 1050, 1000, "+",
                                                dist=50))


class AccessoryBlindByIdentity(unittest.TestCase):
    """R2-IS-1: a locus is read-route blind only if H37Rv carries 0.9 of its
    sequence at 95% identity or more (h37rv_cov95), not at any identity."""

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, os.path.join(ROOT, "accessory", "bin"))
        import locus_presence
        cls.blind = staticmethod(locus_presence.read_route_blind)

    def test_rule(self):
        old = dict(novelty="copy_number", h37rv_cov="1.0")
        self.assertTrue(self.blind(old))                   # no column: old rule
        self.assertFalse(self.blind(dict(old, h37rv_cov95="0.0")))
        self.assertFalse(self.blind(dict(old, h37rv_cov95="0.706")))
        self.assertTrue(self.blind(dict(old, h37rv_cov95="0.95")))
        self.assertFalse(self.blind(dict(novelty="novel", h37rv_cov="0.0",
                                         h37rv_cov95="")))

    def test_merge_catalogues_measures_identity(self):
        blastn = os.path.join(os.environ.get("MTB_QC_BIN", ""), "blastn")
        if not os.access(blastn, os.X_OK):
            self.skipTest("no blastn at $MTB_QC_BIN; source config/project_env.sh")
        import random
        rnd = random.Random(11)
        h37 = "".join(rnd.choice("ACGT") for _ in range(5000))
        same = h37[3000:3400]
        div = "".join(c if rnd.random() > 0.15 else
                      rnd.choice([x for x in "ACGT" if x != c]) for c in same)
        with tempfile.TemporaryDirectory() as d:
            write(f"{d}/h.fa", f">NC_000962.3\n{h37}\n")
            hdr = ("locus_id\tpos\tklass\trep_len\tn_alleles\tcarriers_any"
                   "\tcarrier_frac\th37rv_cov\tnovelty\n")
            with open(f"{d}/loci.tsv", "w") as fh:
                fh.write(hdr)
                fh.write("ACC_0000500\t500\tpolymorphic\t400\t1\t1\t0.5\t1.0\tcopy_number\n")
                fh.write("ACC_0001500\t1500\tpolymorphic\t400\t1\t1\t0.5\t1.0\tcopy_number\n")
            with open(f"{d}/g.vcf", "w") as fh:
                fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL"
                         "\tFILTER\tINFO\tFORMAT\tA\tB\n")
                fh.write(f"{PANSN}\t500\t.\t{h37[499]}\t{h37[499]}{same}\t.\t.\t.\tGT\t1\t0\n")
                fh.write(f"{PANSN}\t1500\t.\t{h37[1499]}\t{h37[1499]}{div}\t.\t.\t.\tGT\t1\t0\n")
            r = run([sys.executable, "accessory/bin/merge_catalogues.py",
                     "--loci", f"{d}/loci.tsv", "--graph-vcf", f"{d}/g.vcf",
                     "--out", f"{d}/c.tsv", "--out-fasta", f"{d}/c.fa",
                     "--h37rv", f"{d}/h.fa", "--blastn", blastn])
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            import csv
            rows = {x["locus_id"]: x for x in
                    csv.DictReader(open(f"{d}/c.tsv"), delimiter="\t")}
        self.assertGreaterEqual(float(rows["ACC_0000500"]["h37rv_cov95"]), 0.9)
        self.assertLess(float(rows["ACC_0001500"]["h37rv_cov95"]), 0.9)
        self.assertTrue(self.blind(rows["ACC_0000500"]))
        self.assertFalse(self.blind(rows["ACC_0001500"]))


class Level2UnmeasuredLocus(unittest.TestCase):
    """R2-INT-2: a variant inside an accessory locus no sample was measured at
    is left unconditional (it used to have every leaf blanked and leave the
    scan); a locus that was measured and has no carrier still makes its
    variants inapplicable."""

    HDR = ("##fileformat=VCFv4.2\n"
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
           "S1\tS2\tS3\tS4\n")
    TREE = "(O:1,((S1:1,S2:1)x:1,(S3:1,S4:1)y:1)ing:1);\n"

    def test_unmeasured_locus_is_not_conditioned_on(self):
        py = os.environ.get("MTB_PY_VT", "")
        if not (py and os.access(py, os.X_OK) and run(
                [py, "-c", "import mtbvartools, dendropy"]).returncode == 0):
            self.skipTest("MTB_PY_VT with mtbvartools not set; source "
                          "config/project_env.sh")
        with tempfile.TemporaryDirectory() as d:
            write(f"{d}/m.vcf", self.HDR
                  + "NC_000962.3\t100\th37rv:100:A>G\tA\tG\t.\tPASS\t"
                    "CLASS=small;FRAME=h37rv\tGT\t1\t1\t0\t0\n"
                  + "node_7\t1\tnode:7:0:A>G\tA\tG\t.\tPASS\t"
                    "CLASS=small;FRAME=node\tGT\t1\t0\t0\t0\n"
                  + "node_8\t1\tnode:8:0:C>T\tC\tT\t.\tPASS\t"
                    "CLASS=small;FRAME=node\tGT\t1\t0\t0\t0\n")
            write(f"{d}/t.nwk", self.TREE)
            write(f"{d}/og.fa", ">O\nA\n")
            write(f"{d}/og.tsv", "column\tchrom\tpos\tref\talt\n"
                  "0\tNC_000962.3\t100\tA\tG\n")
            write(f"{d}/nl.tsv", "node\tlocus\n7\tLBLIND\n8\tLMEAS\n")
            os.makedirs(f"{d}/pres")
            for smp in ("S1", "S2", "S3", "S4"):
                write(f"{d}/pres/{smp}.presence.tsv",
                      "sample\tlocus\tstate\n"
                      f"{smp}\tLBLIND\tUNMEASURABLE\n{smp}\tLMEAS\tABSENT\n")
            r = run([py, "assoc/bin/write_event_matrix.py",
                     "--vcf", f"{d}/m.vcf", "--tree", f"{d}/t.nwk",
                     "--out", f"{d}/ev", "--outgroup-name", "O",
                     "--outgroup-fasta", f"{d}/og.fa",
                     "--outgroup-sites", f"{d}/og.tsv", "--panel-polarity", "",
                     "--accessory-presence", f"{d}/pres",
                     "--node-locus", f"{d}/nl.tsv"])
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            import csv
            v = {x["id"]: x for x in
                 csv.DictReader(open(f"{d}/ev/variants.tsv"), delimiter="\t")}
        blind = v["node:7:0:A>G"]
        self.assertEqual(blind["acc_locus"], "")
        self.assertEqual(blind["acc_locus_unmeasured"], "LBLIND")
        # S1 alone carries it: one gain on S1's branch, not every leaf unknown
        self.assertEqual((blind["n_derived_leaves"], blind["n_gain"]), ("1", "1"))
        meas = v["node:8:0:C>T"]
        self.assertEqual(meas["acc_locus"], "LMEAS")
        self.assertEqual(meas["n_applicable"], "0")
        self.assertEqual((meas["n_derived_leaves"], meas["n_gain"]), ("0", "0"))


class SameDeletionRule(unittest.TestCase):
    """R2-GENO-2: one rule -- near in position and length AND half the bases
    shared -- for clustering graph deletions, adding caller deletions to the
    catalogue, and the merge's "the catalogue supersedes this caller row"."""

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, os.path.join(ROOT, "bin"))
        import sv_intervals
        cls.S = sv_intervals

    def test_rule(self):
        S = self.S
        # two 58 bp deletions 150 bp apart: same_event alone said yes
        self.assertTrue(S.same_event(1000, 58, 1150, 58))
        self.assertFalse(S.same_deletion(1000, 58, 1150, 58))
        self.assertTrue(S.same_deletion(1000, 58, 1020, 58))     # 38 of 58 shared
        self.assertFalse(S.same_deletion(1000, 58, 1040, 58))    # 18 of 58

    def test_caller_deletion_beside_an_interval_is_added(self):
        h = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER"
             "\tINFO\tFORMAT\tA\tB\n")
        with tempfile.TemporaryDirectory() as d:
            # graph: a 58 bp deletion of bases 1001-1058 (anchor 1000)
            write(f"{d}/g.vcf", h + f"{PANSN}\t1000\t.\t" + "A" * 59
                  + "\tA\t.\t.\t.\tGT\t1\t0\n")
            # caller: one 58 bp deletion 150 bp along (no shared base) and one
            # that is the graph's own (anchor 1000)
            write(f"{d}/m.tsv", "key\tsvtype\th37rv_pos\tsvlen\tn_alt\n"
                  "a\tDEL\t1150\t58\t3\nb\tDEL\t1000\t58\t3\n")
            write(f"{d}/is.gff", "")
            r = run([sys.executable, "bin/sv_intervals.py", "--graph-vcf",
                     f"{d}/g.vcf", "--sv-matrix", f"{d}/m.tsv",
                     "--is6110-gff", f"{d}/is.gff", "--out", f"{d}/iv.tsv"])
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            import csv
            iv = list(csv.DictReader(open(f"{d}/iv.tsv"), delimiter="\t"))
        got = sorted((x["source"], x["start"]) for x in iv)
        self.assertEqual(got, [("caller", "1151"), ("graph", "1001")])

    def test_merge_uses_the_rule(self):
        self.assertIn("_ivm.same_deletion(p0 + 1, L, s0, l0)",
                      open("bin/merge_cohort_vcf.py").read())


if __name__ == "__main__":
    unittest.main(warnings="ignore")
