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


if __name__ == "__main__":
    unittest.main(warnings="ignore")
