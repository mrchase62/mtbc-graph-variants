#!/usr/bin/env python3
"""Regression tests for the tree and polarity findings of the 2026-10-05 audit
(analysis/audit/tree_polarity.md, graphvcf.md): Faults A and B, GRAPHVCF-2/3,
TP-1, TP-2, TP-3, TP-6, TP-7, TP-8, TP-9 and TP-10.

Same pattern as tests/run_tests.py: small synthetic inputs, the pipeline's own
scripts run as subprocesses, a few seconds in all.

    $MTB_PY tests/test_audit_polarity.py -v

Tests that need bcftools (the outgroup and the polarity table) read
MTB_BCFTOOLS, and the event-writer tests need MTB_PY_VT (mtbvartools and
dendropy); source config/project_env.sh first or they are skipped.
The writer tests take most of the time: importing mtbvartools costs about
ten seconds per run.
"""
import csv
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


def rows_of(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def run(args, **kw):
    return subprocess.run(args, capture_output=True, text=True, **kw)


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        else:
            seqs[name] += line.strip()
    return seqs


def bcftools():
    b = os.environ.get("MTB_BCFTOOLS", "")
    return b if b and os.access(b, os.X_OK) else ""


def panel_vcf(d, body, samples=("O", "X")):
    """bgzipped, indexed panel VCF on the PanSN H37Rv contig."""
    txt = ("##fileformat=VCFv4.2\n"
           f"##contig=<ID={PANSN},length=10000>\n"
           '##INFO=<ID=AC,Number=A,Type=Integer,Description="x">\n'
           '##INFO=<ID=AN,Number=1,Type=Integer,Description="x">\n'
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="x">\n'
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
           + "\t".join(samples) + "\n")
    for line in textwrap.dedent(body).strip().splitlines():
        txt += PANSN + "\t" + "\t".join(line.split()) + "\n"
    plain = write(os.path.join(d, "panel.vcf"), txt)
    gz = plain + ".gz"
    r = run([bcftools(), "view", "-Oz", "-o", gz, plain])
    assert r.returncode == 0, r.stderr
    r = run([bcftools(), "index", "-t", gz])
    assert r.returncode == 0, r.stderr
    return gz


class AddOutgroup(unittest.TestCase):
    """Fault A / GRAPHVCF-1, GRAPHVCF-2 and TP-1: assoc/bin/add_outgroup.py."""

    def setUp(self):
        if not bcftools():
            self.skipTest("MTB_BCFTOOLS not set; source config/project_env.sh")

    def outgroup_column(self, panel_body, sites):
        with tempfile.TemporaryDirectory() as d:
            vcf = panel_vcf(d, panel_body)
            write(f"{d}/s.tsv", "column\tchrom\tpos\tref\talt\n" + "".join(
                f"{i}\t{c}\t{p}\t{r}\t{a}\n"
                for i, (c, p, r, a) in enumerate(sites)))
            write(f"{d}/a.fa", ">S1\n" + "N" * len(sites) + "\n")
            r = run([sys.executable, "assoc/bin/add_outgroup.py",
                     "--alignment", f"{d}/a.fa", "--sites", f"{d}/s.tsv",
                     "--out", f"{d}/o.fa", "--panel-vcf", vcf,
                     "--outgroup", "O"],
                    env=dict(os.environ, MTB_BCFTOOLS=bcftools()))
            self.assertEqual(r.returncode, 0, r.stderr)
            return read_fasta(f"{d}/o.fa")["O"]

    def test_any_exact_duplicate_with_gt1_is_alt(self):
        # Fault A: three copies of 100 G>C, the outgroup's ALT in the first;
        # last-duplicate-wins wrote REF
        got = self.outgroup_column("""
            100 . G C . . . GT 1 0
            100 . G C . . . GT 0 1
            100 . G C . . . GT 0 0
            """, [("NC_000962.3", 100, "G", "C")])
        self.assertEqual(got, "C")

    def test_all_exact_missing_is_n(self):
        got = self.outgroup_column("""
            200 . A T . . . GT . 1
            200 . A T . . . GT . 0
            """, [("NC_000962.3", 200, "A", "T")])
        self.assertEqual(got, "N")

    def test_upstream_mnp_and_deletion(self):
        # an MNP starting at 299 puts T at 301 (the ALT); a deletion starting
        # at 398 removes 399-401; nothing at all covers 500; an exact GT 0 at
        # 600; a node-frame site is off the panel contig
        got = self.outgroup_column("""
            299 . CAG AAT . . . GT 1 0
            398 . CATG C . . . GT 1 0
            600 . C T . . . GT 0 1
            """, [("NC_000962.3", 301, "G", "T"),
                  ("NC_000962.3", 299, "C", "G"),
                  ("NC_000962.3", 400, "T", "C"),
                  ("NC_000962.3", 398, "C", "A"),
                  ("NC_000962.3", 500, "A", "G"),
                  ("NC_000962.3", 600, "C", "T"),
                  ("node_7", 1, "A", "G")])
        # 299 C>G: the MNP puts A there, a third base -> N
        # 398 C>A: the deletion keeps its anchor base C -> REF
        self.assertEqual(got, "TNNCACN")

    def test_missing_covering_record_is_n(self):
        # the outgroup has no call over a deletion spanning the site
        got = self.outgroup_column("""
            700 . CATG C . . . GT . 1
            """, [("NC_000962.3", 702, "T", "A")])
        self.assertEqual(got, "N")


class AlignmentStarAndSiblings(unittest.TestCase):
    """TP-2 and TP-7: bin/vcf_to_alignment.py."""

    HDR = ("##fileformat=VCFv4.2\n##contig=<ID=c,length=1000>\n"
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="x">\n'
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
           "S1\tS2\tS3\tS4\n")

    def align(self, body, *extra):
        with tempfile.TemporaryDirectory() as d:
            txt = self.HDR + "".join("\t".join(l.split()) + "\n" for l in
                                     textwrap.dedent(body).strip().splitlines())
            write(f"{d}/v.vcf", txt)
            r = run([sys.executable, "bin/vcf_to_alignment.py", "--vcf",
                     f"{d}/v.vcf", "--out", f"{d}/a.fa", "--sites-out",
                     f"{d}/s.tsv", "--max-missing", "0.5", *extra])
            self.assertEqual(r.returncode, 0, r.stderr)
            return read_fasta(f"{d}/a.fa"), rows_of(f"{d}/s.tsv")

    def test_star_record_is_a_column(self):
        seqs, sites = self.align("""
            c 10 . T C,* . . . GT 0 1 2 1
            """, "--ref-sample", "R")
        self.assertEqual([(s["pos"], s["ref"], s["alt"]) for s in sites],
                         [("10", "T", "C")])
        self.assertEqual({k: v for k, v in seqs.items()},
                         {"S1": "T", "S2": "C", "S3": "N", "S4": "C", "R": "T"})

    def test_sibling_allele_carrier_is_not_ref(self):
        # split multi-allelic site: S3 carries G, and is 0 in the C record
        seqs, sites = self.align("""
            c 20 . A C . . . GT 0 1 0 1
            c 20 . A G . . . GT 0 0 1 1
            """)
        self.assertEqual(len(sites), 2)
        self.assertEqual(seqs["S3"][0], "N")      # was A
        self.assertEqual(seqs["S2"][1], "N")      # was A
        self.assertEqual(seqs["S1"], "AA")


class AncestralFaultB(unittest.TestCase):
    """Fault B: bin/ancestral_alleles.py's ingroup is the MRCA of every
    non-outgroup leaf, not the root's other child."""

    def aa(self, tree, fasta, sites, *extra):
        with tempfile.TemporaryDirectory() as d:
            write(f"{d}/t.nwk", tree)
            write(f"{d}/a.fa", fasta)
            write(f"{d}/s.tsv", "chrom\tpos\tref\talt\n" + sites)
            r = run([sys.executable, "bin/ancestral_alleles.py", "--tree",
                     f"{d}/t.nwk", "--alignment", f"{d}/a.fa", "--sites",
                     f"{d}/s.tsv", "--out", f"{d}/o.tsv", *extra])
            return r, (rows_of(f"{d}/o.tsv") if r.returncode == 0 else None)

    def test_second_outgroup_inside_roots_child(self):
        # The panel tree's shape, with its leaf names and the default
        # outgroups: GCF_035581225 roots it, and GCF_000253375 sits in the
        # root's other child. The MTBC (A, B, C) is all C at site 1, so its
        # ancestor is C; the old ingroup (GCF_000253375 + MTBC) tied and the
        # root outgroup made it T. Site 2: the MTBC ties {C, T}; the nearer
        # outgroup says C, the further one T.
        r, rows = self.aa("(GCF_035581225,(GCF_000253375,(A,(B,C))));\n",
                          ">GCF_035581225\nTT\n>GCF_000253375\nTC\n"
                          ">A\nCC\n>B\nCT\n>C\nCT\n",
                          "c\t1\tC\tT\nc\t2\tC\tT\n")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual([(x["AA"], x["polarity"]) for x in rows],
                         [("C", "ref_ancestral"), ("C", "ref_ancestral")])

    def test_mrca_containing_an_outgroup_is_fatal(self):
        r, _ = self.aa("((O1,A),(O2,(B,C)));\n",
                       ">O1\nT\n>O2\nT\n>A\nC\n>B\nC\n>C\nC\n",
                       "c\t1\tC\tT\n", "--outgroups", "O2,O1")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("MRCA", r.stdout + r.stderr)


class PanelPolarityTable(unittest.TestCase):
    """GRAPHVCF-3 / TP-8: bin/panel_polarity.py pools duplicate records."""

    def setUp(self):
        if not bcftools():
            self.skipTest("MTB_BCFTOOLS not set; source config/project_env.sh")

    def test_duplicates_pooled(self):
        with tempfile.TemporaryDirectory() as d:
            vcf = panel_vcf(d, """
                100 . G C . . AC=1;AN=4 GT 1 0 0 .
                100 . G C . . AC=1;AN=4 GT 0 1 0 0
                200 . A T . . AC=1;AN=4 GT 0 1 0 0
                """, samples=("O", "X", "Y", "Z"))
            r = run([sys.executable, "bin/panel_polarity.py", "--panel-vcf",
                     vcf, "--outgroup", "O", "--bcftools", bcftools(),
                     "--h37rv", self.h37rv(d, "T" * 10000),
                     "--out", f"{d}/p.tsv"])
            self.assertEqual(r.returncode, 0, r.stderr)
            got = [(x["chrom"], x["pos"], x["ancestral"], x["panel_af"])
                   for x in rows_of(f"{d}/p.tsv")]
        # 100: the outgroup is ALT in one copy (the last copy said REF);
        # carriers O and X of the four genomes, Z called in the second copy
        self.assertEqual(got, [(PANSN, "100", "ALT", "0.5000"),
                               (PANSN, "200", "REF", "0.2500")])

    @staticmethod
    def h37rv(d, seq):
        return write(os.path.join(d, "h37rv.fasta"), f">NC_000962.3\n{seq}\n")

    def test_keys_left_aligned_like_the_cohort_vcf(self):
        """Review 2, R2-INT-1. The collapsed VCF trims but does not left-align;
        the cohort VCF's keys go through mtb_norm.normalise. An insertion in a
        homopolymer written at the run's right end must be keyed where the
        cohort VCF keys it, and two right-shifted copies of the one event must
        pool into a single row (the outgroup is ALT in the second copy)."""
        sys.path.insert(0, os.path.join(ROOT, "bin"))
        import mtb_norm
        # 1-based: 299 C, 300-305 AAAAAA, 306 G, T elsewhere
        seq = list("T" * 1000)
        seq[298] = "C"
        for i in range(299, 305):
            seq[i] = "A"
        seq[305] = "G"
        seq = "".join(seq)
        with tempfile.TemporaryDirectory() as d:
            vcf = panel_vcf(d, """
                302 . A AA . . AC=1;AN=4 GT 0 1 0 0
                305 . A AA . . AC=1;AN=4 GT 1 0 1 0
                """, samples=("O", "X", "Y", "Z"))
            r = run([sys.executable, "bin/panel_polarity.py", "--panel-vcf",
                     vcf, "--outgroup", "O", "--bcftools", bcftools(),
                     "--h37rv", self.h37rv(d, seq), "--out", f"{d}/p.tsv"])
            self.assertEqual(r.returncode, 0, r.stderr)
            got = rows_of(f"{d}/p.tsv")
        cohort_key = mtb_norm.normalise(seq, 305, "A", "AA")
        self.assertEqual(cohort_key, (299, "C", "CA"))
        self.assertEqual(len(got), 1, got)
        x = got[0]
        self.assertEqual((int(x["pos"]), x["ref"], x["alt"]), cohort_key)
        self.assertEqual((x["ancestral"], x["n_records"], x["panel_af"]),
                         ("ALT", "2", "0.7500"))

    def test_h37rv_required(self):
        r = run([sys.executable, "bin/panel_polarity.py", "--panel-vcf", "x",
                 "--outgroup", "O", "--out", "o"])
        self.assertEqual(r.returncode, 2)
        self.assertIn("--h37rv", r.stderr)


_VT = []


def vt_python():
    """MTB_PY_VT if it can import the writer's libraries, checked once: the
    import alone takes about ten seconds on this filesystem."""
    if not _VT:
        p = os.environ.get("MTB_PY_VT", "")
        ok = p and os.access(p, os.X_OK) and \
            run([p, "-c", "import mtbvartools, dendropy"]).returncode == 0
        _VT.append(p if ok else "")
    return _VT[0]


class EventWriter(unittest.TestCase):
    """TP-2, TP-3, TP-6, TP-8 and TP-9: assoc/bin/write_event_matrix.py."""

    VCF_HDR = ("##fileformat=VCFv4.2\n"
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
               "S1\tS2\tS3\tS4\n")
    # S1,S2 and S3,S4 are sister pairs: a site splitting them evenly is tied
    # inside the ingroup, and only the outgroup can say which state is old.
    TREE = "(O:1,((S1:1,S2:1)x:1,(S3:1,S4:1)y:1)ing:1);\n"

    def setUp(self):
        self.py = vt_python()
        if not self.py:
            self.skipTest("MTB_PY_VT with mtbvartools not set; source "
                          "config/project_env.sh")
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.d, ignore_errors=True)

    def writer(self, records, og_sites, og_seq, tree=None, extra=(),
               polarity=""):
        d = self.d
        write(f"{d}/m.vcf", self.VCF_HDR + "".join(
            "\t".join(r.split()) + "\n" for r in records))
        write(f"{d}/t.nwk", tree or self.TREE)
        write(f"{d}/og.fa", f">O\n{og_seq}\n")
        write(f"{d}/og.tsv", "column\tchrom\tpos\tref\talt\n" + "".join(
            f"{i}\t{c}\t{p}\t{r}\t{a}\n"
            for i, (c, p, r, a) in enumerate(og_sites)))
        r = run([self.py, "assoc/bin/write_event_matrix.py",
                 "--vcf", f"{d}/m.vcf", "--tree", f"{d}/t.nwk",
                 "--out", f"{d}/ev", "--outgroup-name", "O",
                 "--outgroup-fasta", f"{d}/og.fa",
                 "--outgroup-sites", f"{d}/og.tsv",
                 "--panel-polarity", polarity, *extra])
        return r, (rows_of(f"{d}/ev/variants.tsv") if r.returncode == 0
                   else None)

    def test_star_record_reads_the_outgroup(self):
        # TP-2: 500 T>C,* ties inside the ingroup; the outgroup is T (REF)
        r, v = self.writer(
            ["NC_000962.3 500 h37rv:500:T>C T C,* . PASS "
             "CLASS=small;FRAME=h37rv GT 1 1 0 2"],
            [("NC_000962.3", 500, "T", "C")], "T")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(v[0]["root_state"], "ancestral")
        self.assertEqual(v[0]["n_gain"], "1")

    def test_star_record_reads_the_polarity_table(self):
        # TP-2 / TP-8: the table row is keyed on the SNP, with its contig
        pol = write(f"{self.d}/pol.tsv",
                    "chrom\tpos\tref\talt\toutgroup_gt\tancestral\tpanel_af\n"
                    f"{PANSN}\t500\tT\tC\t1\tALT\t0.5000\n"
                    f"{PANSN}\t1\tA\tG\t1\tALT\t0.5000\n")
        r, v = self.writer(
            ["NC_000962.3 500 h37rv:500:T>C T C,* . PASS "
             "CLASS=small;FRAME=h37rv GT 1 1 0 2",
             "node_7 1 node:7:0:A>G A G . PASS "
             "CLASS=small;FRAME=node GT 1 1 0 0"],
            [("NC_000962.3", 500, "T", "C")], "T", polarity=pol)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(v[0]["polarity"], "alt_ancestral_outgroup")
        # a node-frame record at POS=1 is not the H37Rv position-1 row
        self.assertEqual(v[1]["polarity"], "unpolarised")

    def test_reference_tip_unknown_at_node_frame(self):
        # TP-3: the H37Rv tip's fabricated REF at a node-frame record is what
        # resolved the root to ancestral, and so made the gain on (S3,S4)
        tree = "(O:1,(H37Rv:1,((S1:1,S2:1)p:1,(S3:1,S4:1)c:1)q:1)ing:1);\n"
        r, v = self.writer(
            ["node_7 1 node:7:0:A>G A G . PASS "
             "CLASS=small;FRAME=node GT 0 0 1 1",
             "NC_000962.3 100 h37rv:100:A>G A G . PASS "
             "CLASS=small;FRAME=h37rv GT 0 0 1 1"],
            [("NC_000962.3", 100, "A", "G")], "N", tree=tree,
            extra=("--ref-sample", "H37Rv"))
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(v[0]["n_gain"], "0")
        # on the H37Rv frame the reference tip is still REF
        self.assertEqual(v[1]["n_gain"], "1")

    def test_dedupe_first_keeps_outgroup_aligned(self):
        # TP-6: 100 is duplicated; 200 ties inside the ingroup and its
        # outgroup allele is ALT, so its root is derived
        rec100 = ("NC_000962.3 100 h37rv:100:A>G A G . PASS "
                  "CLASS=small;FRAME=h37rv GT 1 0 0 0")
        r, v = self.writer(
            [rec100, rec100,
             "NC_000962.3 200 h37rv:200:C>T C T . PASS "
             "CLASS=small;FRAME=h37rv GT 1 1 0 0"],
            [("NC_000962.3", 100, "A", "G"), ("NC_000962.3", 200, "C", "T")],
            "AT", extra=("--dedupe", "first"))
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual([x["pos"] for x in v], ["100", "200"])
        self.assertEqual(v[1]["root_state"], "derived")

    def test_missing_polarity_table_is_fatal(self):
        # TP-9: a table that was asked for and is absent stops the run
        rec = ["NC_000962.3 100 h37rv:100:A>G A G . PASS "
               "CLASS=small;FRAME=h37rv GT 1 0 0 0"]
        site = [("NC_000962.3", 100, "A", "G")]
        r, _ = self.writer(rec, site, "A",
                           polarity=f"{self.d}/no_such_table.tsv")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("panel-polarity", r.stdout + r.stderr)
        r, _ = self.writer(rec, site, "A", polarity="")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("panel polarity    none",
                      open(f"{self.d}/ev/summary.txt").read())


class MergeHeaderCounts(unittest.TestCase):
    """TP-10: the merged VCF's AA header carries no counts from one table."""

    def test_no_hard_coded_counts(self):
        src = open("bin/merge_cohort_vcf.py").read()
        for stale in ("72,986", "6,567", "371 of"):
            self.assertNotIn(stale, src)


if __name__ == "__main__":
    unittest.main(warnings="ignore")
