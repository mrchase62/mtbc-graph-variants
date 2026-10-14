#!/usr/bin/env python3
"""Tests for bin/panel_ancestral.py, the one Fitch reconstruction over the
panel tree that replaced ancestral_alleles.py and panel_polarity.py (user,
2026-10-13: "yes, use parsimony over all outgroups for every variant class"),
and for P0 step ancestral, which runs it.

Same pattern as tests/run_tests.py: small synthetic inputs, the pipeline's own
scripts run as subprocesses. The three tests the decision asked for:

  TatcHelyElement      the tatC-helY element is ancestral-present under every
                       outgroup order (and, with the kit's SV rows still
                       missing, tied under every order -- never absent);
  FaultB               the ancestral node is the MRCA of the ingroup, never a
                       node containing an outgroup;
  TieOnlyOutgroupsSolve a tie inside lineages 1-4 that only the outgroups
                       resolve.

Then format and compatibility with the tables' readers, the kit, the alignment
helper and the P0 step. Needs MTB_BCFTOOLS (source config/project_env.sh) and
numpy; the event-writer check needs MTB_PY_VT.

    $MTB_PY tests/test_audit_parsimony.py -v
"""
import csv
import hashlib
import itertools
import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest
import warnings

import numpy as np

warnings.simplefilter("ignore", ResourceWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
PANSN = "GCF_000195955#1#NC_000962.3"
H37 = "GCF_000195955"
# the element: 68 bp inserted after H37Rv 2,352,065 (the graph VCF's ALT)
ELEMENT_POS = 2352065
ELEMENT = "C" + "CTCGCCTGGGCTGGCGAGCAGACGCAAAATCCCCCGCACGCCCGGCGTGTCGGGGGATTTTGCGTCTG"
L51, L61, L8, ET = "GCF_022870265", "GCF_000253355", "GCF_012923765", "GCF_035581225"


def write(path, text):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write(textwrap.dedent(text).lstrip("\n"))
    return path


def rows_of(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def bcftools():
    b = os.environ.get("MTB_BCFTOOLS") or shutil.which("bcftools") or ""
    if not b:
        c = "/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/mtb_isolates/bin/bcftools"
        b = c if os.path.exists(c) else ""
    return b if b and os.access(b, os.X_OK) else ""


def env(**kw):
    e = dict(os.environ)
    e.update(MTB_SITE_FILE="/dev/null", MTB_PY=sys.executable)
    for k in ("MTB_BUILD_DIR", "OG", "GRAPH_DIR", "BUILD_ID", "ANC_TREE",
              "ANC_ALN", "ANC_SITES", "ANC_OUTGROUPS", "ANC_KIT", "ANC_KIT_EVENTS", "OUTGROUP",
              "MTB_OUTGROUP", "BUILD_ROOT", "ACCESSORY_DIR"):
        e.pop(k, None)
    if bcftools():
        e["MTB_BCFTOOLS"] = bcftools()
    e.update({k: str(v) for k, v in kw.items()})
    return e


def run(cmd, cwd=None, **kw):
    return subprocess.run(cmd, cwd=cwd, env=env(**kw), capture_output=True,
                          text=True)


def panel_vcf(path, samples, body, length):
    """bgzipped, indexed panel VCF on the PanSN H37Rv contig. `body` lines:
    POS REF ALT GT... (one GT per sample)."""
    txt = ("##fileformat=VCFv4.2\n"
           f"##contig=<ID={PANSN},length={length}>\n"
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="x">\n'
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
           + "\t".join(samples) + "\n")
    for line in textwrap.dedent(body).strip().splitlines():
        f = line.split()
        txt += "\t".join([PANSN, f[0], ".", f[1], f[2], ".", ".", ".", "GT"]
                         + f[3:]) + "\n"
    plain = write(path, txt)
    r = run([bcftools(), "view", "-Oz", "-o", plain + ".gz", plain])
    assert r.returncode == 0, r.stderr
    r = run([bcftools(), "index", "-t", "-f", plain + ".gz"])
    assert r.returncode == 0, r.stderr
    return plain + ".gz"


def h37rv(path, seq):
    return write(path, f">NC_000962.3\n{seq}\n")


def states(path, seq, overrides):
    """A genome_states.py npz: H37Rv's bases with `overrides` {pos: base}."""
    code = {"A": 1, "C": 2, "G": 3, "T": 4, "N": 0}
    s = np.frombuffer(seq.encode(), dtype=np.uint8).copy()
    lut = np.zeros(256, dtype=np.uint8)
    for b, c in code.items():
        lut[ord(b)] = c
    s = lut[s]
    for p, b in overrides.items():
        s[p - 1] = code[b]
    np.savez(path, state=s)
    return path


def ancestral(d, tree, vcf, ref, outgroups, kit=(), extra=()):
    write(f"{d}/t.nwk", tree)
    cmd = [sys.executable, "bin/panel_ancestral.py", "--tree", f"{d}/t.nwk",
           "--panel-vcf", vcf, "--h37rv", ref, "--ref-sample", H37,
           "--outgroups", ",".join(outgroups), "--bcftools", bcftools(),
           "--out-ancestral", f"{d}/anc.tsv", "--out-polarity", f"{d}/pol.tsv",
           *extra]
    for acc, p in kit:
        cmd += ["--kit", f"{acc}={p}"]
    r = run(cmd)
    if r.returncode != 0:
        return r, None, None
    return r, rows_of(f"{d}/anc.tsv"), rows_of(f"{d}/pol.tsv")


def by_pos(rows):
    return {(int(x["pos"]), x["ref"], x["alt"]): x for x in rows}


class _Base(unittest.TestCase):
    def setUp(self):
        if not bcftools():
            self.skipTest("MTB_BCFTOOLS not set; source config/project_env.sh")
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)


# ==========================================================================
# 1. The tatC-helY element
# ==========================================================================
class TatcHelyElement(_Base):
    """The 68 bp element at H37Rv 2,352,065. The shape of the P200 panel tree:
    canettii roots it, then lineage 8, then L5.1 and L6.1 joining the lineage
    1-4,7 clade at ONE node. L5.1 carries the element, L6.1 does not, and
    inside lineages 1-4,7 the carriers and non-carriers tie (one clade each;
    H37Rv, a non-carrier, among them). Nearest-first let the outgroup list's
    order decide this; Fitch over all outgroups must not."""

    TREE = (f"({ET},({L8},(({L51},{L61}),((A1,A2),((B1,B2),{H37})))));\n")
    OUTGROUPS = [L51, L61, L8, ET]

    def seq(self):
        s = ["T"] * 2352200
        s[ELEMENT_POS - 1] = "C"
        s[2352100 - 1] = "A"
        return "".join(s)

    def inputs(self, kit_in_graph):
        seq = self.seq()
        ref = h37rv(f"{self.d}/h37rv.fasta", seq)
        if kit_in_graph:
            # canettii and L8 in the graph, carrying the element (CX333; or
            # the kit once OUT-1 gives it event rows)
            samples = ["A1", "A2", "B1", "B2", L51, L61, L8, ET]
            body = (f"{ELEMENT_POS} C {ELEMENT} 1 1 0 0 1 0 1 1\n"
                    "2352100 A G 1 1 0 0 0 0 0 0\n")
            kit = []
        else:
            # the P200 design: L8 and canettii are kit rows, SNP states only
            samples = ["A1", "A2", "B1", "B2", L51, L61]
            body = (f"{ELEMENT_POS} C {ELEMENT} 1 1 0 0 1 0\n"
                    "2352100 A G 1 1 0 0 0 0\n")
            kit = [(L8, states(f"{self.d}/{L8}.npz", seq, {})),
                   (ET, states(f"{self.d}/{ET}.npz", seq, {}))]
        vcf = panel_vcf(f"{self.d}/panel.vcf", samples, body, len(seq))
        return vcf, ref, kit

    def run_orders(self, kit_in_graph):
        vcf, ref, kit = self.inputs(kit_in_graph)
        outs = set()
        states_ = set()
        for order in itertools.permutations(self.OUTGROUPS):
            r, anc, pol = ancestral(self.d, self.TREE, vcf, ref, list(order), kit)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            outs.add(open(f"{self.d}/pol.tsv").read() + open(f"{self.d}/anc.tsv").read())
            p = by_pos(pol).get((ELEMENT_POS, "C", ELEMENT))
            states_.add(p["ancestral"] if p else "TIED")
        return outs, states_

    def test_present_under_every_outgroup_order(self):
        outs, st = self.run_orders(kit_in_graph=True)
        self.assertEqual(len(outs), 1)              # byte-identical, all 24 orders
        self.assertEqual(st, {"ALT"})               # the element is ancestral

    def test_kit_without_event_rows_ties_never_absent(self):
        # until OUT-1, the kit puts no constraint on an SV; L5.1 and L6.1
        # disagree, so the element is TIED (left unpolarised) under every
        # order -- not "absent", which nearest-first gave with L6.1 first
        outs, st = self.run_orders(kit_in_graph=False)
        self.assertEqual(len(outs), 1)
        self.assertEqual(st, {"TIED"})


# ==========================================================================
# 2. Audit fault B
# ==========================================================================
class FaultB(_Base):
    """The ancestral node is the MRCA of every non-outgroup leaf. CX333's
    tree has two canettii, the second inside the root's other child, and
    taking "the root's child that is not the outgroup" reported the
    canettii-plus-MTBC node."""

    def test_second_outgroup_inside_roots_child(self):
        # site 100: the MTBC is all C, the two canettii T: the MTBC ancestor
        # is C (the old node, canettii2 + MTBC, tied and the root made it T).
        # Site 200: the MTBC ties {C, T}; canettii2 says C, canettii1 T.
        seq = "A" * 300
        ref = h37rv(f"{self.d}/h.fasta", seq[:99] + "C" + seq[100:199] + "C" + seq[200:])
        vcf = panel_vcf(f"{self.d}/p.vcf", ["GCF_035581225", "GCF_000253375", "A", "B", "C"],
                        "100 C T 1 1 0 0 0\n200 C T 1 0 0 1 1\n", 300)
        tree = f"(GCF_035581225,(GCF_000253375,((A,{H37}),(B,C))));\n"
        r, anc, pol = ancestral(self.d, tree, vcf, ref, ["GCF_035581225", "GCF_000253375"])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual([(x["AA"], x["polarity"]) for x in anc],
                         [("C", "ref_ancestral"), ("C", "ref_ancestral")])

    def test_mrca_containing_an_outgroup_is_fatal(self):
        seq = "A" * 300
        ref = h37rv(f"{self.d}/h.fasta", seq)
        vcf = panel_vcf(f"{self.d}/p.vcf", ["O1", "O2", "A", "B", "C"],
                        "100 A T 1 1 0 0 0\n", 300)
        r, _, _ = ancestral(self.d, f"((O1,A),(O2,(B,(C,{H37}))));\n", vcf, ref, ["O2", "O1"])
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("MRCA", r.stdout + r.stderr)

    def test_outgroup_left_off_the_list_is_fatal(self):
        # P200: lineage 8 omitted from the outgroups makes the "ancestral
        # node" the one above L5.1 and L6.1
        seq = "A" * 300
        ref = h37rv(f"{self.d}/h.fasta", seq)
        vcf = panel_vcf(f"{self.d}/p.vcf", ["A", "B", L51, L61, L8, ET],
                        "100 A T 1 0 0 0 0 0\n", 300)
        tree = f"({ET},({L8},(({L51},{L61}),(A,(B,{H37})))));\n"
        r, _, _ = ancestral(self.d, tree, vcf, ref, [L51, L61, ET])
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("fault B", r.stdout + r.stderr)


# ==========================================================================
# 3. A tie inside lineages 1-4 that only the outgroups resolve
# ==========================================================================
class TieOnlyOutgroupsSolve(_Base):
    """Two ingroup clades disagree, so Fitch inside lineages 1-4,7 alone
    ties; the outgroups decide, together."""

    TREE = f"({ET},({L8},(({L51},{L61}),((A1,A2),((B1,B2),{H37})))));\n"
    OG = [L51, L61, L8, ET]

    def test_outgroups_resolve(self):
        seq = list("A" * 2000)
        for p, b in ((100, "C"), (200, "G"), (300, "T"), (400, "A")):
            seq[p - 1] = b
        seq = "".join(seq)
        ref = h37rv(f"{self.d}/h.fasta", seq)
        # A clade carries ALT, B clade + H37Rv REF everywhere
        body = ("100 C T 1 1 0 0 1 1\n"     # L5.1 and L6.1 both ALT -> ALT
                "200 G A 1 1 0 0 0 0\n"     # both REF -> REF
                "300 T C 1 1 0 0 1 0\n"     # they disagree; the kit says ALT
                "400 A AGG 1 1 0 0 1 0\n"   # they disagree; kit has no indels
                "500 A G 1 0 0 0 0 0\n")    # decided inside: REF
        vcf = panel_vcf(f"{self.d}/p.vcf", ["A1", "A2", "B1", "B2", L51, L61], body, 2000)
        kit = [(L8, states(f"{self.d}/l8.npz", seq, {300: "C"})),
               (ET, states(f"{self.d}/et.npz", seq, {300: "C"}))]
        r, anc, pol = ancestral(self.d, self.TREE, vcf, ref, self.OG, kit)
        self.assertEqual(r.returncode, 0, r.stderr)
        p = by_pos(pol)
        self.assertEqual((p[(100, "C", "T")]["ancestral"], p[(100, "C", "T")]["method"]),
                         ("ALT", "outgroups"))
        self.assertEqual((p[(200, "G", "A")]["ancestral"], p[(200, "G", "A")]["method"]),
                         ("REF", "outgroups"))
        self.assertEqual((p[(300, "T", "C")]["ancestral"], p[(300, "T", "C")]["method"]),
                         ("ALT", "outgroups"))
        self.assertNotIn((400, "A", "AGG"), p)         # tied: left unpolarised
        self.assertEqual((p[(500, "A", "G")]["ancestral"], p[(500, "A", "G")]["method"]),
                         ("REF", "ingroup"))
        a = by_pos(anc)
        self.assertEqual(a[(300, "T", "C")]["AA"], "C")
        self.assertEqual(a[(300, "T", "C")]["polarity"], "alt_ancestral")
        self.assertNotIn((400, "A", "AGG"), a)         # ancestral.tsv: SNPs only

    def test_without_the_kit_the_snp_ties(self):
        seq = "A" * 299 + "T" + "A" * 1700
        ref = h37rv(f"{self.d}/h.fasta", seq)
        vcf = panel_vcf(f"{self.d}/p.vcf", ["A1", "A2", "B1", "B2", L51, L61, L8, ET],
                        "300 T C 1 1 0 0 1 0 . .\n", 2000)
        r, anc, pol = ancestral(self.d, self.TREE, vcf, ref, self.OG)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(pol, [])
        # the A clade and L5.1 carry C apart: 2 changes, so HOMOPLASTIC as
        # well (AA-1: the flags combine)
        self.assertEqual((anc[0]["AA"], anc[0]["flag"], anc[0]["polarity"]),
                         (".", "TIED:CT,HOMOPLASTIC", "unknown"))


# ==========================================================================
# 4. Formats, the readers, pooling, normalisation, the kit, refusals
# ==========================================================================
class TableFormats(_Base):

    def test_columns_the_readers_need(self):
        seq = "A" * 99 + "C" + "A" * 900
        ref = h37rv(f"{self.d}/h.fasta", seq)
        vcf = panel_vcf(f"{self.d}/p.vcf", ["O", "X", "Y"],
                        "100 C T 1 1 0\n100 C CA 0 1 0\n", 1000)
        r, anc, pol = ancestral(self.d, f"(O,(X,(Y,{H37})));\n", vcf, ref, ["O"])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(open(f"{self.d}/anc.tsv").readline().rstrip("\n").split("\t"),
                         ["chrom", "pos", "ref", "alt", "AA", "flag", "polarity",
                          "changes", "method"])
        self.assertEqual(open(f"{self.d}/pol.tsv").readline().rstrip("\n").split("\t"),
                         ["chrom", "pos", "ref", "alt", "outgroup_gt", "ancestral",
                          "panel_af", "n_records", "method", "changes", "flag"])
        self.assertEqual([(x["chrom"], x["pos"], x["ref"], x["alt"]) for x in anc],
                         [(PANSN, "100", "C", "T")])
        # merge_cohort_vcf.py's join: (pos, ref, alt) -> (AA, flag, changes).
        # The ingroup (X,(Y,H37)) ties and O decides: OUTGROUP (AA-1), 1 change
        j = {(int(x["pos"]), x["ref"], x["alt"]): (x["AA"], x["flag"], x["changes"])
             for x in anc}
        self.assertEqual(j[(100, "C", "T")], ("T", "OUTGROUP", "1"))
        x = by_pos(pol)[(100, "C", "T")]
        self.assertEqual((x["outgroup_gt"], x["ancestral"], x["panel_af"], x["n_records"]),
                         ("1", "ALT", "0.6667", "1"))

    def test_duplicates_pooled(self):
        # GRAPHVCF-3 / TP-8: copies of one key are pooled (the last copy said
        # REF for O); panel_af counts each genome once
        vcf = panel_vcf(f"{self.d}/p.vcf", ["O", "X", "Y", "Z"],
                        "100 G C 1 0 0 .\n100 G C 0 1 0 0\n200 A T 0 1 0 0\n", 10000)
        seq = "T" * 99 + "G" + "T" * 99 + "A" + "T" * 9800
        ref = h37rv(f"{self.d}/h.fasta", seq)
        r, _, pol = ancestral(self.d, f"(O,(X,(Y,(Z,{H37}))));\n", vcf, ref, ["O"])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual([(x["chrom"], x["pos"], x["ancestral"], x["panel_af"], x["n_records"])
                          for x in pol],
                         [(PANSN, "100", "ALT", "0.5000", "2"),
                          (PANSN, "200", "REF", "0.2500", "1")])

    def test_keys_left_aligned_like_the_cohort_vcf(self):
        # review 2, R2-INT-1: two right-shifted copies of one insertion in a
        # homopolymer pool into the cohort VCF's key
        sys.path.insert(0, os.path.join(ROOT, "bin"))
        import mtb_norm
        seq = list("T" * 1000)
        seq[298] = "C"
        for i in range(299, 305):
            seq[i] = "A"
        seq[305] = "G"
        seq = "".join(seq)
        ref = h37rv(f"{self.d}/h.fasta", seq)
        vcf = panel_vcf(f"{self.d}/p.vcf", ["O", "X", "Y", "Z"],
                        "302 A AA 0 1 0 0\n305 A AA 1 0 1 0\n", 1000)
        r, _, pol = ancestral(self.d, f"(O,(X,(Y,(Z,{H37}))));\n", vcf, ref, ["O"])
        self.assertEqual(r.returncode, 0, r.stderr)
        key = mtb_norm.normalise(seq, 305, "A", "AA")
        self.assertEqual(key, (299, "C", "CA"))
        self.assertEqual(len(pol), 1, pol)
        x = pol[0]
        self.assertEqual((int(x["pos"]), x["ref"], x["alt"]), key)
        self.assertEqual((x["ancestral"], x["n_records"], x["panel_af"]),
                         ("ALT", "2", "0.7500"))

    def test_sibling_allele_is_no_call(self):
        # TP-7 / R2-BUILD-3: O is GT 0 at A>G because it carries A>C, a
        # sibling record changing the same base. Read as REF, it would make
        # the tied ingroup REF-ancestral; it is no call, so the key ties
        seq = "T" * 99 + "A" + "T" * 900
        ref = h37rv(f"{self.d}/h.fasta", seq)
        vcf = panel_vcf(f"{self.d}/p.vcf", ["O", "X", "Y"],
                        "100 A G 0 1 0\n100 A C 1 0 0\n", 1000)
        r, anc, pol = ancestral(self.d, f"(O,(X,(Y,{H37})));\n", vcf, ref, ["O"])
        self.assertEqual(r.returncode, 0, r.stderr)
        a = by_pos(anc)
        self.assertEqual(a[(100, "A", "G")]["flag"], "TIED:AG")
        self.assertEqual(a[(100, "A", "C")]["AA"], "A")

    def test_ancestral_inside_an_ingroup_tie_from_the_ingroup_alone(self):
        # review 3.7 (was run_tests.AncestralAtIngroup): the root ties, the
        # ingroup does not -- the ingroup's state is reported
        seq = "A" * 99 + "C" + "A" * 900
        ref = h37rv(f"{self.d}/h.fasta", seq)
        vcf = panel_vcf(f"{self.d}/p.vcf", ["O", "A", "B", "C"],
                        "100 C T 1 0 0 0\n", 1000)
        r, anc, pol = ancestral(self.d, f"(O,(A,(B,(C,{H37}))));\n", vcf, ref, ["O"])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual((anc[0]["AA"], pol[0]["method"]), ("C", "ingroup"))


class Refusals(_Base):

    def setUp(self):
        super().setUp()
        self.ref = h37rv(f"{self.d}/h.fasta", "A" * 1000)
        self.vcf = panel_vcf(f"{self.d}/p.vcf", ["O", "X", "Y"], "100 A T 1 1 0\n", 1000)

    def test_required_arguments(self):
        r = run([sys.executable, "bin/panel_ancestral.py", "--tree", "t",
                 "--panel-vcf", "v", "--out-ancestral", "a", "--out-polarity", "p"])
        self.assertEqual(r.returncode, 2)
        for n in ("--h37rv", "--ref-sample", "--outgroups"):
            self.assertIn(n, r.stderr)

    def test_tree_and_vcf_must_name_the_same_genomes(self):
        r, _, _ = ancestral(self.d, f"(O,(X,(Q,{H37})));\n", self.vcf, self.ref, ["O"])
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("exactly the VCF's samples", r.stderr)

    def test_kit_must_be_an_outgroup(self):
        k = states(f"{self.d}/k.npz", "A" * 1000, {})
        r, _, _ = ancestral(self.d, f"(O,(K,(X,(Y,{H37}))));\n", self.vcf, self.ref,
                            ["O"], kit=[("K", k)])
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not in --outgroups", r.stderr)

    def test_kit_of_another_length_refused(self):
        k = states(f"{self.d}/k.npz", "A" * 999, {})
        r, _, _ = ancestral(self.d, f"(K,(O,(X,(Y,{H37}))));\n", self.vcf, self.ref,
                            ["O", "K"], kit=[("K", k)])
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("H37Rv's length", r.stderr)

    def test_accessions_checked(self):
        acc = write(f"{self.d}/acc.txt", f"O\nX\n{H37}\n")
        r, _, _ = ancestral(self.d, f"(O,(X,(Y,{H37})));\n", self.vcf, self.ref, ["O"],
                            extra=["--accessions", acc])
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("accessions", r.stderr)


# ==========================================================================
# 5. The alignment helper and the tree-taxa check
# ==========================================================================
class KitRows(unittest.TestCase):

    def test_rows_from_states(self):
        with tempfile.TemporaryDirectory() as d:
            seq = "ACGTACGTAC"
            aln = write(f"{d}/a.fa", ">X\nGTA\n>Y\nCTA\n")
            sites = write(f"{d}/s.tsv", "column\tchrom\tpos\tref\talt\n"
                          f"0\t{PANSN}\t2\tC\tG\n1\t{PANSN}\t4\tT\tC\n"
                          f"2\t{PANSN}\t5\tA\tG\n")
            # K: pos 2 G (ALT), pos 4 A (a third base -> N), pos 5 uncalled
            k = states(f"{d}/k.npz", seq, {2: "G", 4: "A", 5: "N"})
            r = run([sys.executable, "bin/append_kit_rows.py", "--alignment", aln,
                     "--sites", sites, "--kit", f"K={k}", "--out", f"{d}/o.fa"])
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(open(f"{d}/o.fa").read(), ">X\nGTA\n>Y\nCTA\n>K\nGNN\n")
            r = run([sys.executable, "bin/append_kit_rows.py", "--alignment", aln,
                     "--sites", sites, "--kit", f"X={k}", "--out", f"{d}/o.fa"])
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("already has a row", r.stderr)

    def test_tree_taxa(self):
        with tempfile.TemporaryDirectory() as d:
            t = write(f"{d}/t.nwk", f"(ET,(L8,(A,(B,{H37}))));\n")
            acc = write(f"{d}/acc.txt", f"A\nB\n{H37}\n")
            ok = run([sys.executable, "bin/p0_check.py", "tree-taxa", "--tree", t,
                      "--accessions", acc, "--kit", "ET,L8"])
            self.assertEqual(ok.returncode, 0, ok.stderr)
            bad = run([sys.executable, "bin/p0_check.py", "tree-taxa", "--tree", t,
                       "--accessions", acc, "--kit", "ET"])
            self.assertNotEqual(bad.returncode, 0)
            self.assertIn("taxa of", bad.stderr)


# ==========================================================================
# 6. P0 step ancestral
# ==========================================================================
class P0StepAncestral(_Base):
    """One step writes both tables; the build's stamped outgroup (P200:
    GCF_022870265, an in-graph L5.1) must be one of the tree's in-graph
    outgroups; the old panel_polarity step refuses with a pointer."""

    def setUp(self):
        super().setUp()
        d = self.d
        os.symlink(os.path.join(ROOT, "bin"), f"{d}/bin")
        self.og = write(f"{d}/g/g.smooth.final.og", "graph bytes\n")
        gsha = hashlib.sha256(b"graph bytes\n").hexdigest()
        self.b = f"{d}/build/{gsha[:12]}"
        write(f"{self.b}/build_info.tsv",
              f"build_id\t{gsha[:12]}\ngraph_sha256\t{gsha}\noutgroup\t{L51}\n")
        samples = ["A1", "A2", L51, L61]
        write(f"{self.b}/assets/accessions.txt", "".join(s + "\n" for s in samples + [H37]))
        seq = "A" * 99 + "C" + "A" * 900
        h37rv(f"{self.b}/refs/{H37}.fasta", seq)
        vcf = panel_vcf(f"{d}/v/panel.vcf", samples,
                        "100 C T 1 0 1 1\n100 C CGG 1 0 1 0\n", 1000)
        shutil.copy(vcf, f"{self.b}/assets/graph_collapsed.vcf.gz")
        for s in ("accessions", "assets", "refs"):
            write(f"{self.b}/logs/{s}.done", "2026-10-14\n")
        self.tree = write(f"{d}/tree.nwk",
                          f"({ET},({L8},(({L51},{L61}),(A1,(A2,{H37})))));\n")
        self.kit = (f"{L8}={states(f'{d}/l8.npz', seq, {100: 'T'})},"
                    f"{ET}={states(f'{d}/et.npz', seq, {100: 'T'})}")

    def p0(self, *args, **kw):
        kw.setdefault("BUILD_ROOT", f"{self.d}/build")
        kw.setdefault("MTB_GATK_SIF", "/nonexistent.sif")
        return run(["bash", os.path.join(ROOT, "bin/p0_prepare.sh"), *args],
                   cwd=self.d, OG=self.og, **kw)

    def step(self, **kw):
        kw.setdefault("ANC_TREE", self.tree)
        kw.setdefault("ANC_OUTGROUPS", f"{L51},{L61},{L8},{ET}")
        kw.setdefault("ANC_KIT", self.kit)
        kw.setdefault("OUTGROUP", L51)
        return self.p0("--step", "ancestral", **kw)

    def test_writes_both_tables_and_records_inputs(self):
        r = self.step()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        anc = rows_of(f"{self.b}/assets/ancestral.tsv")
        pol = by_pos(rows_of(f"{self.b}/assets/panel_polarity.tsv"))
        self.assertEqual([(x["AA"], x["polarity"]) for x in anc], [("T", "alt_ancestral")])
        self.assertEqual(pol[(100, "C", "T")]["ancestral"], "ALT")
        self.assertNotIn((100, "C", "CGG"), pol)       # L5.1/L6.1 disagree: tied
        done = open(f"{self.b}/logs/ancestral.done").read()
        for k in ("tree_sha256", "graph_vcf_sha256", "h37rv_sha256", "code_sha256",
                  f"kit_{L8}_sha256", f"kit_{ET}_sha256", "outgroups\t", f"outgroup\t{L51}"):
            self.assertIn(k, done)
        self.assertTrue(os.path.exists(f"{self.b}/assets/ancestral_inputs/kit/{ET}.npz"))
        # done: a second run is a no-op; a changed kit is refused
        self.assertEqual(self.step().returncode, 0)
        states(f"{self.d}/l8.npz", "A" * 1000, {})
        r = self.step()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn(f"kit_{L8}_sha256", r.stderr)

    def test_stamped_outgroup_must_match(self):
        r = self.step(OUTGROUP=L61)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn(f"stamped with outgroup {L51}", r.stderr)

    def test_outgroup_must_be_an_in_graph_tree_outgroup(self):
        r = self.step(ANC_OUTGROUPS=f"{L61},{L8},{ET}", ANC_TREE=write(
            f"{self.d}/t2.nwk", f"({ET},({L8},({L61},({L51},(A1,(A2,{H37}))))));\n"))
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not one of ANC_OUTGROUPS", r.stderr)

    def test_taxa_checked(self):
        r = self.step(ANC_KIT=f"{L8}={self.d}/l8.npz")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("taxa of", r.stderr)

    def test_needs_tree_and_outgroups(self):
        r = self.p0("--step", "ancestral", OUTGROUP=L51)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("needs ANC_TREE and ANC_OUTGROUPS", r.stderr)

    def test_panel_polarity_step_points_to_ancestral(self):
        r = self.p0("--step", "panel_polarity", OUTGROUP=L51)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("part of step 'ancestral'", r.stderr)

    def kit_events(self):
        # route A files for both kit genomes: the indel 100 C>CGG is carried
        # by L8 and lacked (callable, nothing there) by ET1291
        write(f"{self.d}/ev/{L8}.var", "R\tNC_000962.3\t0\t1000\n"
              "V\tNC_000962.3\t100\t100\t1\t60\t-\tgg\tq\t1\t3\t+\n")
        write(f"{self.d}/ev/{ET}.var", "R\tNC_000962.3\t0\t1000\n")
        for g in (L8, ET):
            write(f"{self.d}/ev/{g}.sv.tsv",
                  "type\th37rv_start\th37rv_end\tref_len\talt_len\tsource\n")
        return f"{L8}={self.d}/ev/{L8},{ET}={self.d}/ev/{ET}"

    def test_kit_events_used_and_recorded(self):
        r = self.step(ANC_KIT_EVENTS=self.kit_events())
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        pol = by_pos(rows_of(f"{self.b}/assets/panel_polarity.tsv"))
        # L5.1 1, L6.1 0, the ingroup ties; L8 1, ET1291 0 -> the root ties,
        # but the L8 node below it resolves: ALT
        self.assertEqual(pol[(100, "C", "CGG")]["ancestral"], "ALT")
        done = open(f"{self.b}/logs/ancestral.done").read()
        for k in (f"kit_{L8}_var_sha256", f"kit_{ET}_sv_sha256"):
            self.assertIn(k, done)
        self.assertTrue(os.path.exists(f"{self.b}/assets/ancestral_inputs/kit/{L8}.var"))

    def test_kit_events_genome_must_be_a_kit_genome(self):
        r = self.step(ANC_KIT_EVENTS=f"{L51}={self.d}/ev/x")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not in ANC_KIT", r.stderr)


# ==========================================================================
# 6b. Kit events (route A, option A of 2026-10-14)
# ==========================================================================
def var_file(path, R, V):
    """A paftools call file: R = [(s, e)] 0-based half-open, V = [(s, e,
    ref, alt)] in paftools' form ('-' for an empty allele, lower case)."""
    txt = "".join(f"R\tNC_000962.3\t{s}\t{e}\n" for s, e in R)
    txt += "".join(f"V\tNC_000962.3\t{s}\t{e}\t1\t60\t{r}\t{a}\tq\t0\t0\t+\n"
                   for s, e, r, a in V)
    return write(path, txt)


def sv_file(path, rows):
    return write(path, "type\th37rv_start\th37rv_end\tref_len\talt_len\tsource\n"
                 + "".join("\t".join(map(str, r)) + "\n" for r in rows))


class KitEventRules(unittest.TestCase):
    """kit_event_states(): exact allele -> 1; callable and untouched -> 0;
    SVs by type, start within 50 bp and length within 20%; ambiguity -> no
    call."""

    def setUp(self):
        sys.path.insert(0, os.path.join(ROOT, "bin"))
        import panel_ancestral
        self.P = panel_ancestral
        self.d = tempfile.mkdtemp()
        # 1-based: 200-205 AAAAAA homopolymer, G elsewhere
        s = ["G"] * 20000
        for i in range(199, 205):
            s[i] = "A"
        self.seq = "".join(s)

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def states(self, keys, R, V, svs=()):
        pre = f"{self.d}/k"
        var_file(pre + ".var", R, V)
        sv_file(pre + ".sv.tsv", svs)
        ev = self.P.read_kit_events(pre, self.seq, "NC_000962.3")
        keys = [(PANSN, *k) for k in keys]
        kinds = [self.P.kind_of(k[2], k[3]) for k in keys]
        return list(self.P.kit_event_states(keys, kinds, self.seq, "NC_000962.3", ev))

    def test_small_indel_rules(self):
        keys = [(199, "GA", "G"),          # 1-bp deletion in the run, normalised
                (500, "G", "GT"),          # callable, nothing there: 0
                (700, "G", "GT"),          # a kit SNP touches it: no call
                (19990, "G", "GT")]        # not callable: no call
        # paftools writes the deletion at the run's right end (0-based 204)
        st = self.states(keys, R=[(0, 15000)],
                         V=[(204, 205, "a", "-"), (700, 701, "g", "c")])
        self.assertEqual(st, [1, 0, -1, -1])

    def test_two_alignments_is_not_callable(self):
        st = self.states([(500, "G", "GT")], R=[(0, 500), (500, 15000)], V=[])
        self.assertEqual(st, [-1])

    def test_mnp_never_alt(self):
        st = self.states([(800, "GG", "TC"), (900, "GG", "TC")], R=[(0, 15000)],
                         V=[(799, 800, "g", "t"), (800, 801, "g", "c")])
        self.assertEqual(st, [-1, 0])

    def test_sv_by_position_and_length(self):
        ins = "G" + "T" * 100
        keys = [(3000, "G", ins),               # break INS 110 bp at +30: 1
                (6000, "G" * 101, "G"),         # break DEL 140 bp: off by 29%
                (9000, "G", ins),               # callable, no kit SV: 0
                (12000, "G", ins)]              # a 12 bp kit indel 40 bp away
        st = self.states(keys, R=[(0, 15000)],
                         V=[(12040, 12040, "-", "cccccccccccc")],
                         svs=[("INS", 3031, 3031, 0, 110, "break"),
                              ("DEL", 6001, 6140, 140, 0, "break")])
        self.assertEqual(st, [1, -1, 0, -1])

    def test_sv_two_candidates(self):
        a = "G" + "T" * 100
        b = "G" + "T" * 99 + "C"
        # one kit SV without sequence, two graph alleles: no call for both
        st = self.states([(3000, "G", a), (3000, "G", b)], R=[(0, 15000)], V=[],
                         svs=[("INS", 3001, 3001, 0, 100, "break")])
        self.assertEqual(st, [-1, -1])
        # the kit's exact allele (a V line) picks one; the other is no call
        st = self.states([(3000, "G", a), (3000, "G", b)], R=[(0, 15000)],
                         V=[(3000, 3000, "-", "t" * 100)])
        self.assertEqual(st, [1, -1])

    def test_no_callable_lines_refused(self):
        pre = f"{self.d}/k"
        var_file(pre + ".var", [], [])
        sv_file(pre + ".sv.tsv", [])
        r = run([sys.executable, "-c",
                 "import sys; sys.path.insert(0, 'bin'); import panel_ancestral as P; "
                 f"P.read_kit_events('{pre}', 'G' * 100, 'NC_000962.3')"])
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("no R (callable) lines", r.stderr)


class KitEventsElement(_Base):
    """The tatC-helY element on the P200 design with the REAL route A files
    of L8 and ET1291 (analysis/outgroup_kit/out): both carry it, so it is
    ancestral-present (ALT) under every outgroup order. A second, one-carrier
    element allele at the same site (one base different) is not the kit's
    allele: no call there, so it stays as the panel alone says."""

    H37 = ("/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants/runroot/"
           "refbias/build/7713a8d71d8e-fix1/refs/GCF_000195955.fasta")
    STATES = ("/n/boslfs02/LABS/sfortune_lab/Lab/mchase/mtbc-graph-variants/"
              "analysis/panel_checks/states")
    EV = os.path.join(ROOT, "analysis", "outgroup_kit", "out")

    def setUp(self):
        super().setUp()
        need = [self.H37] + [f"{self.STATES}/{g}.npz" for g in (L8, ET)] + \
            [f"{self.EV}/{g}{s}" for g in (L8, ET) for s in (".var", ".sv.tsv")]
        if not all(os.path.exists(p) for p in need):
            self.skipTest("the real H37Rv, kit states or route A files are not here")

    def test_element_present_under_every_order(self):
        seq = "".join(l.strip() for l in open(self.H37) if not l.startswith(">"))
        self.assertEqual(seq[ELEMENT_POS - 1], "C")
        other = ELEMENT.replace("GTGTCGGG", "GTGTCTGG")   # one base, as in CX333
        samples = ["A1", "A2", "B1", "B2", L51, L61]
        vcf = panel_vcf(f"{self.d}/panel.vcf", samples,
                        f"{ELEMENT_POS} C {ELEMENT} 1 1 0 0 1 0\n"
                        f"{ELEMENT_POS} C {other} 0 0 1 0 0 0\n", len(seq))
        kit = [(L8, f"{self.STATES}/{L8}.npz"), (ET, f"{self.STATES}/{ET}.npz")]
        ev = ["--kit-events", f"{L8}={self.EV}/{L8}", "--kit-events", f"{ET}={self.EV}/{ET}"]
        tree = TatcHelyElement.TREE
        got, outs = set(), set()
        for order in itertools.permutations([L51, L61, L8, ET]):
            r, anc, pol = ancestral(self.d, tree, vcf, self.H37, list(order), kit, extra=ev)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            p = by_pos(pol)
            got.add((p[(ELEMENT_POS, "C", ELEMENT)]["ancestral"],
                     p.get((ELEMENT_POS, "C", other), {}).get("ancestral", "TIED")))
            outs.add(open(f"{self.d}/pol.tsv").read())
        self.assertEqual(len(outs), 1)
        # the main allele: present; the one-carrier variant: absent (REF)
        # from the panel's own genomes, the kit putting no constraint on it
        self.assertEqual(got, {("ALT", "REF")})

    def test_kit_events_need_kit_states(self):
        seq = "".join(l.strip() for l in open(self.H37) if not l.startswith(">"))
        vcf = panel_vcf(f"{self.d}/panel.vcf", ["A1", "A2", "B1", "B2", L51, L61, L8, ET],
                        f"{ELEMENT_POS} C {ELEMENT} 1 1 0 0 1 0 1 1\n", len(seq))
        r, _, _ = ancestral(self.d, TatcHelyElement.TREE, vcf, self.H37,
                            [L51, L61, L8, ET], extra=["--kit-events", f"{L8}={self.EV}/{L8}"])
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not a --kit genome", r.stderr)


# ==========================================================================
# 7. The event writer reads the new table
# ==========================================================================
class EventWriterReadsNewTable(unittest.TestCase):
    """The method column is tolerated, and a row that is not REF or ALT is
    ignored (left unpolarised), never read as REF-ancestral."""

    def setUp(self):
        p = os.environ.get("MTB_PY_VT", "")
        ok = p and os.access(p, os.X_OK) and \
            run([p, "-c", "import mtbvartools, dendropy"]).returncode == 0
        if not ok:
            self.skipTest("MTB_PY_VT with mtbvartools not set")
        self.py = p

    def test_rows(self):
        with tempfile.TemporaryDirectory() as d:
            write(f"{d}/m.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\t"
                  "QUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
                  "NC_000962.3\t500\th37rv:500:T>TA\tT\tTA\t.\tPASS\t"
                  "CLASS=small;FRAME=h37rv\tGT\t1\t1\t0\t0\n"
                  "NC_000962.3\t600\th37rv:600:G>GA\tG\tGA\t.\tPASS\t"
                  "CLASS=small;FRAME=h37rv\tGT\t1\t1\t0\t0\n")
            write(f"{d}/t.nwk", "(O:1,((S1:1,S2:1)x:1,(S3:1,S4:1)y:1)ing:1);\n")
            write(f"{d}/og.fa", ">O\nT\n")
            write(f"{d}/og.tsv", "column\tchrom\tpos\tref\talt\n0\tNC_000962.3\t500\tT\tC\n")
            pol = write(f"{d}/pol.tsv",
                        "chrom\tpos\tref\talt\toutgroup_gt\tancestral\tpanel_af\t"
                        "n_records\tmethod\n"
                        f"{PANSN}\t500\tT\tTA\t1\tALT\t0.5000\t1\toutgroups\n"
                        f"{PANSN}\t600\tG\tGA\t.\t.\t0.5000\t1\ttied\n")
            r = run([self.py, "assoc/bin/write_event_matrix.py", "--vcf", f"{d}/m.vcf",
                     "--tree", f"{d}/t.nwk", "--out", f"{d}/ev", "--outgroup-name", "O",
                     "--outgroup-fasta", f"{d}/og.fa", "--outgroup-sites", f"{d}/og.tsv",
                     "--panel-polarity", pol])
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            v = rows_of(f"{d}/ev/variants.tsv")
        self.assertEqual([x["polarity"] for x in v], ["alt_ancestral_outgroup", "unpolarised"])


if __name__ == "__main__":
    unittest.main()
