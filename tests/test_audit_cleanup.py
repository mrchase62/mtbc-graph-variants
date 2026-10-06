#!/usr/bin/env python3
"""Regression tests for the clean-up group (2026-10-06): no input taken
silently from CX333 or from a hand-made file, and one node key per base.

Same pattern as tests/run_tests.py: unittest, small synthetic inputs, the
standard library plus the pipeline's own scripts, seconds to run.

    $MTB_PY tests/test_audit_cleanup.py -v
"""
import csv
import hashlib
import importlib.util
import os
import shutil
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

H = "GCF_000195955#1#NC_000962.3"
VCFHDR = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t"
          "INFO\tFORMAT\tS\n")
NODEHDR = "node\taccession\tstart\tstrand\tn_occurrences\tlength\n"


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, rel))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def write(path, text, mode=None):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write(textwrap.dedent(text).lstrip("\n"))
    if mode:
        os.chmod(path, mode)
    return path


def rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def text(rel):
    return open(os.path.join(ROOT, rel)).read()


def code(rel):
    """The file without its comment lines, for checks on what it RUNS."""
    return "\n".join(l for l in text(rel).splitlines()
                     if not l.lstrip().startswith("#"))


def env(**kw):
    e = dict(os.environ)
    e.update(MTB_SITE_FILE="/dev/null", MTB_PY=sys.executable)
    for k in ("MTB_BUILD_DIR", "OG", "GRAPH_DIR", "BUILD_ID", "SLURM_ARRAY_TASK_ID",
              "REFMAP", "MTB_GRAPH_FRAMES", "ACC_CATALOGUE", "ANC_TREE", "ANC_ALN",
              "ANC_SITES", "ANC_OUTGROUPS"):
        e.pop(k, None)
    e.update({k: str(v) for k, v in kw.items()})
    return e


def run(cmd, cwd=None, **kw):
    return subprocess.run(cmd, cwd=cwd, env=env(**kw), capture_output=True,
                          text=True)


# ==========================================================================
# 2. Node offset orientation
# ==========================================================================
class ForwardOffset(unittest.TestCase):
    """odgi counts a node offset along the walk; the key uses the node's own
    forward offset, so R 201016,83,+ and H37Rv 201016,46,- (130 bp) agree."""

    def test_helper(self):
        n = load("mtb_norm_c", "bin/mtb_norm.py")
        self.assertEqual(n.forward_offset(83, "+", 130), 83)
        self.assertEqual(n.forward_offset(46, "-", 130), 83)
        self.assertEqual(n.forward_offset(0, "-", 1), 0)
        self.assertIsNone(n.forward_offset(46, "-", None))     # unknown, not 46
        self.assertEqual(n.forward_offset(n.forward_offset(5, "-", 9), "-", 9), 5)


class P4NodeKey(unittest.TestCase):
    """bin/p4_place.py keyed on odgi's walking offset: the same base got
    node:201016:83 from a reference walking the node + and node:201016:46
    from one walking it -."""

    def place(self, d, tag, node_field, lengths=True):
        # one direct-arm record in core sequence, so the table is never empty
        write(f"{d}/{tag}/direct.vcf", VCFHDR + "c\t50\t.\tA\tG\t50\t.\t.\tGT\t1\n")
        write(f"{d}/{tag}/matched.vcf", VCFHDR + "c\t200\t.\tA\tG\t50\t.\t.\tGT\t1\n")
        # off the H37Rv path (dist 100): a node-frame record
        write(f"{d}/{tag}/hpos.tsv", f"#src\ttgt\tdist\nR#1#c,199,+\t{H},599,+\t100\t+\t+\n")
        write(f"{d}/{tag}/npos.tsv", f"#src\tnode\nR#1#c,199,+\t{node_field}\n")
        write(f"{d}/mask.bed", "c\t1\t2\tx\n")
        write(f"{d}/loci.tsv", "pos\tlocus_id\n")
        write(f"{d}/graph.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\t"
              "QUAL\tFILTER\tINFO\tFORMAT\tOTHER\n")
        write(f"{d}/nodes.tsv", NODEHDR + "201016\tR\t1\t+\t1\t130\n")
        write(f"{d}/frames.tsv", "accession\tstrand\toffset\tpanel_len\n"
              "R\t+\t0\t1000\n")
        r = subprocess.run(
            [sys.executable, "bin/p4_place.py", "--sample", tag,
             "--reference", "R", "--build-id", "b",
             "--direct", f"{d}/{tag}/direct.vcf", "--matched", f"{d}/{tag}/matched.vcf",
             "--h37rv-pos", f"{d}/{tag}/hpos.tsv", "--node-pos", f"{d}/{tag}/npos.tsv",
             "--mask", f"{d}/mask.bed", "--loci", f"{d}/loci.tsv",
             "--graph-vcf", f"{d}/graph.vcf", "--out", f"{d}/{tag}/placed.tsv",
             "--frames", f"{d}/frames.tsv"]
            + (["--node-lengths", f"{d}/nodes.tsv"] if lengths else []),
            capture_output=True, text=True)
        return r

    def test_both_walks_one_key(self):
        with tempfile.TemporaryDirectory() as d:
            r1 = self.place(d, "fwd", "201016,83,+")
            r2 = self.place(d, "rev", "201016,46,-")
            self.assertEqual(r1.returncode, 0, r1.stderr + r1.stdout)
            self.assertEqual(r2.returncode, 0, r2.stderr + r2.stdout)
            k1 = [(q["key"], q["node_offset"]) for q in rows(f"{d}/fwd/placed.tsv")
                  if q["frame"] == "node"]
            k2 = [(q["key"], q["node_offset"]) for q in rows(f"{d}/rev/placed.tsv")
                  if q["frame"] == "node"]
        self.assertEqual(k1, [("node:201016:83", "83")])
        self.assertEqual(k2, [("node:201016:83", "83")])        # was node:201016:46

    def test_reverse_walk_without_length_is_not_keyed(self):
        # a forward-walked node needs no length; a reverse-walked one is
        # dropped and counted, never keyed on the walking offset
        with tempfile.TemporaryDirectory() as d:
            r = self.place(d, "rev", "201016,46,-", lengths=False)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            self.assertIn("no node length for a reverse-walked node: 1", r.stdout)
            self.assertEqual([q for q in rows(f"{d}/rev/placed.tsv")
                              if q["frame"] == "node"], [])

    def test_shell_passes_the_builds_node_table(self):
        s = code("bin/p4_place.sh")
        self.assertIn('NODES="${BUILD}/assets/node_positions.tsv"', s)
        self.assertIn('--node-lengths "$NODES"', s)


class IS6110NodeKey(unittest.TestCase):
    """is6110/bin/is6110_project_sites.py keyed an off-path site on odgi's
    walking offset, so carriers whose references walk the node in opposite
    directions got two keys for one insertion site."""

    def run_sites(self, d, nodes=True):
        write(f"{d}/frames.tsv",
              "accession\tpanel_len\trefs_len\tstrand\toffset\tagree\n"
              "R1\t100000\t100000\t+\t0\tok\n"
              "R2\t100000\t100000\t+\t0\tok\n"
              "GCF_000195955\t100000\t100000\t+\t0\tok\n")
        # R1 walks node 5 forward (offset 3), R2 in reverse (offset 6 of 10):
        # the same base. Both are 100 bp off the H37Rv path.
        write(f"{d}/odgi", f"""\
            #!{sys.executable}
            import sys
            a = sys.argv[1:]
            if a[0] == "paths" and "-L" in a:
                print("R1#1#c"); print("R2#1#c"); print("{H}")
            elif a[0] == "paths" and "-H" in a:
                print("path.name\\tpath.length\\tpath.step.count\\tnode.5")
                print("R1#1#c\\t10\\t5\\t1"); print("R2#1#c\\t10\\t5\\t1")
            elif a[0] == "position":
                q = [l.strip() for l in open(a[a.index("-F") + 1]) if l.strip()]
                for l in q:
                    if "-v" in a:
                        print(l + ("\\t5,3,+" if l.startswith("R1") else "\\t5,6,-"))
                    else:
                        print(l + "\\t{H},7000,+\\t100\\t+")
            """, mode=0o755)
        write(f"{d}/recon.tsv",
              "sample\tgeometry\tchrom_side\tclean_pos\torig_pos\treads\treads_q\n"
              "A\ttsd\tagreed\t5000\t5000\t9\t9\n"
              "B\ttsd\tagreed\t6000\t6000\t9\t9\n")
        write(f"{d}/refmap.tsv", "sample\treference\nA\tR1\nB\tR2\n")
        write(f"{d}/nodes.tsv", NODEHDR + "5\tR1\t11\t+\t1\t10\n5\tR2\t21\t-\t1\t10\n")
        e = dict(os.environ, MTB_GRAPH_FRAMES=f"{d}/frames.tsv")
        e.pop("MTB_BUILD_DIR", None)
        return subprocess.run(
            [sys.executable, "is6110/bin/is6110_project_sites.py",
             "--reconcile", f"{d}/recon.tsv", "--refmap", f"{d}/refmap.tsv",
             "--graph", "g.og", "--odgi", f"{d}/odgi", "--workdir", d,
             "--out", f"{d}/out.tsv", "--threads", "1"]
            + (["--node-lengths", f"{d}/nodes.tsv"] if nodes else []),
            capture_output=True, text=True, env=e)

    def test_both_walks_one_key(self):
        with tempfile.TemporaryDirectory() as d:
            r = self.run_sites(d)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            got = {q["sample"]: (q["key"], q["node_offset"]) for q in rows(f"{d}/out.tsv")}
        self.assertEqual(got, {"A": ("node:5:3", "3"),
                               "B": ("node:5:3", "3")})       # B was node:5:6

    def test_reverse_walk_without_node_table_refused(self):
        with tempfile.TemporaryDirectory() as d:
            r = self.run_sites(d, nodes=False)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("need their length", r.stderr)


class P5ReadsForwardOffset(unittest.TestCase):
    """bin/p5_states.py read a node key at start + offset along the sample's
    own reference, but the offset ran along the CARRIER's walk. For a
    reference walking the node in reverse the key's base is start + L-1-off.
    And the node tables were optional, with a CX333 default."""

    @classmethod
    def setUpClass(cls):
        cls.d = d = tempfile.mkdtemp()
        seq = "A" * 4_000_100
        os.makedirs(f"{d}/refs")
        for name in ("GCF_000195955", "RREF"):
            with open(f"{d}/refs/{name}.fasta", "w") as fh:
                fh.write(f">{name}\n")
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i + 80] + "\n")
        # node 7, 10 bp, walked `-` by RREF from panel position 101: forward
        # offset 2 is path position 101 + (10-1-2) = 108. The sample carries
        # a non-reference call at 103, where start + 2 would have read.
        write(f"{d}/keys.tsv",
              "key\tframe\th37rv_pos\tnode\tnode_offset\tcanonical_ref\t"
              "canonical_alt\tregion\tkind\tacc_locus\n"
              "node:7:2:A>C\tnode\t\t7\t2\tA\tC\toff_path_accessory\tSNP\t\n")
        write(f"{d}/placed.tsv", "sample\tframe\th37rv_pos\tref\talt\tr_pos\tnode\tnode_offset\n")
        write(f"{d}/s.g.vcf", VCFHDR
              + "c\t1\t.\tA\t<NON_REF>\t.\t.\tEND=102\tGT:DP\t0:30\n"
              + "c\t103\t.\tA\tC,<NON_REF>\t.\t.\t.\tGT:DP\t1:30\n"
              + "c\t104\t.\tA\t<NON_REF>\t.\t.\tEND=4000100\tGT:DP\t0:30\n")
        write(f"{d}/proj.pos", "H#1#c,0,+\tRREF#1#c,0,+\t0\t+\t+\n")
        write(f"{d}/frames.tsv", "accession\tstrand\toffset\tpanel_len\n"
              "RREF\t+\t0\t4000100\n")
        write(f"{d}/node_paths.tsv", "node\tpaths\n7\tRREF\n")
        write(f"{d}/node_positions.tsv", NODEHDR + "7\tRREF\t101\t-\t1\t10\n")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d, ignore_errors=True)

    def states(self, *extra):
        d = self.d
        return subprocess.run(
            [sys.executable, "bin/p5_states.py", "--sample", "S",
             "--reference", "RREF", "--keys", f"{d}/keys.tsv",
             "--placed", f"{d}/placed.tsv", "--gvcf", f"{d}/s.g.vcf",
             "--projected", f"{d}/proj.pos",
             "--h37rv", f"{d}/refs/GCF_000195955.fasta",
             "--out", f"{d}/s.states.tsv", *extra],
            capture_output=True, text=True,
            env=dict(os.environ, MTB_GRAPH_FRAMES=f"{d}/frames.tsv"))

    def test_reverse_walk_read_at_forward_offset(self):
        r = self.states("--node-paths", f"{self.d}/node_paths.tsv",
                        "--node-positions", f"{self.d}/node_positions.tsv")
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        cells = [l.split("\t")[1] for l in open(f"{self.d}/s.states.tsv")
                 if l[:1].isdigit()]
        # read at 108, a reference block: REF (sparse form omits REF).
        # Read at 103 it was NOCALL.
        self.assertEqual(cells, [])

    def test_node_tables_required(self):
        r = self.states()
        self.assertEqual(r.returncode, 2)
        self.assertIn("--node-positions", r.stderr)
        r = self.states("--node-paths", f"{self.d}/none.tsv",
                        "--node-positions", f"{self.d}/none.tsv")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("no node table", r.stderr)            # was only a note


def _rc(x):
    return x.translate(str.maketrans("ACGT", "TGCA"))[::-1]


class NodeForwardAlleles(unittest.TestCase):
    """D41: node-frame alleles were written in the reading reference's
    orientation, so one event read by references walking the node in
    opposite directions was two keys (C>T and G>A). They are now written on
    the node's forward strand."""

    def setUp(self):
        self.n = load("mtb_norm_d41", "bin/mtb_norm.py")

    def test_forward_unchanged(self):
        self.assertEqual(self.n.node_forward_alleles(83, "+", "C", "T"),
                         (83, "C", "T", "forward"))

    def test_snp_and_mnp_reversed(self):
        f = self.n.node_forward_alleles
        self.assertEqual(f(83, "-", "G", "A"), (83, "C", "T", "reversed"))
        # an MNP's first base as R reads it is its LAST on the forward strand
        self.assertEqual(f(84, "-", "GA", "TC"), (83, "TC", "GA", "reversed"))

    def test_indels_reanchored_on_the_forward_left(self):
        fwd = "ACGTACCGTAGGCTAACGTT"         # the node, forward
        rseq = "NNNNN" + _rc(fwd) + "NNNNN"  # R walks it `-`, from R pos 6
        L = len(fwd)
        def rpos(o):                        # R position of forward offset o
            return 6 + (L - 1 - o)
        # deletion of fwd[9] (forward: fwd[8] anchors, REF fwd[8:10])
        # R left-anchors it on fwd[10], complemented, at R pos rpos(10)
        ref, alt = _rc(fwd[9:11]), _rc(fwd[10])
        self.assertEqual(
            self.n.node_forward_alleles(10, "-", ref, alt, rseq, rpos(10)),
            (8, fwd[8:10], fwd[8], "reanchored"))
        # insertion of "GG" between fwd[8] and fwd[9]
        ref, alt = _rc(fwd[9]), _rc(fwd[9]) + "CC"
        self.assertEqual(
            self.n.node_forward_alleles(9, "-", ref, alt, rseq, rpos(9)),
            (8, fwd[8], fwd[8] + "GG", "reanchored"))

    def test_record_that_would_leave_the_node_is_unchanged(self):
        f = self.n.node_forward_alleles
        self.assertEqual(f(0, "-", "GT", "G", "ACGT", 1)[3], "off_node")
        self.assertEqual(f(5, "-", "GT", "G", "", 1), (5, "GT", "G",
                                                       "no_anchor_sequence"))


class P4NodeAllelesOneKey(unittest.TestCase):
    """End to end through bin/p4_place.py: the same SNP and the same deletion
    read by a reference walking the node + and one walking it - (or stored
    flipped in the panel) give one key and one allele."""

    FWD = "ACGTACCGTAGGCTAACGTTGCAT"          # node 9, forward

    def place(self, d, tag, ref_seq, rec_pos, ref, alt, walk_off, walk,
              flipped=False):
        write(f"{d}/{tag}/direct.vcf", VCFHDR + "c\t50\t.\tA\tG\t50\t.\t.\tGT\t1\n")
        write(f"{d}/{tag}/matched.vcf",
              VCFHDR + f"c\t{rec_pos}\t.\t{ref}\t{alt}\t50\t.\t.\tGT\t1\n")
        write(f"{d}/{tag}/hpos.tsv",
              f"#src\ttgt\tdist\nR#1#c,{rec_pos - 1},+\t{H},599,+\t100\t+\t+\n")
        write(f"{d}/{tag}/npos.tsv",
              f"#src\tnode\nR#1#c,{rec_pos - 1},+\t9,{walk_off},{walk}\n")
        write(f"{d}/{tag}/r.fasta", f">c\n{ref_seq}\n")
        write(f"{d}/{tag}/frames.tsv", "accession\tstrand\toffset\tpanel_len\n"
              f"R\t{'-' if flipped else '+'}\t0\t{len(ref_seq)}\n")
        write(f"{d}/mask.bed", "c\t1\t2\tx\n")
        write(f"{d}/loci.tsv", "pos\tlocus_id\n")
        write(f"{d}/graph.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\t"
              "QUAL\tFILTER\tINFO\tFORMAT\tOTHER\n")
        write(f"{d}/nodes.tsv", NODEHDR + f"9\tR\t1\t+\t1\t{len(self.FWD)}\n")
        r = subprocess.run(
            [sys.executable, "bin/p4_place.py", "--sample", tag,
             "--reference", "R", "--build-id", "b",
             "--direct", f"{d}/{tag}/direct.vcf", "--matched", f"{d}/{tag}/matched.vcf",
             "--h37rv-pos", f"{d}/{tag}/hpos.tsv", "--node-pos", f"{d}/{tag}/npos.tsv",
             "--mask", f"{d}/mask.bed", "--loci", f"{d}/loci.tsv",
             "--graph-vcf", f"{d}/graph.vcf", "--out", f"{d}/{tag}/placed.tsv",
             "--ref-fasta", f"{d}/{tag}/r.fasta", "--node-lengths", f"{d}/nodes.tsv",
             "--frames", f"{d}/{tag}/frames.tsv"],
            capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        return [(q["key"], q["ref"], q["alt"]) for q in rows(f"{d}/{tag}/placed.tsv")
                if q["frame"] == "node"]

    def both(self, o_fwd, ref_f, alt_f, o_rev, ref_r, alt_r, flipped=False):
        F, L = self.FWD, len(self.FWD)
        plus = "N" * 10 + F + "N" * 10           # R walks `+` from pos 11
        minus = "N" * 10 + _rc(F) + "N" * 10     # R walks `-` from pos 11
        with tempfile.TemporaryDirectory() as d:
            a = self.place(d, "fwd", plus, 11 + o_fwd, ref_f, alt_f, o_fwd, "+")
            if flipped:
                # stored flipped in the panel, walked `+` there: R's own
                # sequence still reads the node reverse complemented
                b = self.place(d, "rev", minus, 11 + o_rev, ref_r, alt_r,
                               L - 1 - o_rev, "+", flipped=True)
            else:
                b = self.place(d, "rev", minus, 11 + o_rev, ref_r, alt_r,
                               o_rev, "-")
        return a, b

    def test_snp(self):
        F, L = self.FWD, len(self.FWD)
        o = 7                                     # forward offset of the SNP
        a, b = self.both(o, F[o], "T", L - 1 - o, _rc(F[o]), _rc("T"))
        self.assertEqual(a, [(f"node:9:{o}", F[o], "T")])
        self.assertEqual(b, a)                    # was ("node:9:7", "C", "A")

    def test_snp_flipped_reference(self):
        F, L = self.FWD, len(self.FWD)
        o = 7
        a, b = self.both(o, F[o], "T", L - 1 - o, _rc(F[o]), _rc("T"),
                         flipped=True)
        self.assertEqual(b, a)

    def test_deletion(self):
        F, L = self.FWD, len(self.FWD)
        # forward: anchor F[8], F[9] deleted. Reverse: anchored on F[10].
        w = L - 1 - 10                            # R's walk offset of F[10]
        a, b = self.both(8, F[8:10], F[8], w, _rc(F[9:11]), _rc(F[10]))
        self.assertEqual(a, [("node:9:8", F[8:10], F[8])])
        self.assertEqual(b, a)


class P5NodeKeyOnForwardStrand(P5ReadsForwardOffset):
    """D41 in P5: a node key's REF is on the node's forward strand, and a
    reference walking the node `-` reads its complement. The sample's own
    call at the key's base equals the key's REF there only after
    complementing: a reversion, so REF, not NOCALL."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        d = cls.d
        # the key's base is path position 108 (see the parent); the sample
        # calls T there, the complement of the key's forward REF A
        write(f"{d}/s.g.vcf", VCFHDR
              + "c\t1\t.\tA\t<NON_REF>\t.\t.\tEND=107\tGT:DP\t0:30\n"
              + "c\t108\t.\tA\tT,<NON_REF>\t.\t.\t.\tGT:DP\t1:30\n"
              + "c\t109\t.\tA\t<NON_REF>\t.\t.\tEND=4000100\tGT:DP\t0:30\n")

    def test_reverse_walk_read_at_forward_offset(self):
        pass                                      # the parent's gVCF only

    def test_reversion_compared_on_the_reference_strand(self):
        r = self.states("--node-paths", f"{self.d}/node_paths.tsv",
                        "--node-positions", f"{self.d}/node_positions.tsv")
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        cells = [l.split("\t")[1] for l in open(f"{self.d}/s.states.tsv")
                 if l[:1].isdigit()]
        self.assertEqual(cells, [])               # REF; was NOCALL


# ==========================================================================
# 1. Old-graph inputs
# ==========================================================================
class P4P5ReadTheBuildsCollapsedVcf(unittest.TestCase):
    """P4 and P4b read all_variants.nolab.vcf.gz beside the graph, which a new
    build does not produce; P5's SV catalogue and sanity check defaulted to
    CX333 files; P4, P4b and P5 svgt globbed graphs/CX333 for the graph."""

    def test_shell_defaults(self):
        for p in ("bin/p4_place.sh", "bin/p4b_place_sv.sh", "bin/p5_svgt.sh"):
            s = code(p)
            self.assertNotIn("nolab", s, p)
            self.assertNotIn("graphs/CX333", s, p)
            self.assertIn('GRAPH_VCF="${GRAPH_VCF:-${BUILD}/assets/graph_collapsed.vcf.gz}"', s, p)
            self.assertIn("""OG="${OG:-$(awk -F'\\t' '$1=="graph"{print $2}' "${BUILD}/build_info.tsv")}\"""", s, p)

    def test_p0_writes_the_asset_they_read(self):
        self.assertIn('"${BUILD}/assets/graph_collapsed.vcf.gz"', code("bin/p0_prepare.sh"))

    def test_sanity_needs_the_graph_vcf_and_finish_passes_it(self):
        r = run([sys.executable, "bin/p5_sanity.py", "--matrix", "x"])
        self.assertEqual(r.returncode, 2)
        self.assertIn("--graph-vcf", r.stderr)
        self.assertIn('--graph-vcf "${BUILD}/assets/graph_collapsed.vcf.gz"',
                      code("bin/p5_finish.sh"))


class AssocChainRequiredInputs(unittest.TestCase):
    """add_outgroup.py and write_event_matrix.py defaulted to CX333's graph
    VCF, panel polarity table and node-locus table."""

    def test_add_outgroup_needs_panel_vcf(self):
        r = run([sys.executable, "assoc/bin/add_outgroup.py", "--alignment", "a",
                 "--sites", "s", "--out", "o"])
        self.assertEqual(r.returncode, 2)
        self.assertIn("--panel-vcf", r.stderr)
        self.assertNotIn("graphs/CX333", code("assoc/bin/add_outgroup.py"))

    def test_event_writer_needs_panel_polarity(self):
        r = run([sys.executable, "assoc/bin/write_event_matrix.py", "--vcf", "v",
                 "--tree", "t", "--out", "o"])
        if "No module named" in r.stderr:
            self.skipTest("the writer's imports are not available here")
        self.assertEqual(r.returncode, 2)
        self.assertIn("--panel-polarity", r.stderr)

    def test_event_writer_has_no_default_node_locus(self):
        s = code("assoc/bin/write_event_matrix.py")
        self.assertNotIn("refbias/assets/panel_polarity.tsv", s)
        self.assertNotIn("accessory/assets/node_locus.tsv", s)


class IS6110NoPilotIsmapper(unittest.TestCase):
    """is6110_place_by_flank.py and is6110_reconcile.py joined against the
    PILOT's ISMapper tables (refbias/p1f) unless told otherwise."""

    def test_defaults_empty(self):
        for p in ("is6110/bin/is6110_place_by_flank.py",
                  "is6110/bin/is6110_reconcile.py"):
            s = code(p)
            self.assertNotIn("refbias/p1f", s, p)
            self.assertIn('"--ismapper-dir", default=""', s, p)

    def test_reconcile_blank_when_join_off(self):
        # one element-side stack and a chromosome side that agrees with it;
        # a pilot-shaped ISMapper table in the working directory must not be
        # read, and the column is blank, not 0
        with tempfile.TemporaryDirectory() as d:
            write(f"{d}/el/A.elstacks.tsv",
                  "clean_pos\torig_pos\treads\treads_q\tboth_el_termini\tspan\t"
                  "el_start\tel_end\tsa_mapq_mean\n"
                  "5000\t5000\t20\t20\t1\t3\t1\t1\t60\n")
            write(f"{d}/el/A.junctions.tsv",
                  "pos\tis6110\tclips_start\tclips_end\treadthrough\tsa_at_is6110\n"
                  "5000\t1\t10\t10\t0\t20\n")
            write(f"{d}/refbias/p1f/A/A/IS6110/A__NC_000962.3_table.txt",
                  "x\ty\tcall\n5000\t5003\tKnown\n")
            r = subprocess.run(
                [sys.executable, os.path.join(ROOT, "is6110/bin/is6110_reconcile.py"),
                 "--elside-dir", f"{d}/el", "--junc-dir", f"{d}/el",
                 "--crossmap", "", "--check-against", "",
                 "--out", f"{d}/out.tsv", "--summary-out", f"{d}/sum.tsv"],
                cwd=d, capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            got = {q["ismapper"] for q in rows(f"{d}/out.tsv")}
        self.assertEqual(got, {""})                      # was {"1"}, the pilot's


class IscleanEveryPanelGenome(unittest.TestCase):
    """is6110/bin/p1i_build_matched.sh hard-coded CX333's refs/ and built only
    the refmap's references; P0's is6110_intervals step needs every panel
    genome of the build."""

    def test_builds_every_accession_from_the_build(self):
        with tempfile.TemporaryDirectory() as d:
            b = f"{d}/build/B1"
            write(f"{b}/build_info.tsv", "build_id\tB1\n")
            write(f"{b}/assets/accessions.txt", "A1\nA2\n")
            for a, s in (("A1", "ACGT" * 10), ("A2", "TTGA" * 10)):
                write(f"{b}/refs/{a}.fasta", f">{a}_c\n{s}\n")
                write(f"{d}/gff/{a}.is6110.gff", "##gff-version 3\n")
            # A1 already built from these bytes; A2 built from other bytes
            out = f"{d}/iso"
            for a in ("A1", "A2"):
                for ext in (".isclean.fasta", ".isclean.fasta.bwt", ".crossmap.tsv"):
                    write(f"{out}/{a}{ext}", "x\n")
            write(f"{out}/A1.ref.sha256", hashlib.sha256(
                open(f"{b}/refs/A1.fasta", "rb").read()).hexdigest() + "\n")
            write(f"{out}/A2.ref.sha256", "0" * 64 + "\n")
            log = f"{d}/built.txt"
            fake_py = write(f"{d}/fakepy", f"""\
                #!/bin/bash
                # stands in for is6110_build_isclean.py: record and write outputs
                echo "$@" >> {log}
                while [[ $# -gt 0 ]]; do
                  case "$1" in --out-fasta) printf '>c_isclean\\nAC\\n' > "$2";;
                               --out-crossmap) printf 'clean_pos\\tcum_deleted\\n' > "$2";; esac
                  shift
                done
                """, mode=0o755)
            fake = write(f"{d}/faketool", """\
                #!/bin/bash
                # bwa index / samtools faidx stand-in
                [[ $1 == index ]] && touch "$2.bwt"
                [[ $1 == faidx && $# -ge 3 ]] && printf '>c\\nAC\\n'
                exit 0
                """, mode=0o755)
            r = run(["bash", os.path.join(ROOT, "is6110/bin/p1i_build_matched.sh")],
                    cwd=d, SLURM_SUBMIT_DIR=d, MTB_BUILD_DIR=b, MTB_PY=fake_py,
                    MTB_BWA=fake, MTB_SAMTOOLS=fake, GFFDIR=f"{d}/gff", OUTDIR=out)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            built = open(log).read() if os.path.exists(log) else ""
            man = [q["reference"] for q in rows(f"{out}/manifest.tsv")]
            sha2 = open(f"{out}/A2.ref.sha256").read().strip()
        self.assertEqual(man, ["A1", "A2"])
        self.assertIn(f"{b}/refs/A2.fasta", built)          # rebuilt: other bytes
        self.assertNotIn("A1.fasta", built)                 # kept: same bytes
        self.assertEqual(sha2, hashlib.sha256(b">A2_c\n" + b"TTGA" * 10 + b"\n").hexdigest())


class P0AncestralAndCatalogue(unittest.TestCase):
    """P0 'ancestral' defaulted to data/trees/cx333.*; the accessory catalogue
    was copied from a hand-made accessory/assets file with no producer."""

    def setUp(self):
        self.d = d = tempfile.mkdtemp()
        for x in ("bin", "accessory"):
            os.symlink(os.path.join(ROOT, x), os.path.join(d, x))
        self.og = write(f"{d}/g/g.smooth.final.og", "graph bytes\n")
        gsha = hashlib.sha256(b"graph bytes\n").hexdigest()
        self.b = f"{d}/build/{gsha[:12]}"
        write(f"{self.b}/build_info.tsv", f"build_id\t{gsha[:12]}\ngraph_sha256\t{gsha}\n")
        write(f"{self.b}/assets/accessions.txt", f"A1\nA2\n{H.split('#')[0]}\n")
        write(f"{self.b}/logs/accessions.done", "2026-10-06\n")

    def tearDown(self):
        shutil.rmtree(self.d)

    def p0(self, *args, **kw):
        kw.setdefault("BUILD_ROOT", f"{self.d}/build")
        kw.setdefault("MTB_GATK_SIF", "/nonexistent.sif")
        return run(["bash", os.path.join(ROOT, "bin/p0_prepare.sh"), *args],
                   cwd=self.d, OG=self.og, **kw)

    def test_ancestral_has_no_default_inputs(self):
        # data/trees/cx333.* present, as in the working tree: not taken
        for f in ("cx333.rooted.nwk", "cx333.snps.fasta", "cx333.sites.tsv"):
            write(f"{self.d}/data/trees/{f}", "x\n")
        r = self.p0("--step", "ancestral")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("needs ANC_TREE, ANC_ALN and ANC_SITES", r.stderr)
        self.assertNotIn("cx333", code("bin/p0_prepare.sh"))

    def test_catalogue_made_from_the_builds_files(self):
        import gzip
        a = f"{self.b}/assets"
        write(f"{a}/accessory_loci.tsv",
              "locus_id\tpos\tklass\trep_len\tn_alleles\tcarriers_any\t"
              "carrier_frac\th37rv_cov\tnovelty\n"
              "ACC_5000\t5000\tpolymorphic\t400\t1\t2\t0.6\t0\tnovel\n")
        ins = "G" + "".join("ACGT"[(i * 7) % 4] for i in range(400))
        with gzip.open(f"{a}/graph_collapsed.vcf.gz", "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\t"
                     "FILTER\tINFO\tFORMAT\tA1\tA2\n"
                     f"c\t5000\t.\tG\t{ins}\t60\t.\t.\tGT\t1\t0\n")
        write(f"{self.b}/logs/assets.done", "2026-10-06\n")
        # review 2, R2-IS-1: the step measures h37rv_cov95 against the build's
        # H37Rv with blastn, so it needs refs (a stub blastn finding nothing)
        write(f"{self.b}/refs/GCF_000195955.fasta", ">NC_000962.3\n" + "ACGT" * 2000 + "\n")
        write(f"{self.b}/logs/refs.done", "2026-10-06\n")
        blastn = write(f"{self.d}/fakeblastn", "#!/bin/bash\nexit 0\n", mode=0o755)
        # a hand-made catalogue where the old step copied it from: not used
        write(f"{self.d}/accessory_assets/accessory_catalogue.tsv", "locus_id\nACC_9\n")
        fake = write(f"{self.d}/fakebwa", "#!/bin/bash\ntouch \"$2.bwt\"\n", mode=0o755)
        r = self.p0("--step", "catalogue", MTB_BWA=fake, MTB_BLASTN=blastn,
                    ACC_CATALOGUE_DIR=f"{self.d}/accessory_assets")
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        c = rows(f"{a}/accessory_catalogue.tsv")
        self.assertEqual([(x["locus_id"], x["graph_len"], x["panel_carriers"]) for x in c],
                         [("ACC_5000", "400", "A1")])
        self.assertEqual(c[0]["h37rv_cov95"], "0.0")
        self.assertTrue(open(f"{a}/accessory_catalogue.fasta").read().startswith(">ACC_5000"))
        self.assertTrue(os.path.exists(f"{a}/accessory_catalogue.fasta.bwt"))
        self.assertTrue(os.path.exists(f"{self.b}/logs/catalogue.done"))


class PresenceReadsTheBuildsCatalogue(unittest.TestCase):
    """The presence step's genotyper read accessory/assets/ (its own default);
    it now reads the build's catalogue, and the merge is passed it too."""

    def test_presence_gets_the_builds_catalogue(self):
        with tempfile.TemporaryDirectory() as d:
            b = f"{d}/build/B1"
            write(f"{b}/build_info.tsv", "build_id\tB1\n")
            for ext, t in (("tsv", "locus_id\nACC_1\n"), ("fasta", ">ACC_1\nAC\n")):
                write(f"{b}/assets/accessory_catalogue.{ext}", t)
            write(f"{d}/refbias/coh.crams.tsv", "S1\tx/S1.cram\n")
            refmap = write(f"{d}/refmap.tsv", "sample\treference\nS1\tR1\n")
            fake = write(f"{d}/fakepy", """\
                #!/bin/bash
                while [[ $# -gt 0 ]]; do
                  case "$1" in --out) out="$2";; esac
                  args="$args $1"; shift
                done
                echo "$args" > "$out"
                """, mode=0o755)
            os.makedirs(f"{d}/accessory")
            os.symlink(os.path.join(ROOT, "accessory", "bin"), f"{d}/accessory/bin")
            os.symlink(os.path.join(ROOT, "config"), f"{d}/config")
            r = run(["bash", os.path.join(ROOT, "accessory/bin/locus_presence_array.sh")],
                    cwd=d, SLURM_SUBMIT_DIR=d, SLURM_ARRAY_TASK_ID=1, COHORT_TAG="coh",
                    REFMAP=refmap, MTB_BUILD_DIR=b, MTB_PY_VT=fake,
                    MTB_CRAM_ROOT=d, MTB_H37RV=f"{d}/h.fa", OUTDIR=f"{d}/out")
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            args = open(f"{d}/out/S1.presence.tsv").read() \
                if os.path.exists(f"{d}/out/S1.presence.tsv") else ""
        self.assertIn(f"--catalogue {b}/assets/accessory_catalogue.tsv", args)
        self.assertIn(f"--fasta {b}/assets/accessory_catalogue.fasta", args)

    def test_merge_passed_the_builds_catalogue(self):
        self.assertIn('--accessory-catalogue "$_acccat"', code("bin/p5_finish.sh"))
        self.assertIn('_acccat="${BUILD}/assets/accessory_catalogue.tsv"',
                      code("bin/p5_finish.sh"))


class FrameConvertDocstring(unittest.TestCase):
    """frame_convert.py said odgi already accounts for inverted steps, the
    claim the leftovers group disproved."""

    def test_no_longer_claims_it(self):
        s = text("graphframe/bin/frame_convert.py")
        self.assertNotIn("already accounts for a locally inverted step", s)
        self.assertIn("does NOT already account for a locally inverted step", s)


if __name__ == "__main__":
    unittest.main()
