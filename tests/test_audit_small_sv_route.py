#!/usr/bin/env python3
"""Regression tests for two small-variant routing faults found by the 1K
evaluation (analysis/sv_eval_1k/README.md, problems 3 and 4; 2026-10-16).

PROBLEM 3. R's own SV-sized differences entered the small-variant matrix
through P4's composed arm (its inherited half) and through p5_states.py's
reading of a reference block as "the sample carries what R carries", and the
harmonization pass then wrote them CLASS=sv SVSOURCE=small_allele as if
HaplotypeCaller had called them: 486 insertions of 500 bp or more on
new1k_p200, 21,459 ALT cells. The SV block's inherited route (P4b, HANDOFF
ABS-2) owns them. Pinned: p5_keys.py makes no key of an inherited allele of
50 bp or more, p5_states.py does not read one as the sample's own record, and
an SV-sized insertion R carries is NOCALL `r_large_insertion`, not ALT.

PROBLEM 4. A deletion of 50 bp or more the sample did not call, with no gVCF
coverage at its anchor, was NOCALL `no_coverage` even where the flanks are
covered and nothing is inside: the deletion itself. It is now read in the
sample's P1 alignment to H37Rv (inherited_check.check_deletion).

    $MTB_PY tests/test_audit_small_sv_route.py -v
"""
import csv
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
pk = load("p5_keys", "bin/p5_keys.py")

KEY_COLS = ["key", "frame", "h37rv_pos", "node", "node_offset",
            "canonical_ref", "canonical_alt", "region", "kind", "acc_locus"]
PLACED_COLS = ["sample", "reference", "build_id", "arm", "component", "region",
               "frame", "key", "r_pos", "h37rv_pos", "dist_to_ref", "node",
               "node_offset", "frame_strand", "ref", "alt", "kind", "size",
               "acc_locus", "qual"]


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
    """reads: (start0, cigar tuples, mapq); sequence from `seq`."""
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"},
           "SQ": [{"SN": ctg, "LN": length}]}
    tmp = path + ".unsorted.bam"
    with pysam.AlignmentFile(tmp, "wb", header=hdr) as out:
        for i, (start, cig, mq) in enumerate(reads):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
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
    """Gapless reads of length rl starting every `step` bases in [lo, hi - rl)."""
    return [(p, [(0, rl)], mq) for p in range(lo, hi - rl, step)]


def write_tsv(path, cols, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    return path


class Case:
    """A one-sample p5_states.py run on a synthetic 4 Mb genome."""

    def __init__(self, d, h, r, keys, placed, proj, gvcf_blocks):
        self.d = d
        os.makedirs(f"{d}/refs")
        write_fasta(f"{d}/refs/GCF_000195955.fasta", "GCF_000195955", h)
        write_fasta(f"{d}/refs/RREF.fasta", "RREF", r)
        write_tsv(f"{d}/keys.tsv", KEY_COLS, keys)
        write_tsv(f"{d}/placed.tsv", PLACED_COLS, placed)
        with open(f"{d}/proj.pos", "w") as fh:
            for src, tgt in proj:
                fh.write(f"H#1#c,{src - 1},+\tRREF#1#c,{tgt - 1},+\t0\t+\t+\n")
        with open(f"{d}/s.g.vcf", "w") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\t"
                     "FILTER\tINFO\tFORMAT\tS\n")
            for s, e, dp in gvcf_blocks:
                fh.write(f"c\t{s}\t.\t{r[s - 1]}\t<NON_REF>\t.\t.\tEND={e}\tGT:DP\t0:{dp}\n")
        with open(f"{d}/frames.tsv", "w") as fh:
            fh.write("accession\tstrand\toffset\tpanel_len\nRREF\t+\t0\t4000000\n")
        write_tsv(f"{d}/np.tsv", ["node", "paths"], [])
        write_tsv(f"{d}/npos.tsv", ["node", "accession", "start", "strand",
                                    "n_occurrences", "length"], [])
        with open(f"{d}/p1.vcf", "w") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\t"
                     "FILTER\tINFO\tFORMAT\tS\n")
        with open(f"{d}/m.bed", "w") as fh:
            fh.write("c\t100\t200\tpe_ppe|x|y\n")
        self.h = h

    def run(self, extra, out="out.tsv"):
        d = self.d
        cmd = [sys.executable, "bin/p5_states.py", "--sample", "S",
               "--reference", "RREF", "--keys", f"{d}/keys.tsv",
               "--placed", f"{d}/placed.tsv", "--gvcf", f"{d}/s.g.vcf",
               "--projected", f"{d}/proj.pos",
               "--h37rv", f"{d}/refs/GCF_000195955.fasta",
               "--node-paths", f"{d}/np.tsv", "--node-positions", f"{d}/npos.tsv",
               "--out", f"{d}/{out}"] + extra
        p = subprocess.run(cmd, capture_output=True, text=True,
                           env=dict(os.environ, MTB_GRAPH_FRAMES=f"{d}/frames.tsv"))
        if p.returncode != 0:
            raise AssertionError(p.stderr + p.stdout)
        return {int(l.split("\t")[0]): l.rstrip("\n").split("\t")[1:]
                for l in open(f"{d}/{out}") if l[:1].isdigit()}, p.stdout

    def checked(self, bam):
        d = self.d
        return ["--h37rv-bam", bam, "--h37rv-vcf", f"{d}/p1.vcf",
                "--panel-diffs", f"{d}/none.tsv.gz", "--graph-vcf", f"{d}/none.vcf",
                "--repeat-mask", f"{d}/m.bed", "--samtools", SAMTOOLS,
                "--absent-detail", f"{d}/detail.tsv"]


def genome(seed=11):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(4_000_100))


# ---------------------------------------------------------------------------
class InheritedSvAllele(unittest.TestCase):
    """Problem 3: R's SV-sized alleles are the SV block's, not small variants."""

    def test_rule(self):
        f = pk.inherited_sv_allele
        self.assertTrue(f(dict(component="inherited", ref="A", alt="A" + "C" * 50)))
        self.assertTrue(f(dict(component="inherited", ref="A" + "C" * 60, alt="A")))
        self.assertFalse(f(dict(component="inherited", ref="A", alt="A" + "C" * 49)))
        # the sample's own HaplotypeCaller call is kept at any size
        self.assertFalse(f(dict(component="called", ref="A", alt="A" + "C" * 400)))
        # a long allele with a small net change is not SV-sized
        self.assertFalse(f(dict(component="inherited", ref="A" * 300, alt="C" * 290)))

    def test_keys_skip_inherited_sv_alleles(self):
        h = genome()
        d = tempfile.mkdtemp()
        try:
            write_fasta(f"{d}/h.fa", "GCF_000195955", h)
            ins = "G" * 80
            rows = [
                # R's 80 bp insertion, inherited: no key
                dict(sample="S", component="inherited", frame="h37rv",
                     h37rv_pos=3000, ref=h[2999], alt=h[2999] + ins, region="pe_ppe"),
                # R's 30 bp insertion, inherited: a key (a small variant)
                dict(sample="S", component="inherited", frame="h37rv",
                     h37rv_pos=7000, ref=h[6999], alt=h[6999] + "T" * 30, region="pe_ppe"),
                # a called 60 bp insertion: a key
                dict(sample="S", component="called", frame="h37rv",
                     h37rv_pos=9000, ref=h[8999], alt=h[8999] + "A" * 60, region="core"),
            ]
            os.makedirs(f"{d}/p4")
            write_tsv(f"{d}/p4/S.placed.tsv", PLACED_COLS, rows)
            write_tsv(f"{d}/refmap.tsv", ["sample", "reference"],
                      [dict(sample="S", reference="RREF")])
            p = subprocess.run([sys.executable, "bin/p5_keys.py", "--refmap",
                                f"{d}/refmap.tsv", "--dir", f"{d}/p4", "--h37rv",
                                f"{d}/h.fa", "--out", f"{d}/keys.tsv"],
                               capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr + p.stdout)
            got = {int(r["h37rv_pos"]) for r in csv.DictReader(
                open(f"{d}/keys.tsv"), delimiter="\t")}
            self.assertNotIn(3000, got)
            self.assertEqual(len(got), 2)
            self.assertIn("1 inherited records of >= 50 bp make no key", p.stdout)
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_states_r_carried_insertion_is_nocall_not_alt(self):
        h = genome()
        ins = "".join(random.Random(3).choice("ACGT") for _ in range(80))
        r = h[:3000] + ins + h[3000:]            # R carries 80 bp after 3000
        small = "TTAGG"
        r = r[:7080] + small + r[7080:]          # and 5 bp after H37Rv 7000
        d = tempfile.mkdtemp()
        try:
            keys = [dict(key=f"h37rv:3000:{h[2999]}>{h[2999] + ins}", frame="h37rv",
                         h37rv_pos=3000, canonical_ref=h[2999],
                         canonical_alt=h[2999] + ins, region="pe_ppe", kind="INDEL"),
                    dict(key=f"h37rv:7000:{h[6999]}>{h[6999] + small}", frame="h37rv",
                         h37rv_pos=7000, canonical_ref=h[6999],
                         canonical_alt=h[6999] + small, region="pe_ppe", kind="INDEL")]
            # the sample's placed table holds R's 80 bp allele as an inherited row
            placed = [dict(sample="S", reference="RREF", arm="composed",
                           component="inherited", frame="h37rv", h37rv_pos=3000,
                           ref=h[2999], alt=h[2999] + ins, region="pe_ppe",
                           kind="INDEL", size=80)]
            c = Case(d, h, r, keys, placed, [(3000, 3000), (7000, 7080)],
                     [(1, len(r), 30)])
            cells, out = c.run(["--no-absent-check"])
            # before the fix: ALT (the placed inherited row)
            self.assertEqual(cells[0], ["NOCALL", "", "r_large_insertion"])
            # a small insertion R carries is still ALT from its covered block
            self.assertEqual(cells[1][0], "ALT")
            self.assertIn("inherited records of >= 50 bp) not stated ALT: 1", out)
            # no placed row: R's allele read off the reference block, which
            # before the fix was ALT too (r_carries_alt)
            write_tsv(f"{d}/placed.tsv", PLACED_COLS, [])
            cells, out = c.run(["--no-absent-check"], out="out2.tsv")
            self.assertEqual(cells[0], ["NOCALL", "", "r_large_insertion"])
            self.assertEqual(cells[1][0], "ALT")
        finally:
            shutil.rmtree(d, ignore_errors=True)


# ---------------------------------------------------------------------------
class DeletionRuleUnit(unittest.TestCase):
    """check_deletion on a stub frame: flanks relative to the sample's depth."""

    class Frame:
        def __init__(self, flank, inside, baseline):
            import numpy as np
            self.depth = np.full(20_000, flank, dtype=np.int32)
            self.depth[10_000:10_100] = inside
            self.baseline = baseline
            self.min_bridge = max(5.0, 0.25 * baseline)
            self.bam, self.ctg = None, "c"

        def junction(self, s, e):
            return "unresolved", dict(bridge_start=0, bridge_end=0, clip_start=0,
                                      clip_end=0, joined=0)

    def setUp(self):
        self._md = ic._med_depth
        ic._med_depth = lambda *a, **k: 0.0           # no MAPQ-0 reads inside

    def tearDown(self):
        ic._med_depth = self._md

    def test_normal_flanks_no_depth_is_alt(self):
        st, why, _ = ic.check_deletion(self.Frame(60, 0, 80), ic.Mask([]), 10_000, 10_099, 5)
        self.assertEqual((st, why), ("ALT", "del_no_depth_inside"))

    def test_low_flanks_relative_to_the_sample_are_dropout(self):
        # flanks of 8x in an 80x sample: sequence that is barely covered anyway
        st, why, _ = ic.check_deletion(self.Frame(8, 0, 80), ic.Mask([]), 10_000, 10_099, 5)
        self.assertEqual((st, why), ("NOCALL", "del_depth_dropout"))

    def test_flanks_uncovered_stay_no_coverage(self):
        st, why, _ = ic.check_deletion(self.Frame(2, 0, 80), ic.Mask([]), 10_000, 10_099, 5)
        self.assertEqual((st, why), ("NOCALL", "no_coverage"))


# ---------------------------------------------------------------------------
class UncalledDeletionCheck(unittest.TestCase):
    """Problem 4: a deletion with covered flanks and nothing inside is ALT."""

    def setUp(self):
        need_reads(self)

    def run_case(self, reads, mask_bed=None):
        h = genome(5)
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        dele = h[4999:5100]                      # deletes H37Rv 5001..5100
        keys = [dict(key=f"h37rv:5000:{dele}>{h[4999]}", frame="h37rv",
                     h37rv_pos=5000, canonical_ref=dele, canonical_alt=h[4999],
                     region="core", kind="INDEL")]
        # the anchor's gVCF block reaches into the deleted bases: DP 0
        c = Case(d, h, h, keys, [], [(5000, 5000)],
                 [(1, 4995, 30), (4996, 5102, 0), (5103, len(h), 30)])
        if mask_bed:
            with open(f"{d}/m.bed", "w") as fh:
                fh.write(mask_bed)
        bam = write_bam(f"{d}/p1.bam", "c", len(h),
                        reads + tiled(200_000, 201_000), h)   # some depth elsewhere
        old, _ = c.run(["--no-absent-check"], out="old.tsv")
        self.assertEqual(old[0], ["NOCALL", "", "no_coverage"])
        cells, out = c.run(c.checked(bam))
        self.assertIn("#absent_check\tv2", open(f"{d}/out.tsv").read())
        return cells[0]

    def test_covered_flanks_nothing_inside_is_alt(self):
        st = self.run_case(tiled(3800, 4986) + tiled(5100, 6300))
        self.assertEqual((st[0], st[2]), ("ALT", "del_no_depth_inside"))

    def test_split_reads_joining_the_ends_are_alt(self):
        # reads with the 100 bp deletion in the CIGAR, from end to end
        rd = [(p, [(0, 4999 - p), (2, 100), (0, 100 - (4999 - p))], 60)
              for p in range(4905, 4995, 5)]
        st = self.run_case(tiled(3800, 4986) + tiled(5100, 6300) + rd)
        self.assertEqual((st[0], st[2]), ("ALT", "del_junction_joined"))

    def test_uncovered_flanks_stay_nocall(self):
        st = self.run_case([])
        self.assertEqual(st, ["NOCALL", "", "no_coverage"])

    def test_spanned_with_depth_is_ref(self):
        st = self.run_case(tiled(3800, 6300))
        self.assertEqual(st[0], "REF")
        self.assertEqual(st[2], "del_spanned_ref")

    def test_depth_alone_in_the_mask_is_not_a_deletion(self):
        # GC-rich PE_PGRS: no depth inside, low flanks, in the repeat mask --
        # REF in every marin36_p200 assembly where the old rule said ALT
        st = self.run_case(tiled(3800, 4986) + tiled(5100, 6300),
                           mask_bed="c\t5020\t5080\tpe_ppe|PE_PGRS|x\n")
        self.assertEqual(st, ["NOCALL", "", "del_depth_dropout"])

    def test_split_reads_decide_in_the_mask(self):
        rd = [(p, [(0, 4999 - p), (2, 100), (0, 100 - (4999 - p))], 60)
              for p in range(4905, 4995, 5)]
        st = self.run_case(tiled(3800, 4986) + tiled(5100, 6300) + rd,
                           mask_bed="c\t5020\t5080\tpe_ppe|PE_PGRS|x\n")
        self.assertEqual((st[0], st[2]), ("ALT", "del_junction_joined"))

    def test_mapq0_inside_is_a_repeat(self):
        inside = [(p, [(0, 50)], 0) for p in range(5000, 5050, 2)]
        st = self.run_case(tiled(3800, 4986) + tiled(5100, 6300) + inside)
        self.assertEqual(st, ["NOCALL", "", "del_repeat_mapq0"])


if __name__ == "__main__":
    unittest.main()
