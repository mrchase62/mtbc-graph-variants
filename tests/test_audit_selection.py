#!/usr/bin/env python3
"""Regression tests for the 2026-10-05 audit, group "selection".

GRAPHVCF-4 / P0P2-5 (reference selection keyed by position), GRAPHVCF-5 (the
collapse pads REF and loses carriers), PGB-9/PGB-10 (the collapsed VCF as the
product; production's vcfwave settings), GRAPHVCF-6 (panel AF by position),
P0P2-8 (P2's summary reads the pilot's P1 folder) and P0P2-11 (CRLF tables).

Same pattern as tests/run_tests.py: small synthetic inputs, standard library
plus the pipeline's own scripts, seconds to run.

    $MTB_PY tests/test_audit_selection.py
"""
import csv
import gzip
import os
import re
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


def write(path, text):
    with open(path, "w") as fh:
        fh.write(textwrap.dedent(text).lstrip("\n"))
    return path


def vcf(path, samples, records, extra_header=""):
    """records: (pos, ref, alt, info, [gt per sample])"""
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##contig=<ID=c,length=1000>\n")
        fh.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="">\n')
        fh.write('##INFO=<ID=AN,Number=1,Type=Integer,Description="">\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="">\n')
        fh.write(extra_header)
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                 + "\t".join(samples) + "\n")
        for pos, ref, alt, info, gts in records:
            fh.write(f"c\t{pos}\t.\t{ref}\t{alt}\t60\t.\t{info or '.'}\tGT\t"
                     + "\t".join(gts) + "\n")
    return path


def rows(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def py():
    return sys.executable


# ---------------------------------------------------------------- GRAPHVCF-4

class SelectionByAllele(unittest.TestCase):
    """t8_select_reference.py keyed the panel by position and ignored the ALT."""

    SAMPLES = ["gA", "gG"]

    def select(self, panel_records, isolate_records, depth=None, top=2,
               full=False, check=True, h37rv=""):
        """`depth`: {pos: reads}; default 30 at every panel site.
        isolate_records: (pos, ref, alt, gt[, filter]). `h37rv`: the name
        of the all-REF H37Rv candidate (D21); empty leaves it out."""
        with tempfile.TemporaryDirectory() as d:
            panel = vcf(os.path.join(d, "panel.vcf"), self.SAMPLES, panel_records)
            iso = os.path.join(d, "iso.vcf")
            with open(iso, "w") as fh:
                fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
                for rec in isolate_records:
                    pos, ref, alt, gt = rec[:4]
                    flt = rec[4] if len(rec) > 4 else "."
                    fh.write(f"c\t{pos}\t.\t{ref}\t{alt}\t50\t{flt}\t.\tGT:DP\t{gt}:30\n")
            dep = os.path.join(d, "depth.tsv")
            if depth is None:
                depth = {r[0]: 30 for r in panel_records}
            with open(dep, "w") as fh:
                for p, n in sorted(depth.items()):
                    fh.write(f"c\t{p}\t{n}\n")
            out = os.path.join(d, "cand.tsv")
            r = subprocess.run([py(), "bin/t8_select_reference.py", "--vcf", iso,
                                "--panel-snps", panel, "--depth", dep,
                                "--out", out, "--top", str(top),
                                "--h37rv-name", h37rv],
                               capture_output=True, text=True)
            if check:
                self.assertEqual(r.returncode, 0, r.stderr)
            else:
                return r
            got = rows(out)
            if full:
                return got
            return {x["reference"]: int(x["snp_distance"]) for x in got}

    def test_split_multiallelic_site_scored_per_allele(self):
        # C>A and C>G at one position, two rows; the isolate carries C>A.
        # Keyed by position the C>G row survived, the isolate "matched" it and
        # gG (the wrong allele) came out nearest.
        d = self.select([(100, "C", "A", "", ["1", "0"]),
                         (100, "C", "G", "", ["0", "1"])],
                        [(100, "C", "A", "1")])
        self.assertEqual(d, {"gA": 0, "gG": 2})

    def test_padded_snp_is_trimmed_and_kept(self):
        # CG>TG is the SNP C>T; it was dropped by the len(REF)==1 filter
        d = self.select([(200, "CG", "TG", "", ["1", "0"])],
                        [(200, "C", "T", "1")])
        self.assertEqual(d, {"gA": 0, "gG": 1})

    def test_isolate_genotype_zero_is_not_a_call(self):
        d = self.select([(300, "C", "T", "", ["1", "0"])],
                        [(300, "C", "T", "0")])
        self.assertEqual(d, {"gA": 1, "gG": 0})


class SelectionOverCalledSites(SelectionByAllele):
    """The user's decisions D20 and D24 (2026-10-07): compare only sites both
    sides called, and rank per compared site. An uncovered isolate site
    counted as REF, and a genome with more missing cells looked closer."""

    def test_uncovered_isolate_site_is_not_ref(self):
        # gA carries 100 and 200, gG neither; the isolate carries 200 and has
        # no reads at 100. Counted as REF at 100 it was 1 from both genomes;
        # over the sites it covers, gA is 0 and nearest.
        d = self.select([(100, "C", "T", "", ["1", "0"]),
                         (200, "C", "T", "", ["1", "0"])],
                        [(200, "C", "T", "1")], depth={100: 0, 200: 30})
        self.assertEqual(d, {"gA": 0, "gG": 1})

    def test_ranked_per_compared_site(self):
        # gA: 3 mismatches over 20 sites (0.15); gG: 2 over the 11 it has
        # genotyped (0.18). The raw count chose gG; per site gA is nearer.
        recs = []
        for i in range(20):
            ga = "1" if i < 3 else "0"
            gg = "1" if i < 2 else ("." if i >= 11 else "0")
            recs.append((100 + i, "C", "T", "", [ga, gg]))
        got = self.select(recs, [], full=True)
        self.assertEqual([x["reference"] for x in got], ["gA", "gG"])
        self.assertEqual([x["snp_distance"] for x in got], ["3", "2"])
        self.assertEqual([x["n_compared"] for x in got], ["20", "11"])
        self.assertEqual(got[0]["distance_per_site"], "0.15")

    def test_filtered_or_indel_record_is_not_ref_evidence(self):
        # the isolate has a filtered call at 100 and a deletion over 300;
        # neither site is REF evidence, so only 200 is compared
        d = self.select([(100, "C", "T", "", ["1", "0"]),
                         (200, "C", "T", "", ["1", "0"]),
                         (301, "C", "T", "", ["1", "0"])],
                        [(100, "C", "T", "1", "LowQual"), (200, "C", "T", "1"),
                         (300, "ACG", "A", "1")])
        self.assertEqual(d, {"gA": 0, "gG": 1})

    def test_sparse_genome_is_not_a_candidate(self):
        # gG is genotyped at 1 of 10 sites, under half gA's: not ranked, even
        # though its one site matches
        recs = [(100 + i, "C", "T", "", ["1" if i else "0", "0" if i == 0 else "."])
                for i in range(10)]
        got = self.select(recs, [], full=True)
        self.assertEqual([x["reference"] for x in got], ["gA"])

    def test_no_coverage_input_is_fatal(self):
        with tempfile.TemporaryDirectory() as d:
            panel = vcf(os.path.join(d, "p.vcf"), self.SAMPLES,
                        [(100, "C", "T", "", ["1", "0"])])
            iso = vcf(os.path.join(d, "i.vcf"), ["S"], [])
            r = subprocess.run([py(), "bin/t8_select_reference.py", "--vcf", iso,
                                "--panel-snps", panel, "--out",
                                os.path.join(d, "o.tsv")],
                               capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("--bam or --depth", r.stderr)

    def test_p1_passes_the_bam(self):
        s = open("bin/p1_select_reference.sh").read()
        self.assertIn('--bam "$H37BAM"', s)

    def test_summary_gap_is_per_site(self):
        import importlib.util
        spec = importlib.util.spec_from_file_location("p1s", "bin/p1_summary.py")
        m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
        # counts say nxt is 1 CLOSER; per site it is 0.0001 * 10,000 = 1 farther
        b = dict(snp_distance="100", n_compared="10000", distance_per_site="0.01")
        n = dict(snp_distance="99", n_compared="9800", distance_per_site="0.0101")
        self.assertEqual(m.snp_gap(b, n), 1)
        self.assertEqual(m.snp_gap(dict(snp_distance="3"), dict(snp_distance="5")), 2)


# ------------------------------------------------------ GRAPHVCF-5, PGB-9

def collapse_py():
    """The collapse program embedded in bin/vcf_collapse.sh."""
    s = open("bin/vcf_collapse.sh").read()
    m = re.search(r"<<'PY' \|\| true\n(.*?)\nPY\n", s, re.S)
    if not m:
        raise AssertionError("no COLLAPSE_PY program in bin/vcf_collapse.sh")
    return m.group(1)


# bcftools for the end-to-end collapse tests: the host's, or the lab's conda
# one (MTB_BCFTOOLS_DIR overrides), put on PATH for the script
BCF_DIR = os.environ.get("MTB_BCFTOOLS_DIR") or next(
    (d for d in ("/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/mtb_isolates/bin",)
     if os.path.exists(os.path.join(d, "bcftools"))), "")
HAVE_BCF = bool(shutil.which("bcftools") or BCF_DIR)


def bcf_env():
    return dict(os.environ, PATH=BCF_DIR + os.pathsep + os.environ["PATH"]
                if BCF_DIR else os.environ["PATH"])


def ref_fasta(path, bases, length=1000, fill="T"):
    """A one-contig FASTA for the collapse's left-alignment; `bases` is
    {1-based pos: sequence}."""
    seq = list(fill * length)
    for p, b in bases.items():
        seq[p - 1:p - 1 + len(b)] = b
    with open(path, "w") as fh:
        fh.write(">anything\n" + "".join(seq) + "\n")
    return path


def parse_vcf_text(text):
    out = {}
    for line in text.splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        out[(int(f[1]), f[3], f[4])] = f[9:]
    return out


class CollapseKeepsAllelesTrimmed(unittest.TestCase):
    """The norm -m +any / -m -any round-trip padded REF and lost carriers."""

    RECORDS = [
        # one SNP split across two vcfwave records: carriers must be unioned
        (10, "C", "A", "", ["1", "0", "0"]),
        (10, "C", "A", "", ["0", "1", "."]),
        # a SNP and an insertion on the same base, both carried by sample s0:
        # a merged haploid record can hold only one of them
        (20, "G", "T", "", ["1", "0", "0"]),
        (20, "G", "GAA", "", ["1", "1", "0"]),
        # a padded SNP written by another allele's REF length
        (30, "CG", "TG", "", ["0", "0", "1"]),
        (30, "CGA", "C", "", ["1", "0", "0"]),
    ]
    WANT = {
        (10, "C", "A"): ["1", "1", "0"],
        (20, "G", "T"): ["1", "0", "0"],
        (20, "G", "GAA"): ["1", "1", "0"],
        (30, "C", "T"): ["0", "0", "1"],
        (30, "CGA", "C"): ["1", "0", "0"],
    }

    def test_collapse_program(self):
        with tempfile.TemporaryDirectory() as d:
            src = vcf(os.path.join(d, "in.vcf"), ["s0", "s1", "s2"], self.RECORDS)
            r = subprocess.run([py(), "-c", collapse_py()], stdin=open(src),
                               capture_output=True, text=True, check=True)
        self.assertEqual(parse_vcf_text(r.stdout), self.WANT)

    @unittest.skipUnless(HAVE_BCF, "no bcftools")
    def test_collapse_script_end_to_end(self):
        # the whole script, as the build runs it: no duplicate key, no padded
        # SNP, every carrier kept, AC refilled
        with tempfile.TemporaryDirectory(dir="/tmp") as d:
            src = vcf(os.path.join(d, "in.vcf"), ["s0", "s1", "s2"], self.RECORDS)
            fa = ref_fasta(os.path.join(d, "h.fa"),
                           {10: "C", 20: "G", 30: "CGA"})
            out = os.path.join(d, "out.vcf.gz")
            subprocess.run(["bash", "bin/vcf_collapse.sh", "TEST", "--in", src,
                            "--out", out, "--ref", fa], check=True,
                           capture_output=True, text=True, env=bcf_env())
            text = gzip.open(out, "rt").read()
        got = parse_vcf_text(text)
        self.assertEqual(got, self.WANT)
        info = [l.split("\t")[7] for l in text.splitlines()
                if l.startswith("c\t10\t")][0]
        self.assertIn("AC=2", info.split(";"))


@unittest.skipUnless(HAVE_BCF, "no bcftools")
class CollapseLeftAligns(unittest.TestCase):
    """The user's decision D22 (2026-10-07): the graph VCF is left-aligned
    against H37Rv before the union, so one event vcfwave wrote at two
    positions is one record, at the position every cohort key uses."""

    # G at 100, then AAAA at 101-104, T at 105
    BASES = {100: "GAAAA"}

    def collapse(self, d, records, ref=True, bases=None):
        src = vcf(os.path.join(d, "in.vcf"), ["s0", "s1", "s2"], records)
        out = os.path.join(d, "out.vcf.gz")
        cmd = ["bash", "bin/vcf_collapse.sh", "TEST", "--in", src, "--out", out]
        if ref:
            cmd += ["--ref", ref_fasta(os.path.join(d, "h.fa"),
                                       bases or self.BASES)]
        env = bcf_env()
        # project_env.sh fills an unset MTB_REF_FASTA with the working tree's
        # H37Rv, so "no reference" is a path that does not exist
        env["MTB_REF_FASTA"] = os.path.join(d, "no_such_ref.fa")
        r = subprocess.run(cmd, capture_output=True, text=True, env=env)
        return r, (parse_vcf_text(gzip.open(out, "rt").read())
                   if r.returncode == 0 else None)

    def test_one_deletion_at_two_positions_is_one_record(self):
        with tempfile.TemporaryDirectory(dir="/tmp") as d:
            r, got = self.collapse(d, [(102, "AA", "A", "", ["1", "0", "0"]),
                                       (103, "AA", "A", "", ["0", "1", "0"])])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(got, {(100, "GA", "G"): ["1", "1", "0"]})

    def test_insertion_left_aligned(self):
        with tempfile.TemporaryDirectory(dir="/tmp") as d:
            r, got = self.collapse(d, [(104, "A", "AA", "", ["0", "0", "1"])])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(got, {(100, "G", "GA"): ["0", "0", "1"]})

    def test_snp_untouched(self):
        with tempfile.TemporaryDirectory(dir="/tmp") as d:
            r, got = self.collapse(d, [(103, "A", "C", "", ["1", "0", "0"])])
        self.assertEqual(got, {(103, "A", "C"): ["1", "0", "0"]})

    def test_no_reference_is_fatal(self):
        with tempfile.TemporaryDirectory(dir="/tmp") as d:
            r, _ = self.collapse(d, [(103, "A", "C", "", ["1", "0", "0"])],
                                 ref=False)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("left-alignment needs H37Rv", r.stderr)

    def test_ref_disagreeing_with_h37rv_is_fatal(self):
        with tempfile.TemporaryDirectory(dir="/tmp") as d:
            r, _ = self.collapse(d, [(103, "C", "T", "", ["1", "0", "0"])])
        self.assertNotEqual(r.returncode, 0)

    def test_wrong_length_reference_is_fatal(self):
        with tempfile.TemporaryDirectory(dir="/tmp") as d:
            src = vcf(os.path.join(d, "in.vcf"), ["s0"],
                      [(103, "A", "C", "", ["1"])])
            fa = os.path.join(d, "h.fa")
            with open(fa, "w") as fh:
                fh.write(">x\n" + "A" * 999 + "\n")
            r = subprocess.run(["bash", "bin/vcf_collapse.sh", "TEST", "--in",
                                src, "--out", os.path.join(d, "o.vcf.gz"),
                                "--ref", fa], capture_output=True, text=True,
                               env=bcf_env())
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("999 bp", r.stderr)

    def test_decompose_passes_the_reference(self):
        s = open("bin/vcf_decompose.sh").read()
        self.assertIn('--in "$OUT" --out "$COLLAPSED" --ref "$REF_FA"', s)


class DecomposeMatchesProduction(unittest.TestCase):
    """PGB-10: vcf_decompose.sh defaulted to -I 64 and had no long-allele
    fallback; PGB-9: its own inline collapse duplicated vcf_collapse.sh."""

    def setUp(self):
        self.s = open("bin/vcf_decompose.sh").read()

    def test_production_inv_min_and_fallback(self):
        self.assertRegex(self.s, r"(?m)^WAVE_I=1000$")
        self.assertRegex(self.s, r"(?m)^MAX_INV_ALLELE=5000$")
        self.assertRegex(self.s, r"(?m)^FALLBACK_I=1000$")
        self.assertIn('big_*) iv="$FALLBACK_I"', self.s)
        self.assertIn("##MTB_decompose=", self.s)

    def test_partition_awk_splits_on_longest_allele(self):
        awk = re.search(r"awk -F'\\t' -v cap=\"\$MAX_INV_ALLELE\" \\\n.*?'\n(.*?)'\n",
                        self.s, re.S)
        self.assertIsNotNone(awk)
        with tempfile.TemporaryDirectory() as d:
            body = ("c\t1\t.\tA\tC\n"
                    f"c\t2\t.\tA\tC,{'G' * 12}\n"
                    f"c\t3\t.\t{'T' * 11}\tA\n")
            subprocess.run(["awk", "-F\t", "-v", "cap=10",
                            "-v", f"small={d}/s", "-v", f"big={d}/b",
                            awk.group(1)], input=body, text=True, check=True)
            small = open(f"{d}/s").read().splitlines()
            big = open(f"{d}/b").read().splitlines()
        self.assertEqual([l.split("\t")[1] for l in small], ["1"])
        self.assertEqual([l.split("\t")[1] for l in big], ["2", "3"])

    def test_collapse_is_delegated(self):
        self.assertNotIn("norm -m +any", self.s)
        self.assertIn('bash "$COLLAPSE_SH"', self.s)


# ---------------------------------------------------------------- GRAPHVCF-6

class PanelAfByAllele(unittest.TestCase):
    """p5_matrix.py gave a key the max AF of any record at its position."""

    def test_sites_panel_af_per_allele(self):
        with tempfile.TemporaryDirectory() as d:
            gv = vcf(os.path.join(d, "graph.vcf"), ["g1", "g2", "g3", "g4"], [
                (100, "C", "A", "AC=1;AN=4", ["1", "0", "0", "0"]),
                # an insertion on the same base, carried by most of the panel
                (100, "C", "CT", "AC=3;AN=4", ["0", "1", "1", "1"]),
                # padded SNP: the key is G>T at 200
                (200, "GA", "TA", "AC=3;AN=4", ["1", "1", "1", "0"]),
            ])
            keys = os.path.join(d, "keys.tsv")
            with open(keys, "w") as fh:
                fh.write("key\tframe\th37rv_pos\tnode\tnode_offset\tcanonical_ref\t"
                         "canonical_alt\tregion\tkind\tacc_locus\n")
                for pos, r, a in ((100, "C", "A"), (100, "C", "G"), (200, "G", "T")):
                    fh.write(f"h37rv:{pos}:{r}>{a}\th37rv\t{pos}\t\t\t{r}\t{a}\t"
                             f"core\tSNP\t\n")
            write(os.path.join(d, "refmap.tsv"), "sample\treference\nS1\tR\n")
            write(os.path.join(d, "S1.states.tsv"), """
                key\tstate
                h37rv:100:C>A\tALT
                h37rv:100:C>G\tALT
                h37rv:200:G>T\tALT
                """)
            sites = os.path.join(d, "sites.tsv")
            subprocess.run([py(), "bin/p5_matrix.py", "--refmap",
                            os.path.join(d, "refmap.tsv"), "--keys", keys,
                            "--dir", d, "--out", os.devnull, "--no-dense",
                            "--graph-vcf", gv, "--sites-out", sites],
                           check=True, capture_output=True, text=True)
            got = {r["key"]: (r["panel_af"], r["h37rv_minor"]) for r in rows(sites)}
        self.assertEqual(got["h37rv:100:C>A"], ("0.25", "0"))   # was 0.75, 1
        self.assertEqual(got["h37rv:100:C>G"], ("", ""))        # not in the panel
        self.assertEqual(got["h37rv:200:G>T"], ("0.75", "1"))   # was found by pos


# -------------------------------------------------------- P0P2-8, P0P2-11

class P2SummaryInputs(unittest.TestCase):

    def test_p2_call_passes_p1_work(self):
        s = open("bin/p2_call.sh").read()
        call = s[s.index("bin/p2_summary.py"):]
        call = call[:call.index("\n    ;;") if "\n    ;;" in call else call.index("\nfi")]
        self.assertIn('--p1-work "$P1WORK"', call)
        self.assertIn('P1WORK="${P1WORK:-', s)

    def test_p2_summary_lf_and_reads_given_p1_work(self):
        with tempfile.TemporaryDirectory() as d:
            p2, p1w, work = (os.path.join(d, x) for x in ("p2", "p1w", "w"))
            for x in (p2, p1w, work):
                os.makedirs(x)
            write(os.path.join(d, "refmap.tsv"),
                  "sample\treference\tsnp_distance\nS1\tR\t5\n")
            rec = "#CHROM\tPOS\nc\t1\n"
            for p in (os.path.join(p2, "S1.vcf.gz"), os.path.join(p1w, "S1.h37rv.vcf.gz")):
                with gzip.open(p, "wt") as fh:
                    fh.write(rec + ("c\t2\n" if "p1w" in p else ""))
            write(os.path.join(p2, "S1.delly.vcf"), "#h\n")
            out = os.path.join(d, "sum.tsv")
            subprocess.run([py(), "bin/p2_summary.py", "--refmap",
                            os.path.join(d, "refmap.tsv"), "--dir", p2,
                            "--work", work, "--p1-work", p1w, "--out", out],
                           check=True, capture_output=True, text=True)
            raw = open(out, "rb").read()
            got = rows(out)[0]["h37rv_small"]
        self.assertNotIn(b"\r", raw)
        self.assertEqual(got, "2")

    def test_p1_summary_lf(self):
        with tempfile.TemporaryDirectory() as d:
            write(os.path.join(d, "cohort.tsv"), "sample\tlineage\nS1\tlineage4\n")
            write(os.path.join(d, "S1.candidates.tsv"),
                  "rank\treference\tsnp_distance\n1\tR1\t5\n2\tR2\t9\n")
            out = os.path.join(d, "refmap.tsv")
            subprocess.run([py(), "bin/p1_summary.py", "--cohort",
                            os.path.join(d, "cohort.tsv"), "--dir", d,
                            "--lineages", os.path.join(d, "none.csv"),
                            "--crossmap-dir", d, "--out", out],
                           check=True, capture_output=True, text=True)
            raw = open(out, "rb").read()
            got = rows(out)[0]["reference"]
        self.assertNotIn(b"\r", raw)
        self.assertEqual(got, "R1")


if __name__ == "__main__":
    unittest.main(warnings="ignore")
