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

    def select(self, panel_records, isolate_records):
        with tempfile.TemporaryDirectory() as d:
            panel = vcf(os.path.join(d, "panel.vcf"), ["gA", "gG"], panel_records)
            iso = os.path.join(d, "iso.vcf")
            with open(iso, "w") as fh:
                fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
                for pos, ref, alt, gt in isolate_records:
                    fh.write(f"c\t{pos}\t.\t{ref}\t{alt}\t50\t.\t.\tGT:DP\t{gt}:30\n")
            out = os.path.join(d, "cand.tsv")
            subprocess.run([py(), "bin/t8_select_reference.py", "--vcf", iso,
                            "--panel-snps", panel, "--out", out, "--top", "2"],
                           check=True, capture_output=True, text=True)
            return {r["reference"]: int(r["snp_distance"]) for r in rows(out)}

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


# ------------------------------------------------------ GRAPHVCF-5, PGB-9

def collapse_py():
    """The collapse program embedded in bin/vcf_collapse.sh."""
    s = open("bin/vcf_collapse.sh").read()
    m = re.search(r"<<'PY' \|\| true\n(.*?)\nPY\n", s, re.S)
    if not m:
        raise AssertionError("no COLLAPSE_PY program in bin/vcf_collapse.sh")
    return m.group(1)


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

    @unittest.skipUnless(
        shutil.which("bcftools") or (
            shutil.which("singularity")
            and os.path.exists(os.environ.get("MTB_PGGB_SIF", "")
                               or "/n/netscratch/sfortune_lab/Lab/mchase/MtbPangenome/containers/pggb_latest.sif")),
        "no bcftools (host or container)")
    def test_collapse_script_end_to_end(self):
        # the whole script, as the build runs it: no duplicate key, no padded
        # SNP, every carrier kept, AC refilled
        with tempfile.TemporaryDirectory(dir="/tmp") as d:
            src = vcf(os.path.join(d, "in.vcf"), ["s0", "s1", "s2"], self.RECORDS)
            out = os.path.join(d, "out.vcf.gz")
            subprocess.run(["bash", "bin/vcf_collapse.sh", "TEST", "--in", src,
                            "--out", out], check=True, capture_output=True,
                           text=True)
            text = gzip.open(out, "rt").read()
        got = parse_vcf_text(text)
        self.assertEqual(got, self.WANT)
        info = [l.split("\t")[7] for l in text.splitlines()
                if l.startswith("c\t10\t")][0]
        self.assertIn("AC=2", info.split(";"))


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
