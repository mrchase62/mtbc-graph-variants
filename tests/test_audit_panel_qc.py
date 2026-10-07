#!/usr/bin/env python3
"""Regression tests for the 2026-10-05 audit, group "panel_qc".

PGB-3 (build_panel.py cannot reproduce the panel and the documented command
overwrites it), PGB-4 (the provenance rule in code is not production's),
PGB-5 (the foreign-DNA screen sums hits, misses alignment-breaking inserts,
lets two engineered genomes vouch for each other, accepts any background, and
misreads IS6110), PGB-7 (no single-base-run, long-read error or sequence-based
technology check), PGB-8 (the SNP-outlier input has no producer; MAD 0 never
flags) and PGB-13 (the RefSeq selection drops M. orygis).

Same pattern as tests/run_tests.py: small synthetic inputs, standard library
plus the pipeline's own scripts, seconds to run.

    $MTB_PY tests/test_audit_panel_qc.py
"""
import csv
import importlib.util
import io
import os
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


def rows(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def py(*args):
    return subprocess.run([sys.executable] + list(args), capture_output=True, text=True)


class BuildPanelFromRecordedDecisions(unittest.TestCase):
    """PGB-3: every decision from a recorded input; no overwrite of inputs."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        d = self.d
        T = "\t"
        # headerless selection manifest: accession, versioned, organism, ftp
        write(f"{d}/manifest.tsv", "".join(
            f"{a}{T}{a}.1{T}Mycobacterium tuberculosis{T}ftp\n" for a in
            ("GCF_000000001", "GCF_000000002", "GCF_000000003", "GCF_000000004",
             "GCF_000000005", "GCF_000195955")))
        write(f"{d}/coll.tsv", "accession\tref_covered\tbreakpoints\n" + "".join(
            f"GCF_00000000{i}\t1.0\t0\n" for i in range(1, 6)) + "GCF_000195955\t1.0\t0\n")
        write(f"{d}/summary.txt", "#assembly_accession\n" + "".join(
            f"GCF_00000000{i}.1\tPRJ\tSAMN{i}\t\t\t\t\t\tstrain=s{i}\n" for i in range(1, 6))
            + "GCF_000195955.1\tPRJ\tSAMN9\t\t\t\t\t\tstrain=H37Rv\n")
        write(f"{d}/barcode.tsv", "accession\tin_panel\tstrain\tn_lineage_roots\tmixed\troots\n"
              "GCF_000000001\tTrue\tlineage4.1\t1\tFalse\tx\n"
              "GCF_000000002\tTrue\tlineage4.1\t2\tTrue\tx\n")
        write(f"{d}/prov.tsv", "accession\tflag\n" + "".join(
            f"GCF_00000000{i}\tok\n" for i in (1, 2, 3, 5))
            + "GCF_000000004\tREFERENCE_STRUCTURE\nGCF_000195955\tREFERENCE_STRUCTURE\n")
        write(f"{d}/eng.tsv", "accession\treason\tnote\nGCF_000000005\tengineered:pJEB_at_attB\tx\n")
        write(f"{d}/ret.tsv", "accession\treason\nGCF_000195955\treference_path\n")

    def run_bp(self, *extra, out=None):
        d = self.d
        return py("bin/build_panel.py", "--manifest", f"{d}/manifest.tsv",
                  "--collinearity", f"{d}/coll.tsv", "--summary", f"{d}/summary.txt",
                  "--barcode", f"{d}/barcode.tsv", "--provenance", f"{d}/prov.tsv",
                  "--exclusions", f"{d}/eng.tsv", "--retained", f"{d}/ret.tsv",
                  "--out", out or f"{d}/panel.tsv", "--excluded", f"{d}/excl.tsv", *extra)

    def test_recorded_decisions_give_the_panel(self):
        r = self.run_bp()
        self.assertEqual(r.returncode, 0, r.stderr)
        kept = [x["accession"] for x in rows(f"{self.d}/panel.tsv")]
        self.assertEqual(kept, ["GCF_000000001", "GCF_000000003", "GCF_000195955"])
        why = {x["accession"]: x["reason"] for x in rows(f"{self.d}/excl.tsv")}
        self.assertEqual(why, {"GCF_000000002": "mixed_lineage",
                               "GCF_000000004": "provenance:REFERENCE_STRUCTURE",
                               "GCF_000000005": "engineered:pJEB_at_attB"})
        self.assertTrue(os.path.exists(f"{self.d}/panel.tsv.inputs.tsv"))

    def test_refuses_fasta_named_output(self):
        r = self.run_bp(out=f"{self.d}/mtb.complex333.fasta.gz")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("refusing", r.stderr)
        self.assertFalse(os.path.exists(f"{self.d}/mtb.complex333.fasta.gz"))

    def test_refuses_to_overwrite_an_input(self):
        before = open(f"{self.d}/prov.tsv").read()
        r = self.run_bp(out=f"{self.d}/prov.tsv")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("refusing to overwrite an input", r.stderr)
        self.assertEqual(open(f"{self.d}/prov.tsv").read(), before)

    def test_header_row_is_not_a_genome(self):
        d = self.d
        write(f"{d}/manifest_h.tsv", "accession\tx\torganism\n"
              "GCF_000000001\tx\tMycobacterium tuberculosis\n")
        r = py("bin/build_panel.py", "--manifest", f"{d}/manifest_h.tsv",
               "--collinearity", f"{d}/coll.tsv", "--summary", f"{d}/summary.txt",
               "--barcode", f"{d}/barcode.tsv", "--provenance", f"{d}/prov.tsv",
               "--out", f"{d}/p.tsv", "--excluded", f"{d}/e.tsv")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual([x["accession"] for x in rows(f"{d}/p.tsv")], ["GCF_000000001"])

    def test_barcode_is_required(self):
        d = self.d
        r = py("bin/build_panel.py", "--manifest", f"{d}/manifest.tsv",
               "--collinearity", f"{d}/coll.tsv", "--summary", f"{d}/summary.txt",
               "--provenance", f"{d}/prov.tsv", "--out", f"{d}/p.tsv",
               "--excluded", f"{d}/e.tsv")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("required: --barcode", r.stderr)

    def test_missing_provenance_row_is_not_a_pass(self):
        write(f"{self.d}/prov.tsv", "accession\tflag\nGCF_000000001\tok\n")
        r = self.run_bp()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not screened", r.stderr)

    def test_unresolved_review_flag_stops_the_build(self):
        write(f"{self.d}/runs.tsv", "accession\tflag\nGCF_000000003\tsingle_base_run\n")
        r = self.run_bp("--review", f"{self.d}/runs.tsv")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("no recorded decision", r.stderr)
        self.assertFalse(os.path.exists(f"{self.d}/panel.tsv"))


class ProvenanceRule(unittest.TestCase):
    """PGB-4: production's rule, IS6110 profile, technology from sequence."""

    def setUp(self):
        self.m = load("assembly_provenance_screen", "bin/assembly_provenance_screen.py")

    def test_rule_is_profile_then_short_read(self):
        f = self.m.provenance_flag
        self.assertEqual(f(True, "long_read"), "REFERENCE_STRUCTURE")
        self.assertEqual(f(False, "short_read"), "SHORT_READ_SUSPECT")
        self.assertEqual(f(False, "unknown"), "ok")
        self.assertEqual(f(False, "hybrid"), "ok")

    def test_displaced_flank_is_not_a_kept_copy(self):
        ref = [(100, 1455), (5000, 6355)]
        # copy 0 both flanks at reference copy 0; copy 1 one flank 617 bp off
        tot, kept, novel = self.m.classify_sites(2, {0: [100, 1455], 1: [5000, 4383]}, ref)
        self.assertEqual((tot, kept, novel), (2, 1, 1))

    def test_unplaced_copy_counts_only_in_total(self):
        self.assertEqual(self.m.classify_sites(1, {}, [(100, 1455)]), (1, 0, 0))

    def test_technology_keywords(self):
        tc = self.m.tech_class
        self.assertEqual(tc("Illumina NovaSeq; ONT"), "hybrid")
        self.assertEqual(tc("HiSeq 2500"), "short_read")
        self.assertEqual(tc("MiSeq/PacBio RSII"), "hybrid")
        self.assertEqual(tc("Ion Torrent"), "short_read")
        self.assertEqual(tc(""), "unknown")

    def test_sequence_technology_and_conflicts(self):
        m = self.m
        self.assertEqual(m.tech_from_sequence(False, 2, 5), "short_read_like")
        self.assertEqual(m.tech_from_sequence(False, 56, 5), "long_read_like")
        self.assertEqual(m.tech_from_sequence(True, 56, 5), "reference_structured")
        self.assertIn("short_read_assembler",
                      m.tech_conflicts("long_read", "long_read_like", "SOAPdenovo v. v2.04"))
        self.assertEqual(m.tech_conflicts("hybrid", "long_read_like", "Canu, SPAdes"), [])
        self.assertIn("metadata_long_read_but_sequence_short_read_like",
                      m.tech_conflicts("long_read", "short_read_like", "Flye"))


def paf(q, qlen, qs, qe, strand, t, ts, te, tp="P"):
    return (f"{q}\t{qlen}\t{qs}\t{qe}\t{strand}\t{t}\t4411532\t{ts}\t{te}\t"
            f"{qe-qs}\t{qe-qs}\t60\ttp:A:{tp}\n")


class ForeignScreen(unittest.TestCase):
    """PGB-5: gaps extracted, best single homologue, no mutual vouching,
    IS6110 re-check, PanSN one-per-sublineage background."""

    def setUp(self):
        self.m = load("foreign_insertion_screen", "bin/foreign_insertion_screen.py")

    def test_gap_between_blocks_is_extracted(self):
        # GCF_039770655's shape: blocks end at 2,001,528 and resume at 2,002,483,
        # with 1,144 bp of assembly between them
        lines = [paf("c", 4400000, 0, 2070000, "+", "h", 0, 2001528),
                 paf("c", 4400000, 2071144, 4400000, "+", "h", 2002483, 4411532),
                 paf("c", 4400000, 2070100, 2070300, "+", "h", 10, 210, tp="S")]
        g = self.m.alignment_gaps(lines, 1000, 5000)
        self.assertEqual(g, [("c", 2070000, 2071144, 2001528)])

    def test_best_single_genome_not_the_sum(self):
        lines = ["A|1|1000|indel\t1000\t0\t300\t+\tB#1#c\t1\t0\t300\t300\t300\t60\n",
                 "A|1|1000|indel\t1000\t300\t600\t+\tC#1#c\t1\t0\t300\t300\t300\t60\n",
                 "A|1|1000|indel\t1000\t0\t1000\t+\tA#1#c\t1\t0\t1000\t1000\t1000\t60\n",
                 "A|1|1000|indel\t1000\t100\t400\t+\tB#1#c\t1\t0\t300\t300\t300\t60\n"]
        cov = self.m.coverage_by_genome(lines)["A|1|1000|indel"]
        self.assertEqual(cov, {"B": 400, "C": 300})          # self excluded, merged
        best, g, n, v = self.m.score(1000, cov, 0.10, 2)
        self.assertEqual((best, g, n, v), (400, "B", 2, "native"))

    def test_two_engineered_genomes_do_not_vouch_for_each_other(self):
        # the pJEB genome's insert found only in the attB-vector genome
        best, g, n, v = self.m.score(5214, {"GCF_021535155": 3493}, 0.10, 2)
        self.assertEqual(v, "REVIEW_one_homologue")
        self.assertEqual(self.m.score(4445, {}, 0.10, 2)[3], "FOREIGN")

    def test_is6110_recheck(self):
        rv = self.m.recheck_verdict
        self.assertEqual(rv(1360, "FOREIGN", {}, 1355, 2), "native_is6110")
        self.assertEqual(rv(4445, "REVIEW_one_homologue", {"X": 3000}, 0, 2),
                         "REVIEW_one_homologue")
        self.assertEqual(rv(2000, "FOREIGN", {"X": 1500, "Y": 1200}, 0, 2), "native_recheck")

    def test_background_one_per_sublineage_and_pansn(self):
        d = tempfile.mkdtemp()
        write(f"{d}/lin.tsv", "accession\tstrain\nG1\tlineage4.1\nG2\tlineage4.1\n"
              "G3\tlineage2.2.1\nG4\t\n")
        self.assertEqual(self.m.background_from_lineages(["G4", "G3", "G2", "G1"], f"{d}/lin.tsv"),
                         ["G1", "G3", "G4"])
        write(f"{d}/bare.fa", ">NC_000962.3\nACGT\n")
        with self.assertRaises(SystemExit):
            self.m.check_background(f"{d}/bare.fa", 100)
        write(f"{d}/big.fa", "".join(f">G{i}#1#c\nACGT\n" for i in range(5)))
        with self.assertRaises(SystemExit):
            self.m.check_background(f"{d}/big.fa", 3)
        self.assertEqual(self.m.check_background(f"{d}/big.fa", 5), {f"G{i}" for i in range(5)})


BLASTN = os.environ.get(
    "MTB_BLASTN", "/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/autocycler/bin/blastn")
MAKEBLASTDB = os.path.join(os.path.dirname(BLASTN), "makeblastdb")


class UniVecCheck(unittest.TestCase):
    """The user's decision D29 (2026-10-07): every insert is also searched
    against UniVec_Core, VecScreen's way, because the homologue rule cannot
    see a construct that two genomes share."""

    def setUp(self):
        self.m = load("foreign_insertion_screen_uv", "bin/foreign_insertion_screen.py")

    def test_vecscreen_thresholds(self):
        st = self.m.vecscreen_strength
        self.assertEqual(st(30, 400, 450, 5000), "strong")      # internal
        self.assertEqual(st(29, 400, 450, 5000), "moderate")
        self.assertEqual(st(24, 1, 30, 5000), "strong")         # terminal
        self.assertEqual(st(24, 4970, 5000, 5000), "strong")
        self.assertEqual(st(24, 400, 450, 5000), "weak")        # internal weak: 23-24
        self.assertEqual(st(22, 400, 450, 5000), "")
        self.assertEqual(st(19, 26, 60, 5000), "")
        self.assertEqual(st(19, 25, 60, 5000), "moderate")

    def test_strong_hits_merged_per_insert(self):
        lines = ["I|1|5000|indel\tv1\t100\t600\t500\t5000\tvector A\n",
                 "I|1|5000|indel\tv2\t500\t900\t843\t5000\tvector B\n",
                 "I|1|5000|indel\tv3\t2000\t2030\t26\t5000\tmoderate only\n",
                 "J|1|3000|gap\tv3\t2000\t2030\t26\t3000\tmoderate only\n"]
        self.assertEqual(self.m.univec_hits(lines),
                         {"I|1|5000|indel": (801, "vector B")})   # 100..900, inclusive

    def test_vector_overrides_native(self):
        rows = {"I": dict(final_verdict="native", univec_strong_bp=0, univec_hit=""),
                "J": dict(final_verdict="native", univec_strong_bp=0, univec_hit="")}
        self.m.apply_univec(rows, {"I": (2295, "Cloning vector pDRIVE")})
        self.assertEqual(rows["I"]["final_verdict"], "VECTOR")
        self.assertEqual(rows["I"]["univec_strong_bp"], 2295)
        self.assertEqual(rows["J"]["final_verdict"], "native")

    @unittest.skipUnless(os.path.exists(BLASTN) and os.path.exists(MAKEBLASTDB),
                         "no BLAST+")
    def test_blastn_finds_an_embedded_vector(self):
        import random
        rnd = random.Random(7)
        seq = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))  # noqa: E731
        vec = seq(400)
        with tempfile.TemporaryDirectory() as d:
            write(f"{d}/uv.fa", f">gnl|uv|TEST:1-400 test vector\n{vec}\n")
            write(f"{d}/ins.fa", f">A|1|3000|indel\n{seq(1300)}{vec}{seq(1300)}\n"
                                 f">B|1|3000|gap\n{seq(3000)}\n")
            subprocess.run([MAKEBLASTDB, "-in", f"{d}/uv.fa", "-dbtype", "nucl",
                            "-out", f"{d}/uv"], check=True, capture_output=True)
            out = subprocess.run([BLASTN] + self.m.VECSCREEN
                                 + ["-db", f"{d}/uv", "-query", f"{d}/ins.fa",
                                    "-outfmt", self.m.BLAST_FMT],
                                 check=True, capture_output=True, text=True).stdout
        hits = self.m.univec_hits(out.splitlines())
        self.assertEqual(list(hits), ["A|1|3000|indel"])
        self.assertEqual(hits["A|1|3000|indel"][0], 400)

    def test_failed_alignment_is_fatal_not_empty(self):
        with tempfile.TemporaryDirectory() as d:
            mm2 = write(f"{d}/mm2", "#!/bin/sh\necho 'no such ref' >&2\nexit 1\n")
            os.chmod(mm2, 0o755)
            with self.assertRaises(SystemExit) as e:
                self.m.extract_inserts("G", "a.fa", "missing.fa", mm2, "k8", "p.js",
                                       1000, 5000, io.StringIO())
        self.assertIn("minimap2 failed on G", str(e.exception))


class SingleBaseRuns(unittest.TestCase):
    """PGB-7a: a pure run of >= 100 bp is flagged, 99 is not."""

    def test_run_length_flag(self):
        d = tempfile.mkdtemp()
        body = "ACGT" * 300
        write(f"{d}/A.fasta", f">a\n{body}{'G' * 100}{body}\n")
        write(f"{d}/B.fasta", f">b\n{body}{'G' * 99}{body}\n")
        r = py("bin/assembly_qc_stats.py", "--dir", d, "--suffix", ".fasta",
               "--out", f"{d}/s.tsv")
        self.assertEqual(r.returncode, 0, r.stderr)
        got = {x["accession"]: x for x in rows(f"{d}/s.tsv")}
        self.assertEqual(got["A"]["flag"], "single_base_run")
        self.assertEqual(got["A"]["longest_base_run"], "100")
        self.assertEqual(got["A"]["longest_base_run_pos"], "1201")
        self.assertEqual(got["B"]["flag"], "")


class LongReadCounts(unittest.TestCase):
    """PGB-7e / PGB-8: the counts' producer is in bin/."""

    def setUp(self):
        self.m = load("variant_counts", "bin/variant_counts.py")

    def test_homopolymer_indels_and_privacy(self):
        h37 = "ACGTAAAAACGTCCCCGT"
        # deletion of an A inside AAAAA (0-based 5); insertion of C beside CCCC;
        # a 1 bp indel outside a run; a SNP
        var = [(5, 6, "A", "-"), (12, 12, "-", "C"), (1, 2, "C", "-"), (2, 3, "G", "T")]
        snps, i1, hp, i2, keys, _ = self.m.count_variants(var, h37)
        self.assertEqual((snps, i1, hp, i2), (1, 3, 2, 0))
        priv = self.m.private_counts({"x": keys, "y": {(5, "A", "-")}})
        self.assertEqual(priv, {"x": 1, "y": 0})


class OutlierScreen(unittest.TestCase):
    """PGB-8: every candidate must have a count; MAD 0 does not mean pass."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        lin = ["accession\tstrain"] + [f"S{i}\tlineage4.1.1" for i in range(6)] + \
              [f"T{i}\tlineage4.2.1" for i in range(6)]
        write(f"{self.d}/lin.tsv", "\n".join(lin) + "\n")

    def test_missing_count_stops(self):
        write(f"{self.d}/c.tsv", "S0\t5\n")
        write(f"{self.d}/acc.txt", "S0\nS1\n")
        r = py("bin/snp_outlier_screen.py", "--counts", f"{self.d}/c.tsv",
               "--lineages", f"{self.d}/lin.tsv", "--accessions", f"{self.d}/acc.txt",
               "--out", f"{self.d}/o.tsv")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not screened", r.stderr)

    def test_mad_zero_falls_back_to_major(self):
        # sublineage 4.1.1: five zeros and one 40 -> MAD 0; major lineage4 has spread
        c = [f"S{i}\t0" for i in range(5)] + ["S5\t40"] + \
            [f"T{i}\t{v}" for i, v in enumerate((1, 2, 3, 2, 1, 2))]
        write(f"{self.d}/c.tsv", "\n".join(c) + "\n")
        r = py("bin/snp_outlier_screen.py", "--counts", f"{self.d}/c.tsv",
               "--lineages", f"{self.d}/lin.tsv", "--column", "indel1_hp_private",
               "--out", f"{self.d}/o.tsv")
        self.assertEqual(r.returncode, 0, r.stderr)
        got = {x["accession"]: x for x in rows(f"{self.d}/o.tsv")}
        self.assertEqual(got["S5"]["flag"], "EXCESS")
        self.assertIn("major lineage4", got["S5"]["group"])
        self.assertIn("indel1_hp_private", got["S5"])


class RefseqSelection(unittest.TestCase):
    """PGB-13: M. orygis is selected; so are canetti and 'complex sp.'."""

    def test_orygis_selected(self):
        d = tempfile.mkdtemp()
        write(f"{d}/env.sh", f"export MTB_DATA={d}/data\n")

        def row(acc, sp, org):
            f = [""] * 20
            f[0], f[6], f[7], f[10], f[11], f[19] = acc, sp, org, "latest", "Complete Genome", "ftp"
            return "\t".join(f) + "\n"
        with open(f"{d}/sum.txt", "w") as fh:
            fh.write("#assembly_accession\n")
            fh.write(row("GCF_000000001.1", "1773", "Mycobacterium tuberculosis H37Rv"))
            fh.write(row("GCF_015265495.1", "1305738", "Mycobacterium orygis"))
            fh.write(row("GCF_035581225.1", "78331", "Mycobacterium canetti"))
            fh.write(row("GCF_030323705.1", "2583589", "Mycobacterium tuberculosis complex sp. N0052"))
            fh.write(row("GCF_000000009.1", "1764", "Mycobacterium avium"))
        r = subprocess.run(["bash", "bin/refresh_assemblies.sh", "--summary-only",
                            "--summary", f"{d}/sum.txt"],
                           capture_output=True, text=True,
                           env=dict(os.environ, MTB_ENV_FILE=f"{d}/env.sh"))
        self.assertEqual(r.returncode, 0, r.stderr)
        sel = [f for f in os.listdir(f"{d}/data/ncbi") if f.startswith("selected.")]
        got = sorted(l.split("\t")[0] for l in open(f"{d}/data/ncbi/{sel[0]}"))
        self.assertEqual(got, ["GCF_000000001", "GCF_015265495", "GCF_030323705",
                               "GCF_035581225"])


class PansnRefusesToOverwrite(unittest.TestCase):
    """PGB-3: make_pansn_fasta.sh never overwrites a built panel."""

    def test_existing_panel_not_overwritten(self):
        d = tempfile.mkdtemp()
        os.makedirs(f"{d}/fastas"); os.makedirs(f"{d}/bin")
        write(f"{d}/bin/bgzip", "#!/bin/sh\nexit 0\n")
        write(f"{d}/env.sh", f"""
            export MTB_DATA={d} MTB_FASTAS={d}/fastas MTB_QC_BIN={d}/bin
            mtb_require_file() {{ [[ -e "$1" ]] || exit 1; }}
            """)
        write(f"{d}/fastas/mtb.complex333.fasta.gz", "precious\n")
        write(f"{d}/panel.tsv", "accession\torganism\tstrain\n")
        env = dict(os.environ, MTB_ENV_FILE=f"{d}/env.sh")
        r = subprocess.run(["bash", "bin/make_pansn_fasta.sh", "--panel", f"{d}/panel.tsv",
                            "--out", "mtb.complex333"], capture_output=True, text=True, env=env)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("refusing to overwrite", r.stderr)
        self.assertEqual(open(f"{d}/fastas/mtb.complex333.fasta.gz").read(), "precious\n")
        r = subprocess.run(["bash", "bin/make_pansn_fasta.sh"], capture_output=True,
                           text=True, env=env)
        self.assertNotEqual(r.returncode, 0)        # no stale default panel
        self.assertIn("usage", r.stderr)


if __name__ == "__main__":
    unittest.main(warnings="ignore")
