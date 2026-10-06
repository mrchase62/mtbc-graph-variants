#!/usr/bin/env python3
"""Regression tests for the second clean-up group (2026-10-06): a new build
runs end to end with every input from the build or an explicit argument, and
every output is stamped with the new build's id.

Same pattern as tests/run_tests.py: unittest, small synthetic inputs, the
standard library plus the pipeline's own scripts, seconds to run.

    $MTB_PY tests/test_audit_cleanup2.py -v
"""
import csv
import gzip
import hashlib
import importlib.util
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

H = "GCF_000195955#1#NC_000962.3"
H37 = "GCF_000195955"


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


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def env(**kw):
    e = dict(os.environ)
    e.update(MTB_SITE_FILE="/dev/null", MTB_PY=sys.executable)
    for k in ("MTB_BUILD_DIR", "OG", "GRAPH_DIR", "BUILD_ID", "SLURM_ARRAY_TASK_ID",
              "REFMAP", "MTB_GRAPH_FRAMES", "ACC_CATALOGUE", "ANC_TREE", "ANC_ALN",
              "ANC_SITES", "ANC_OUTGROUPS", "ACCESSORY_DIR", "OUTGROUP",
              "MTB_OUTGROUP", "IS6110_LOCI", "IS6110_ANCHORS", "BUILD_ROOT",
              "SAMPLES", "QUERY", "OUT", "ACCPRES", "LINTAB", "REFS"):
        e.pop(k, None)
    e.update({k: str(v) for k, v in kw.items()})
    return e


def run(cmd, cwd=None, **kw):
    return subprocess.run(cmd, cwd=cwd, env=env(**kw), capture_output=True,
                          text=True)


def bcftools():
    b = os.environ.get("MTB_BCFTOOLS") or shutil.which("bcftools")
    if not b:
        c = "/n/boslfs02/LABS/sfortune_lab/Lab/conda/envs/mtb_isolates/bin/bcftools"
        b = c if os.path.exists(c) else ""
    return b


def vcf_gz(path, samples):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n##contig=<ID=c>\n#CHROM\tPOS\tID\tREF\t"
                 "ALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
    return path


def logger(path, log):
    """An executable that appends its argv to `log` and exits 0."""
    return write(path, f'#!/bin/bash\necho "$@" >> {log}\n', mode=0o755)


def calls(log, needle):
    if not os.path.exists(log):
        return []
    return [l for l in open(log).read().splitlines() if needle in l]


# ==========================================================================
# 1. The association tail: node -> locus table, and the build's outgroup
# ==========================================================================
class AssocTail(unittest.TestCase):
    """cohort_assoc_tail.sh FATALed whenever presence tables existed (the
    event writer needs --node-locus and it passed none), and rooted every
    tree and named every outgroup GCF_035581225 whatever the build."""

    CHAIN = os.path.join(ROOT, "assoc", "bin", "cohort_assoc_tail.sh")

    def setUp(self):
        self.d = d = tempfile.mkdtemp()
        os.makedirs(f"{d}/assoc")
        os.symlink(os.path.join(ROOT, "assoc", "bin"), f"{d}/assoc/bin")
        os.symlink(os.path.join(ROOT, "bin"), f"{d}/bin")
        with gzip.open(f"{d}/refbias/t/p5/merged.vcf.gz" if os.makedirs(
                f"{d}/refbias/t/p5") is None else "", "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n##MTB_graph_build=b0\n"
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        self.b = f"{d}/refbias/build/b0"
        write(f"{self.b}/build_info.tsv", "build_id\tb0\noutgroup\tGCF_X\n")
        write(f"{self.b}/assets/panel_polarity.tsv", "chrom\n")
        write(f"{self.b}/assets/graph_collapsed.vcf.gz", "x\n")
        write(f"{d}/accessory/t/S1.presence.tsv", "locus\n")
        write(f"{d}/refbias/t/p4/S1.placed.tsv", "frame\n")
        os.makedirs(f"{d}/data/trees")
        os.makedirs(f"{d}/assoc/t/events")
        self.log = f"{d}/calls.log"
        self.stub = logger(f"{d}/stub.sh", self.log)
        os.makedirs(f"{d}/fakebin")
        logger(f"{d}/fakebin/sbatch", self.log)

    def tearDown(self):
        shutil.rmtree(self.d)

    def tail(self, **kw):
        e = env(MTB_PY=self.stub, MTB_PY_VT=self.stub,
                PATH=f"{self.d}/fakebin:" + os.environ["PATH"], **kw)
        return subprocess.run(["bash", self.CHAIN, "t"], cwd=self.d, env=e,
                              capture_output=True, text=True)

    def prov(self, rel):
        # the record the chain itself expects (review 2, R2-TREES-2: the
        # inputs and code of each product, not only the build and VCF)
        write(f"{self.d}/{rel}", "x\n")
        r = self.tail(MTB_CHAIN_PRINT_PROV=rel)
        assert r.returncode == 0, r.stderr + r.stdout
        rec = "".join(l + "\n" for l in r.stdout.splitlines() if "\t" in l)
        write(f"{self.d}/{rel}.prov", rec)

    def test_tree_rooted_on_the_builds_outgroup(self):
        r = self.tail()
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        og = calls(self.log, "combined_alignment.py")
        self.assertTrue(og and "--outgroup GCF_X" in og[0], og)
        tree = calls(self.log, "build_snp_tree.sh")
        self.assertTrue(tree, open(self.log).read())
        # review 2, R2-TREES-1: the tree is the cohort + panel tree
        self.assertIn("data/trees/t.combined.fasta GCF_X data/trees/t.combined",
                      tree[0])
        self.assertNotIn("GCF_035581225", code("assoc/bin/cohort_assoc_tail.sh"))

    def test_replaced_tree_makes_the_event_matrix_stale(self):
        """Review 2, R2-TREES-2: a record holds the code and every input, so
        a tree replaced after the event matrix was made is refused there,
        naming the input that changed."""
        self.prov("data/trees/t.combined.rooted.nwk")
        self.prov("data/trees/t.rooted.nwk")
        self.prov("assoc/t/node_locus.tsv")
        self.prov("assoc/t/events/summary.txt")
        rec = open(f"{self.d}/assoc/t/events/summary.txt.prov").read()
        for key in ("code", "tree", "outgroup_aln", "polarity", "outgroup",
                    "presence", "node_locus"):
            self.assertIn(f"{key}\t", rec)
        # the pruned tree is rebuilt (its own record kept current)
        self.prov("data/trees/t.rooted.nwk")
        write(f"{self.d}/data/trees/t.rooted.nwk", "(changed);\n")
        self.prov_only("data/trees/t.rooted.nwk")
        r = self.tail()
        self.assertNotEqual(r.returncode, 0, r.stdout)
        self.assertIn("events/summary.txt exists but was not made from this "
                      "run's inputs", r.stderr)
        self.assertIn("tree", r.stderr.split("differs in:")[1])

    def prov_only(self, rel):
        r = self.tail(MTB_CHAIN_PRINT_PROV=rel)
        assert r.returncode == 0, r.stderr + r.stdout
        write(f"{self.d}/{rel}.prov",
              "".join(l + "\n" for l in r.stdout.splitlines() if "\t" in l))

    def test_presence_tables_get_the_cohorts_node_locus_table(self):
        self.prov("data/trees/t.combined.rooted.nwk")
        self.prov("data/trees/t.rooted.nwk")
        r = self.tail()
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        nl = calls(self.log, "node_locus_from_p4.py")
        self.assertTrue(nl, open(self.log).read())
        self.assertIn("--p4dir refbias/t/p4 --out assoc/t/node_locus.tsv", nl[0])
        ev = calls(self.log, "write_event_matrix.py")
        self.assertTrue(ev)
        self.assertIn("--accessory-presence accessory/t "
                      "--node-locus assoc/t/node_locus.tsv", ev[0])
        self.assertIn("--outgroup-name GCF_X", ev[0])
        # the node table is made before the event matrix that reads it
        lines = open(self.log).read().splitlines()
        self.assertLess(lines.index(nl[0]), lines.index(ev[0]))

    def test_build_without_outgroup_key_uses_config(self):
        write(f"{self.b}/build_info.tsv", "build_id\tb0\n")
        r = self.tail(MTB_OUTGROUP="GCF_Y")
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        self.assertIn("--outgroup GCF_Y",
                      calls(self.log, "combined_alignment.py")[0])

    def test_build_with_no_outgroup_refused(self):
        write(f"{self.b}/build_info.tsv", "build_id\tb0\noutgroup\tnone\n")
        r = self.tail()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("has no outgroup", r.stderr)


# ==========================================================================
# 2. The outgroup in one place
# ==========================================================================
class OutgroupConfig(unittest.TestCase):
    """GCF_035581225 was written separately into P0, the tail, add_outgroup,
    panel_polarity and build_snp_tree."""

    def test_config_default_unchanged(self):
        r = subprocess.run(
            ["bash", "-c", f"source {ROOT}/config/project_env.sh; "
             "printf '%s|' \"${MTB_OUTGROUP-unset}\"; "
             "MTB_OUTGROUP= bash -c 'source " + ROOT +
             "/config/project_env.sh; printf \"%s|\" \"${MTB_OUTGROUP-unset}\"'"],
            env=env(), capture_output=True, text=True)
        # default the CX333 canettii; set empty, stays empty (no outgroup)
        self.assertEqual(r.stdout, "GCF_035581225||")

    def test_scripts_take_it_from_config_or_argument(self):
        self.assertNotIn("GCF_035581225", code("bin/build_snp_tree.sh"))
        self.assertIn("MTB_OUTGROUP", code("bin/build_snp_tree.sh"))
        self.assertNotIn("GCF_035581225", code("bin/p0_prepare.sh"))
        for p in ("assoc/bin/add_outgroup.py", "bin/panel_polarity.py"):
            self.assertNotRegex(text(p), r'default="GCF_035581225"', p)
        r = run([sys.executable, "bin/panel_polarity.py", "--panel-vcf", "x",
                 "--out", "o"])
        self.assertEqual(r.returncode, 2)
        self.assertIn("--outgroup", r.stderr)


# ==========================================================================
# 3. P0: no CX333 graph, accessory panel checked, no refbias/t11 defaults
# ==========================================================================
class P0Inputs(unittest.TestCase):

    def setUp(self):
        self.d = d = tempfile.mkdtemp()
        os.symlink(os.path.join(ROOT, "bin"), f"{d}/bin")
        self.og = write(f"{d}/g/g.smooth.final.og", "graph bytes\n")
        gsha = sha(self.og)
        self.bid = gsha[:12]
        self.b = f"{d}/build/{self.bid}"
        write(f"{self.b}/build_info.tsv",
              f"build_id\t{self.bid}\ngraph_sha256\t{gsha}\n")
        write(f"{self.b}/assets/accessions.txt", f"A1\nA2\n{H37}\n")
        write(f"{self.b}/logs/accessions.done", "2026-10-06\n")

    def tearDown(self):
        shutil.rmtree(self.d)

    def p0(self, *args, **kw):
        kw.setdefault("BUILD_ROOT", f"{self.d}/build")
        kw.setdefault("MTB_GATK_SIF", "/nonexistent.sif")
        return run(["bash", os.path.join(ROOT, "bin/p0_prepare.sh"), *args],
                   cwd=self.d, **kw)

    def test_no_graph_named_refused(self):
        # a CX333 graph where the old default looked: not taken
        write(f"{self.d}/graphs/CX333.s10k.k23.K15/x.smooth.final.og", "old\n")
        r = self.p0("--list")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("set OG", r.stderr)
        self.assertNotIn("graphs/CX333", code("bin/p0_prepare.sh"))

    def test_stamp_records_the_outgroup(self):
        shutil.rmtree(f"{self.d}/build")
        r = self.p0("--step", "stamp", OG=self.og, MTB_OUTGROUP="GCF_Q")
        self.assertEqual(r.returncode, 0, r.stderr)
        info = dict(l.rstrip("\n").split("\t", 1)
                    for l in open(f"{self.b}/build_info.tsv"))
        self.assertEqual(info["outgroup"], "GCF_Q")

    def test_polarity_with_another_outgroup_refused(self):
        with open(f"{self.b}/build_info.tsv", "a") as fh:
            fh.write("outgroup\tGCF_Q\n")
        write(f"{self.b}/logs/assets.done", "2026-10-06\n")
        r = self.p0("--step", "panel_polarity", OG=self.og, OUTGROUP="GCF_Z")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("stamped with outgroup GCF_Q", r.stderr)

    def assets(self, accdir, **kw):
        vcf_gz(f"{self.d}/g/snps.vcf.gz", ["A1", "A2"])
        vcf_gz(f"{self.d}/g/all_variants.collapsed.vcf.gz", ["A1", "A2"])
        mask = write(f"{self.d}/mask.bed", "c\t0\t5\n")
        kw.setdefault("MTB_BCFTOOLS", bcftools())
        if accdir is not None:
            kw["ACCESSORY_DIR"] = accdir
        return self.p0("--step", "assets", OG=self.og, MASK=mask, **kw)

    def panel(self, genomes):
        a = f"{self.d}/acc"
        write(f"{a}/panel_manifest.tsv", "locus_id\tpos\nACC_0000100\t100\n")
        write(f"{a}/accessory_novel.fasta", ">ACC_0000100\nACGT\n")
        if genomes is not None:
            write(f"{a}/panel_manifest.genomes.txt", "".join(g + "\n" for g in genomes))
        return a

    def test_assets_needs_an_accessory_dir(self):
        r = self.assets(None)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("needs ACCESSORY_DIR", r.stderr)

    @unittest.skipUnless(bcftools(), "bcftools not available")
    def test_panel_without_genome_record_refused(self):
        r = self.assets(self.panel(None))
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("panel_manifest.genomes.txt is missing", r.stderr)
        self.assertFalse(os.path.exists(f"{self.b}/logs/assets.done"))

    @unittest.skipUnless(bcftools(), "bcftools not available")
    def test_panel_of_other_genomes_refused(self):
        r = self.assets(self.panel(["A1", "Z9"]))
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not in the build", r.stderr)
        self.assertFalse(os.path.exists(f"{self.b}/logs/assets.done"))

    @unittest.skipUnless(bcftools(), "bcftools not available")
    def test_matching_panel_copied_and_no_t11_default(self):
        # the pilot's IS6110 files where the old default looked: not copied
        write(f"{self.d}/refbias/t11/loci_clean.tsv", "locus\tpos\n")
        write(f"{self.d}/refbias/t11/anchor_sets.tsv", "locus\tpos\n")
        a = self.panel(["A1", "A2"])
        r = self.assets(a)
        self.assertEqual(r.returncode, 0, r.stderr)
        x = f"{self.b}/assets"
        self.assertEqual(sha(f"{x}/accessory_loci.tsv"), sha(f"{a}/panel_manifest.tsv"))
        self.assertEqual(open(f"{x}/accessory_loci.genomes.txt").read(), "A1\nA2\n")
        self.assertTrue(os.path.exists(f"{x}/accessory_novel.fasta"))
        self.assertFalse(os.path.exists(f"{x}/is6110_loci.tsv"))
        self.assertFalse(os.path.exists(f"{x}/is6110_anchors.tsv"))
        self.assertNotIn("refbias/t11", code("bin/p0_prepare.sh"))
        self.assertNotIn("refbias/panel", code("bin/p0_prepare.sh"))


class P0CheckAccessoryPanel(unittest.TestCase):

    def test_record_against_build(self):
        with tempfile.TemporaryDirectory() as d:
            acc = write(f"{d}/acc.txt", f"A1\nA2\n{H37}\n")
            write(f"{d}/p/panel_manifest.genomes.txt", "A1\nA2\n")
            write(f"{d}/q/panel_manifest.genomes.txt", "A1\n")
            ok = run([sys.executable, "bin/p0_check.py", "accessory-panel",
                      "--dir", f"{d}/p", "--accessions", acc, "--ref-acc", H37])
            bad = run([sys.executable, "bin/p0_check.py", "accessory-panel",
                       "--dir", f"{d}/q", "--accessions", acc, "--ref-acc", H37])
            none = run([sys.executable, "bin/p0_check.py", "accessory-panel",
                        "--dir", f"{d}/r", "--accessions", acc, "--ref-acc", H37])
        self.assertEqual(ok.returncode, 0, ok.stderr)
        self.assertEqual(bad.returncode, 1)
        self.assertEqual(none.returncode, 1)
        self.assertIn("is missing", none.stderr)


# ==========================================================================
# 4. The accessory panel's producers, in the repository
# ==========================================================================
class AccessoryPanelProducers(unittest.TestCase):
    """refbias/panel/ had no producer in the repository; the producers lived
    in the working tree with the pilot's paths and a fixed 332 genomes."""

    def setUp(self):
        self.d = d = tempfile.mkdtemp()
        write(f"{d}/cand.tsv", "allele_id\tpos\tlen\tAC\tcarriers\n"
              "A00001_100_4\t100\t4\t2\tG1,G2\n")
        write(f"{d}/cand.fasta", ">A00001_100_4\nACGT\n")
        write(f"{d}/val.tsv", "allele_id\tpos\tlen\tAC\tcarrier_median_ident\t"
              "carrier_frac_fulllength\tnoncarrier_median_ident\t"
              "noncarrier_frac_fulllength\tseparation\n"
              "A00001_100_4\t100\t4\t2\t100\t1.0\t0\t0\t1.0\n")
        self.blast = write(f"{d}/blastn", "#!/bin/bash\nexit 0\n", mode=0o755)

    def tearDown(self):
        shutil.rmtree(self.d)

    def build(self, genomes):
        g = write(f"{self.d}/genomes.txt", "".join(x + "\n" for x in genomes))
        return run([sys.executable, os.path.join(ROOT, "bin/build_accessory_panel.py"),
                    "--validation", f"{self.d}/val.tsv",
                    "--candidates", f"{self.d}/cand.tsv",
                    "--fasta", f"{self.d}/cand.fasta", "--genomes", g,
                    "--h37rv", "h", "--blastn", self.blast,
                    "--outdir", f"{self.d}/panel"], cwd=self.d)

    def test_carrier_frac_over_the_graphs_genomes_and_recorded(self):
        r = self.build(["G1", "G2", "G3", "G4"])
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        m = rows(f"{self.d}/panel/panel_manifest.tsv")
        self.assertEqual(m[0]["carrier_frac"], "0.5")          # was 2/332
        self.assertEqual(open(f"{self.d}/panel/panel_manifest.genomes.txt").read(),
                         "G1\nG2\nG3\nG4\n")

    def test_carriers_outside_the_genomes_refused(self):
        r = self.build(["G1", "G3"])
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not in", r.stderr)

    def test_no_pilot_defaults(self):
        for p in ("bin/build_accessory_panel.py", "bin/t2_extract_candidates.py",
                  "bin/t2_summary.py", "bin/t2_validate_accessory.sh"):
            c = code(p)
            for pat in ("refbias/t2", "refbias/panel", "refbias/T2", "hely_tatc"):
                self.assertNotIn(pat, c, (p, pat))
        r = run([sys.executable, os.path.join(ROOT, "bin/t2_extract_candidates.py"),
                 "--vcf", "v", "--bcftools", "b"], cwd=self.d)
        self.assertEqual(r.returncode, 2)
        self.assertIn("--out-fasta", r.stderr)
        r = run(["bash", os.path.join(ROOT, "bin/t2_validate_accessory.sh"), "1"],
                cwd=self.d)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("SAMPLES is unset", r.stderr)


# ==========================================================================
# 5. IS6110: p1iv, stage 1 of p1is, and stage-1 discovery from the build
# ==========================================================================
class P1ivFromTheBuild(unittest.TestCase):
    """bin/p1i_vcf.sh fell back to build 7713a8d71d8e and called the writer
    without --refs, --h37rv or --build-id, so every build's IS6110 records
    were stamped 7713a8d71d8e."""

    def setUp(self):
        self.d = d = tempfile.mkdtemp()
        self.og = write(f"{d}/g.og", "graph\n")
        self.b = f"{d}/build/B1"
        write(f"{self.b}/build_info.tsv", f"build_id\tB1\ngraph\t{self.og}\n")
        write(f"{self.b}/logs/manifest.done", "x\n")
        write(f"{self.b}/refs/{H37}.fasta", ">c\nACGT\n")
        write(f"{d}/coh/p1i/S1.elstacks.tsv", "x\n")
        self.log = f"{d}/py.log"
        self.py = logger(f"{d}/py", self.log)

    def tearDown(self):
        shutil.rmtree(self.d)

    def p1iv(self, **kw):
        return run(["bash", os.path.join(ROOT, "bin/p1i_vcf.sh")], cwd=self.d,
                   SLURM_SUBMIT_DIR=self.d, MTB_PY=self.py, P1IDIR=f"{self.d}/coh/p1i",
                   REFMAP=f"{self.d}/refmap.tsv", RES=f"{self.d}/res", **kw)

    def test_writer_gets_the_builds_refs_h37rv_and_id(self):
        self.p1iv(BUILD_ROOT=f"{self.d}/build")
        w = calls(self.log, "is6110_write_vcf.py")
        self.assertTrue(w, open(self.log).read() if os.path.exists(self.log) else "")
        self.assertIn("--build-id B1", w[0])
        self.assertIn(f"--refs {self.b}/refs", w[0])
        self.assertIn(f"--h37rv {self.b}/refs/{H37}.fasta", w[0])
        f = calls(self.log, "is6110_place_by_flank.py")
        self.assertIn(f"--refs {self.b}/refs", f[0])
        p = calls(self.log, "is6110_project_sites.py")
        self.assertIn(f"--graph {self.og}", p[0])
        self.assertIn(f"--node-lengths {self.b}/assets/node_positions.tsv", p[0])

    def test_no_build_refused(self):
        r = self.p1iv(BUILD_ROOT=f"{self.d}/nobuild")
        self.assertNotEqual(r.returncode, 0)
        self.assertEqual(calls(self.log, "is6110_write_vcf.py"), [])
        self.assertNotIn("7713a8d71d8e", code("bin/p1i_vcf.sh"))


class P1isStage1H37Rv(unittest.TestCase):
    """Stage 1 of p1is called is6110_p5_merge.py without --h37rv, whose
    default was build 7713a8d71d8e's H37Rv."""

    def test_stage1_gets_the_builds_h37rv(self):
        with tempfile.TemporaryDirectory() as d:
            b = f"{d}/build/B1"
            write(f"{b}/refs/{H37}.fasta", ">c\nACGT\n")
            write(f"{d}/res/t_p1i_cohort_keys.tsv", "x\n")
            log = f"{d}/py.log"
            r = run(["bash", os.path.join(ROOT, "bin/p1i_p5states.sh"), "--stage1"],
                    cwd=d, SLURM_SUBMIT_DIR=d, MTB_PY=logger(f"{d}/py", log),
                    MTB_BUILD_DIR=b, COHORT_NAME="t", RES=f"{d}/res",
                    REFMAP=f"{d}/refmap.tsv", P1IDIR=f"{d}/p1i")
            self.assertEqual(r.returncode, 0, r.stderr)
            m = calls(log, "is6110_p5_merge.py")
        self.assertIn(f"--h37rv {b}/refs/{H37}.fasta", m[0])


class DiscoverEveryPanelGenome(unittest.TestCase):
    """is6110/bin/p1i_discover_matched.sh read CX333's refs/ and the pilot's
    refmap, and kept any GFF that existed whatever sequence it was found in."""

    def test_every_accession_from_the_build_redone_on_changed_bytes(self):
        with tempfile.TemporaryDirectory() as d:
            b = f"{d}/build/B1"
            write(f"{b}/build_info.tsv", "build_id\tB1\n")
            write(f"{b}/assets/accessions.txt", "A1\nA2\nA3\n")
            for a in ("A1", "A2", "A3"):
                write(f"{b}/refs/{a}.fasta", f">{a}\nACGT{a}\n")
            g = f"{d}/gff"
            # A1: found in these bytes; A2: in other bytes; A3: no GFF
            write(f"{g}/A1.is6110.gff", "##gff-version 3\n")
            write(f"{g}/A1.ref.sha256", sha(f"{b}/refs/A1.fasta") + "\n")
            write(f"{g}/A2.is6110.gff", "##gff-version 3\n")
            write(f"{g}/A2.ref.sha256", "0" * 64 + "\n")
            log = f"{d}/py.log"
            fake = write(f"{d}/py", f"""\
                #!/bin/bash
                echo "$@" >> {log}
                while [[ $# -gt 0 ]]; do
                  [[ $1 == --out-gff ]] && printf '##gff-version 3\\n' > "$2"
                  shift
                done
                """, mode=0o755)
            r = run(["bash", os.path.join(ROOT, "is6110/bin/p1i_discover_matched.sh")],
                    cwd=d, SLURM_SUBMIT_DIR=d, MTB_BUILD_DIR=b, MTB_PY=fake,
                    OUTDIR=g)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            got = open(log).read() if os.path.exists(log) else ""
            sha2 = open(f"{g}/A2.ref.sha256").read().strip()
        self.assertNotIn("A1.fasta", got)                    # same bytes: kept
        self.assertIn(f"{b}/refs/A2.fasta", got)             # other bytes: redone
        self.assertIn(f"{b}/refs/A3.fasta", got)             # every panel genome
        self.assertEqual(sha2, hashlib.sha256(b">A2\nACGTA2\n").hexdigest())
        self.assertNotIn("7713a8d71d8e", code("is6110/bin/p1i_discover_matched.sh"))


# ==========================================================================
# 6. No build, graph or pilot default in a script a production run reaches
# ==========================================================================
PRODUCTION_PY = [
    "is6110/bin/is6110_write_vcf.py", "is6110/bin/is6110_p5_merge.py",
    "is6110/bin/is6110_p5_stage2.py", "accessory/bin/merge_catalogues.py",
    "accessory/bin/locus_presence.py", "bin/merge_cohort_vcf.py",
    "bin/sv_intervals.py", "bin/t8_select_reference.py",
    "graphframe/bin/graph_frame_offsets.py", "graphframe/bin/frame_detect.py",
    "graphframe/bin/graph_frame.py", "bin/panel_polarity.py",
    "assoc/bin/add_outgroup.py", "bin/retier_intervals.py",
    "refbias/bin/io_contract.py", "bin/build_accessory_panel.py",
    "bin/t2_extract_candidates.py", "bin/t2_summary.py"]
PRODUCTION_SH = [
    "bin/p0_prepare.sh", "bin/p1i_vcf.sh", "bin/p1i_p5states.sh",
    "bin/p3_accessory.sh", "is6110/bin/p1i_discover_matched.sh",
    "is6110/bin/p1i_build_matched.sh", "assoc/bin/cohort_assoc_tail.sh",
    "bin/refbias_run.sh", "bin/t2_validate_accessory.sh"]
OLD = re.compile(r"7713a8d71d8e|graphs/CX333|refbias/t11|refbias/p1f|"
                 r"accessory/assets|refbias/panel|refbias/t2|refbias/assets|"
                 r"GCF_035581225|graphframe/results|insgt/assets|complex333|"
                 r"gwas1000_accessory_census")


class NoOldDefaults(unittest.TestCase):

    def test_python_defaults(self):
        for p in PRODUCTION_PY:
            s = text(p)
            for m in re.finditer(r'default=\s*\(?\s*"([^"]*)"', s):
                self.assertIsNone(OLD.search(m.group(1)), (p, m.group(1)))
            for m in re.finditer(r'^\s*[A-Z_]+\s*=\s*\(?\s*"([^"]*)"', s, re.M):
                self.assertIsNone(OLD.search(m.group(1)), (p, m.group(1)))

    def test_shell_code(self):
        for p in PRODUCTION_SH:
            self.assertIsNone(OLD.search(code(p)), p)

    def test_required_where_there_was_a_default(self):
        cases = [
            ("is6110/bin/is6110_write_vcf.py", [], ["--refs", "--h37rv", "--build-id"]),
            ("is6110/bin/is6110_p5_merge.py", [], ["--h37rv"]),
            ("accessory/bin/merge_catalogues.py", ["--out", "o", "--out-fasta", "f"],
             ["--loci", "--graph-vcf"]),
            ("accessory/bin/locus_presence.py", ["--sample", "s", "--cram", "c",
                                                 "--h37rv", "h"],
             ["--catalogue", "--fasta"]),
            ("bin/merge_cohort_vcf.py", ["--cohort-name", "c", "--out", "o"],
             ["--build-id"]),
            ("bin/sv_intervals.py", ["--out", "o"], ["--graph-vcf"]),
            ("bin/t8_select_reference.py", ["--vcf", "v", "--out", "o"],
             ["--panel-snps"]),
            ("graphframe/bin/graph_frame_offsets.py", ["--out", "o"],
             ["--panel", "--refs"]),
            ("graphframe/bin/frame_detect.py", ["--vcf", "v", "--accession", "a"],
             ["--refs"]),
            ("assoc/bin/add_outgroup.py", ["--alignment", "a", "--sites", "s",
                                           "--out", "o", "--panel-vcf", "p"],
             ["--outgroup"]),
            ("bin/retier_intervals.py", ["--iv-states", "i", "--refmap", "r",
                                         "--out", "o"], ["--intervals"]),
            ("refbias/bin/io_contract.py", ["coh"], ["--build"]),
        ]
        for p, args, need in cases:
            r = run([sys.executable, p, *args])
            self.assertEqual(r.returncode, 2, (p, r.stderr))
            for n in need:
                self.assertIn(n, r.stderr, p)

    def test_stage2_projection_needs_the_builds_graph(self):
        with tempfile.TemporaryDirectory() as d:
            write(f"{d}/graphs/CX333.s10k.k23.K15/x.smooth.final.og", "old\n")
            write(f"{d}/refmap.tsv", "sample\treference\nS1\tR1\n")
            for x in ("bin", "graphframe"):
                os.symlink(os.path.join(ROOT, x), f"{d}/{x}")
            r = run([sys.executable, os.path.join(ROOT, "is6110/bin/is6110_p5_stage2.py"),
                     "--mode", "project", "--ref", "R1", "--refmap", f"{d}/refmap.tsv",
                     "--workdir", f"{d}/w"], cwd=d)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("needs --graph and --paths", r.stderr)

    def test_merge_catalogues_reads_no_cx333_side_tables(self):
        with tempfile.TemporaryDirectory() as d:
            # CX333's insgt tables where the old defaults looked: not read
            write(f"{d}/insgt/assets/insertions.tsv", "contig\th37rv_pos\nX\t5000\n")
            write(f"{d}/loci.tsv", "locus_id\tpos\tklass\trep_len\tn_alleles\t"
                  "carriers_any\tcarrier_frac\th37rv_cov\tnovelty\n"
                  "ACC_5000\t5000\tpolymorphic\t400\t1\t2\t0.6\t0\tnovel\n")
            g = vcf_gz(f"{d}/g.vcf.gz", ["A1", "A2"])
            r = run([sys.executable, os.path.join(ROOT, "accessory/bin/merge_catalogues.py"),
                     "--loci", "loci.tsv", "--graph-vcf", g, "--out", "c.tsv",
                     "--out-fasta", "c.fa"], cwd=d)
            self.assertEqual(r.returncode, 0, r.stderr)
            c = rows(f"{d}/c.tsv")
        self.assertEqual((c[0]["cluster"], c[0]["n_clusters"]), ("", "0"))

    def test_merge_presence_needs_a_catalogue(self):
        r = run([sys.executable, "bin/merge_cohort_vcf.py", "--cohort-name", "c",
                 "--build-id", "b", "--accessory-presence", "x", "--out", "o"])
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("needs --accessory-catalogue", r.stderr)

    def test_runner_passes_the_build_to_the_io_check(self):
        self.assertRegex(code("bin/refbias_run.sh"),
                         r'io_contract\.py[^\n]*\\\n[^\n]*--build "\$BUILD"')

    def test_p3_graph_from_the_build(self):
        c = code("bin/p3_accessory.sh")
        self.assertIn('$1=="graph"', c)
        self.assertIn("build_info.tsv", c)


class FrameTableFromTheBuild(unittest.TestCase):
    """graph_frame.Frames read graphframe/results/ (CX333's frame table)
    whenever MTB_GRAPH_FRAMES was unset."""

    def load(self, **kw):
        old = dict(os.environ)
        for k in ("MTB_GRAPH_FRAMES", "MTB_BUILD_DIR"):
            os.environ.pop(k, None)
        os.environ.update(kw)
        try:
            spec = importlib.util.spec_from_file_location(
                "graph_frame_c2", os.path.join(ROOT, "graphframe/bin/graph_frame.py"))
            m = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(m)
            return m
        finally:
            os.environ.clear()
            os.environ.update(old)

    def test_no_table_without_a_build(self):
        with tempfile.TemporaryDirectory() as d:
            write(f"{d}/graphframe/results/graph_frame_offsets.tsv",
                  "accession\tstrand\toffset\tpanel_len\nOLD\t+\t0\t1\n")
            cwd = os.getcwd()
            os.chdir(d)
            try:
                m = self.load()
                with self.assertRaises(SystemExit):
                    m.Frames()
            finally:
                os.chdir(cwd)

    def test_the_builds_table(self):
        with tempfile.TemporaryDirectory() as d:
            write(f"{d}/b/assets/graph_frame_offsets.tsv",
                  "accession\tstrand\toffset\tpanel_len\nNEW\t+\t0\t1\n")
            m = self.load(MTB_BUILD_DIR=f"{d}/b")
            old = dict(os.environ)
            os.environ.pop("MTB_GRAPH_FRAMES", None)
            try:
                self.assertTrue(m.Frames().known("NEW"))
            finally:
                os.environ.clear()
                os.environ.update(old)


class NodeOffsetComment(unittest.TestCase):
    """node_path_membership.py said P4's node offsets run along the path, so
    the base is start + offset on either strand; they are forward offsets."""

    def test_comment_matches_the_forward_offset(self):
        s = text("accessory/bin/node_path_membership.py")
        self.assertNotIn("start + offset on\n        # either strand", s)
        self.assertNotIn("already run along the path", s)
        self.assertIn("forward_offset", s)


if __name__ == "__main__":
    unittest.main(warnings="ignore")
