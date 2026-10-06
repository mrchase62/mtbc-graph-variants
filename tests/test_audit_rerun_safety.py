#!/usr/bin/env python3
"""Regression tests for the 2026-10-05 audit, group "rerun_safety" (section B).

A rerun on a new graph, into new output folders, must not be able to pick up a
previous build's outputs or the CX333 graph's assets:

  P0P2-1  P1 and P2 "already done" ignored the build and the reference
  TP-4    P0 step markers and the association tail's guards tested existence
  P3IS-8  accessory presence tables kept whatever catalogue made them
  P0P2-3  OG alone left GRAPH_DIR (and the panel SNPs) at CX333
  P0P2-2  build assets were links to files that changed underneath
  P0P2-4  the P1 tie-break read crossmaps outside the build, silently
  P4P5-7  node tables and panel AF were CX333 files, not build assets
  PGB-14  stamp_build_id.sh stamped any VCF with "the only" build

Same pattern as tests/run_tests.py: small synthetic inputs, standard library
plus the pipeline's own scripts, seconds to run.

    $MTB_PY tests/test_audit_rerun_safety.py
"""
import csv
import gzip
import hashlib
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

H37 = "GCF_000195955"


def write(path, text):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write(textwrap.dedent(text).lstrip("\n"))
    return path


def py():
    return os.environ.get("MTB_PY", sys.executable)


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def bcftools():
    b = os.environ.get("MTB_BCFTOOLS") or shutil.which("bcftools") or ""
    if not b:
        for line in open(os.path.join(ROOT, "config", "project_env.sh")):
            if "MTB_BCFTOOLS:=" in line:
                b = line.split(":=", 1)[1].split("}")[0]
    return b if b and os.access(b, os.X_OK) else ""


def vcf_text(samples, build=None, extra=""):
    h = "##fileformat=VCFv4.2\n##contig=<ID=c,length=1000>\n"
    h += '##FORMAT=<ID=GT,Number=1,Type=String,Description="">\n'
    if build:
        h += f"##MTB_graph_build={build}\n"
    h += extra
    h += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    h += "".join("\t" + s for s in samples) + "\n"
    h += "c\t10\t.\tA\tG\t60\t.\t.\tGT" + "".join("\t1" for _ in samples) + "\n"
    return h


def write_vcf_gz(path, samples, build=None, extra=""):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write(vcf_text(samples, build, extra))
    return path


def make_build(root, bid, graph_sha="", complete=True):
    b = os.path.join(root, bid)
    for d in ("assets", "logs", "refs", "annotation"):
        os.makedirs(os.path.join(b, d), exist_ok=True)
    write(os.path.join(b, "build_info.tsv"),
          f"build_id\t{bid}\ngraph\t/nonexistent/{bid}.og\n"
          f"graph_sha256\t{graph_sha or bid * 5}\n")
    if complete:
        for s in ("refs", "assets", "manifest"):
            write(os.path.join(b, "logs", f"{s}.done"), "2026-10-05\n")
    return b


def env(**kw):
    e = dict(os.environ)
    e.update(MTB_SITE_FILE="/dev/null", MTB_PY=py())
    for k in ("MTB_BUILD_DIR", "OG", "GRAPH_DIR", "BUILD_ID", "SLURM_ARRAY_TASK_ID",
              "PANEL_SNPS", "GRAPH_VCF", "MTB_GRAPH_FRAMES"):
        e.pop(k, None)
    if bcftools():
        e["MTB_BCFTOOLS"] = bcftools()
    e.update({k: str(v) for k, v in kw.items()})
    return e


def run(cmd, cwd=None, **kw):
    return subprocess.run(cmd, cwd=cwd, env=env(**kw), capture_output=True,
                          text=True)


# --------------------------------------------------------------------------
class P1Guard(unittest.TestCase):
    """P0P2-1: P1 skips only outputs of the current build."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        d = self.d
        self.build = make_build(os.path.join(d, "build"), "newbuild")
        write(os.path.join(self.build, "refs", f"{H37}.fasta"), ">c\nACGT\n")
        write(os.path.join(self.build, "assets", "panel_snps.vcf.gz"), "x")
        os.makedirs(os.path.join(d, "crams"))
        write(os.path.join(d, "crams", "s.cram"), "x")
        write(os.path.join(d, "cramref.fa"), ">c\nA\n")
        write(os.path.join(d, "cohort.tsv"), "sample\tlineage\nS1\tlineage4\n")
        write(os.path.join(d, "crams.tsv"), "S1\ts.cram\n")
        self.out = os.path.join(d, "p1")
        self.work = os.path.join(d, "work")
        os.makedirs(self.out); os.makedirs(self.work)
        write(os.path.join(self.out, "S1.candidates.tsv"),
              "rank\treference\tsnp_distance\n1\tR1\t5\n2\tR2\t9\n")

    def tearDown(self):
        shutil.rmtree(self.d)

    def p1(self):
        return run(["bash", os.path.join(ROOT, "bin/p1_select_reference.sh"), "S1"],
                   cwd=self.d, MTB_BUILD_DIR=self.build, COHORT=f"{self.d}/cohort.tsv",
                   CRAMMAP=f"{self.d}/crams.tsv", MTB_CRAM_ROOT=f"{self.d}/crams",
                   MTB_CRAM_REF=f"{self.d}/cramref.fa", OUTDIR=self.out,
                   WORK=self.work)

    def test_old_build_outputs_refused(self):
        write_vcf_gz(os.path.join(self.work, "S1.h37rv.vcf.gz"), ["S1"], "oldbuild")
        r = self.p1()
        self.assertNotEqual(r.returncode, 0, r.stdout)
        self.assertIn("another build", r.stderr)

    def test_old_marker_refused(self):
        write_vcf_gz(os.path.join(self.work, "S1.h37rv.vcf.gz"), ["S1"], "newbuild")
        write(os.path.join(self.out, "S1.p1.done"), "build_id\toldbuild\n")
        r = self.p1()
        self.assertNotEqual(r.returncode, 0, r.stdout)

    def test_current_build_legacy_outputs_kept_and_recorded(self):
        write_vcf_gz(os.path.join(self.work, "S1.h37rv.vcf.gz"), ["S1"], "newbuild")
        r = self.p1()
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("already done", r.stdout)
        m = open(os.path.join(self.out, "S1.p1.done")).read()
        self.assertIn("build_id\tnewbuild", m)
        self.assertIn("reference\tR1", m)


class P2Guard(unittest.TestCase):
    """P0P2-1: P2 skips only outputs of this build AND this reference."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        d = self.d
        self.build = make_build(os.path.join(d, "build"), "newbuild")
        for ext in ("", ".bwt", ".fai"):
            write(os.path.join(self.build, "refs", f"R1.fasta{ext}"), "x\n")
        write(os.path.join(self.build, "refs", "R1.dict"), "x\n")
        os.makedirs(os.path.join(d, "crams"))
        write(os.path.join(d, "crams", "s.cram"), "x")
        write(os.path.join(d, "cramref.fa"), ">c\nA\n")
        write(os.path.join(d, "crams.tsv"), "S1\ts.cram\n")
        self.refmap = write(os.path.join(d, "refmap.tsv"),
                            "sample\ta\tb\tc\treference\nS1\t\t\t\tR1\n")
        self.out = os.path.join(d, "p2")
        os.makedirs(self.out)

    def tearDown(self):
        shutil.rmtree(self.d)

    def p2(self):
        return run(["bash", os.path.join(ROOT, "bin/p2_call.sh"), "S1"], cwd=self.d,
                   MTB_BUILD_DIR=self.build, REFMAP=self.refmap,
                   CRAMMAP=f"{self.d}/crams.tsv", MTB_CRAM_ROOT=f"{self.d}/crams",
                   MTB_CRAM_REF=f"{self.d}/cramref.fa", OUTDIR=self.out,
                   WORK=f"{self.d}/work", DELLY_ENV=self.d,
                   DYSGU=f"{self.d}/no_dysgu")

    def outputs(self, build, ref):
        extra = ('##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="'
                 f'HaplotypeCaller --reference refbias/build/{build}/refs/{ref}.fasta '
                 '--output x.vcf.gz",Version="4">\n')
        write_vcf_gz(os.path.join(self.out, "S1.vcf.gz"), ["S1"], build, extra)
        write(os.path.join(self.out, "S1.delly.vcf"), vcf_text(["S1"], build))

    def test_date_only_marker_from_old_build_refused(self):
        self.outputs("oldbuild", "R1")
        write(os.path.join(self.out, "S1.p2.done"), "2026-09-30T10:00:00\n")
        r = self.p2()
        self.assertNotEqual(r.returncode, 0, r.stdout)
        self.assertIn("another build or reference", r.stderr)

    def test_other_reference_refused(self):
        self.outputs("newbuild", "R2")
        r = self.p2()
        self.assertNotEqual(r.returncode, 0, r.stdout)
        self.assertIn("called against", r.stderr)

    def test_marker_names_other_reference_refused(self):
        self.outputs("newbuild", "R1")
        write(os.path.join(self.out, "S1.p2.done"),
              "build_id\tnewbuild\nreference\tR2\n")
        self.assertNotEqual(self.p2().returncode, 0)

    def test_current_outputs_kept_and_marker_records_both(self):
        self.outputs("newbuild", "R1")
        r = self.p2()
        self.assertEqual(r.returncode, 0, r.stderr)
        m = open(os.path.join(self.out, "S1.p2.done")).read()
        self.assertIn("build_id\tnewbuild", m)
        self.assertIn("reference\tR1", m)

    def test_refmap_from_other_build_refused(self):
        write(self.refmap + ".build", "build_id\toldbuild\n")
        r = self.p2()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("selected against build oldbuild", r.stderr)


class P1Summary(unittest.TestCase):
    """P0P2-4 and P0P2-1: interval counts from the build, stale P1 refused."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        write(os.path.join(self.d, "cohort.tsv"), "sample\tlineage\nS1\tlineage4\n")
        write(os.path.join(self.d, "S1.candidates.tsv"),
              "rank\treference\tsnp_distance\n1\tR1\t5\n2\tR2\t7\n")
        self.out = os.path.join(self.d, "refmap.tsv")

    def tearDown(self):
        shutil.rmtree(self.d)

    def summ(self, *extra):
        return subprocess.run([py(), "bin/p1_summary.py", "--cohort",
                               f"{self.d}/cohort.tsv", "--dir", self.d,
                               "--lineages", f"{self.d}/none.csv", "--out",
                               self.out, *extra], capture_output=True, text=True)

    def test_build_intervals_drive_tie_break(self):
        iv = write(os.path.join(self.d, "iv.tsv"),
                   "reference\tintervals\nR1\t3\nR2\t9\n")
        write(os.path.join(self.d, "S1.p1.done"), "build_id\tB\n")
        r = self.summ("--intervals", iv, "--build-id", "B")
        self.assertEqual(r.returncode, 0, r.stderr)
        row = next(csv.DictReader(open(self.out), delimiter="\t"))
        self.assertEqual(row["reference"], "R2")
        self.assertEqual(open(self.out + ".build").read(), "build_id\tB\n")

    def test_missing_interval_count_is_fatal(self):
        iv = write(os.path.join(self.d, "iv.tsv"), "reference\tintervals\nR1\t3\n")
        r = self.summ("--intervals", iv)
        self.assertEqual(r.returncode, 2)
        self.assertIn("no IS6110 interval count", r.stderr)
        self.assertFalse(os.path.exists(self.out))

    def test_no_interval_source_is_fatal(self):
        r = self.summ()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("needs --intervals", r.stderr)

    def test_candidates_from_other_build_refused(self):
        iv = write(os.path.join(self.d, "iv.tsv"),
                   "reference\tintervals\nR1\t3\nR2\t9\n")
        write(os.path.join(self.d, "S1.p1.done"), "build_id\told\n")
        r = self.summ("--intervals", iv, "--build-id", "B")
        self.assertEqual(r.returncode, 2)
        self.assertIn("another build", r.stderr)


# --------------------------------------------------------------------------
class P0(unittest.TestCase):
    """TP-4, P0P2-2, P0P2-3: P0 markers, graph/dir consistency, copied assets."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        d = self.d
        os.symlink(os.path.join(ROOT, "bin"), os.path.join(d, "bin"))
        self.gdir = os.path.join(d, "newgraph")
        self.og = write(os.path.join(self.gdir, "g.smooth.final.og"), "graph bytes\n")
        self.gsha = sha(self.og)
        self.bid = self.gsha[:12]
        self.b = os.path.join(d, "build", self.bid)
        os.makedirs(os.path.join(self.b, "assets"))
        os.makedirs(os.path.join(self.b, "logs"))
        write(os.path.join(self.b, "build_info.tsv"),
              f"build_id\t{self.bid}\ngraph_sha256\t{self.gsha}\n")
        write(os.path.join(self.b, "assets", "accessions.txt"), f"A1\nA2\n{H37}\n")
        write(os.path.join(self.b, "logs", "accessions.done"), "2026-10-01\n")

    def tearDown(self):
        shutil.rmtree(self.d)

    def p0(self, *args, **kw):
        kw.setdefault("BUILD_ROOT", os.path.join(self.d, "build"))
        kw.setdefault("MTB_GATK_SIF", "/nonexistent.sif")
        return run(["bash", os.path.join(ROOT, "bin/p0_prepare.sh"), *args],
                   cwd=self.d, **kw)

    def acc(self):
        # an accessory panel made from this graph's genomes (cleanup2: step
        # assets requires one, with its genome record)
        a = os.path.join(self.d, "acc")
        write(os.path.join(a, "panel_manifest.tsv"), "locus_id\tpos\n")
        write(os.path.join(a, "panel_manifest.genomes.txt"), "A1\nA2\n")
        return a

    def test_og_and_other_graph_dir_refused(self):
        os.makedirs(os.path.join(self.d, "oldgraph"))
        r = self.p0("--list", OG=self.og, GRAPH_DIR=os.path.join(self.d, "oldgraph"))
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("is not in GRAPH_DIR", r.stderr)

    def test_ancestral_date_only_marker_refused(self):
        al = write(os.path.join(self.d, "aln.fasta"), f">A1\nA\n>A2\nA\n>{H37}\nA\n")
        t = write(os.path.join(self.d, "t.nwk"), "(A1,A2,GCF_000195955);\n")
        st = write(os.path.join(self.d, "sites.tsv"), "column\tchrom\tpos\tref\talt\n")
        write(os.path.join(self.b, "logs", "ancestral.done"), "2026-10-01T12:40:52\n")
        r = self.p0("--step", "ancestral", OG=self.og, ANC_TREE=t, ANC_ALN=al,
                    ANC_SITES=st)
        self.assertNotEqual(r.returncode, 0, r.stdout)
        self.assertIn("step 'ancestral' was done with", r.stderr)

    def test_ancestral_alignment_from_other_panel_refused(self):
        al = write(os.path.join(self.d, "aln.fasta"), ">X9\nA\n>A2\nA\n")
        t = write(os.path.join(self.d, "t.nwk"), "(X9,A2);\n")
        st = write(os.path.join(self.d, "sites.tsv"), "column\tchrom\tpos\tref\talt\n")
        r = self.p0("--step", "ancestral", OG=self.og, ANC_TREE=t, ANC_ALN=al,
                    ANC_SITES=st)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("taxa of", r.stderr)

    @unittest.skipUnless(bcftools(), "bcftools not available")
    def test_assets_from_ogs_own_dir_copied_and_checked(self):
        snps = write_vcf_gz(os.path.join(self.gdir, "snps.vcf.gz"), ["A1", "A2"])
        write_vcf_gz(os.path.join(self.gdir, "all_variants.collapsed.vcf.gz"),
                     ["A1", "A2"])
        mask = write(os.path.join(self.d, "mask.bed"), "c\t0\t5\n")
        r = self.p0("--step", "assets", OG=self.og, MASK=mask,
                    ACCESSORY_DIR=self.acc())
        self.assertEqual(r.returncode, 0, r.stderr)
        a = os.path.join(self.b, "assets")
        for f in ("panel_snps.vcf.gz", "graph_collapsed.vcf.gz", "repeat_mask.bed"):
            self.assertFalse(os.path.islink(os.path.join(a, f)), f)
        self.assertEqual(sha(os.path.join(a, "panel_snps.vcf.gz")), sha(snps))
        # the source changes; the build's copy does not
        write(mask, "c\t0\t500\n")
        self.assertEqual(open(os.path.join(a, "repeat_mask.bed")).read(), "c\t0\t5\n")
        # and the step is now refused, not silently "already done"
        r = self.p0("--step", "assets", OG=self.og, MASK=mask,
                    ACCESSORY_DIR=self.acc())
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("mask_sha256", r.stderr)

    @unittest.skipUnless(bcftools(), "bcftools not available")
    def test_panel_snps_of_another_graph_refused(self):
        write_vcf_gz(os.path.join(self.gdir, "snps.vcf.gz"), ["A1", "Z7"])
        write_vcf_gz(os.path.join(self.gdir, "all_variants.collapsed.vcf.gz"),
                     ["A1", "A2"])
        mask = write(os.path.join(self.d, "mask.bed"), "c\t0\t5\n")
        r = self.p0("--step", "assets", OG=self.og, MASK=mask,
                    ACCESSORY_DIR=self.acc())
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not in the build", r.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.b, "logs", "assets.done")))


class BuildVerify(unittest.TestCase):
    """P0P2-2: a build is checked against its manifest at use."""

    def test_changed_and_linked_assets_refused(self):
        with tempfile.TemporaryDirectory() as d:
            b = os.path.join(d, "B")
            m = write(os.path.join(b, "assets", "repeat_mask.bed"), "c\t0\t5\n")
            write(os.path.join(b, "manifest.tsv"),
                  "asset\tpath\tbytes\tsha256\n"
                  f"repeat_mask.bed\trefbias/build/B/assets/repeat_mask.bed\t8\t{sha(m)}\n")
            ok = subprocess.run([py(), "bin/p0_check.py", "verify", "--build", b],
                                capture_output=True, text=True)
            self.assertEqual(ok.returncode, 0, ok.stderr)
            write(m, "c\t0\t500\n")
            bad = subprocess.run([py(), "bin/p0_check.py", "verify", "--build", b],
                                 capture_output=True, text=True)
            self.assertNotEqual(bad.returncode, 0)
            self.assertIn("differs from the manifest", bad.stderr)
            src = write(os.path.join(d, "outside.bed"), "c\t0\t5\n")
            os.remove(m); os.symlink(src, m)
            bad = subprocess.run([py(), "bin/p0_check.py", "verify", "--build", b],
                                 capture_output=True, text=True)
            self.assertIn("links outside the build", bad.stderr)

    def test_crossmap_from_other_reference_refused(self):
        with tempfile.TemporaryDirectory() as d:
            write(os.path.join(d, "acc.txt"), f"R1\n{H37}\n")
            write(os.path.join(d, "refs", "R1.fasta.fai"), "NC_1\t1000\t5\t60\t61\n")
            write(os.path.join(d, "is", "R1.crossmap.tsv"),
                  "orig_start\torig_end\tdeleted_len\tclean_junction\tcum_deleted\n"
                  "10\t109\t100\t9\t100\n")
            write(os.path.join(d, "is", "R1.isclean.fasta.fai"), "NC_1_isclean\t900\t5\t60\t61\n")
            args = [py(), "bin/p0_check.py", "is6110-intervals", "--isclean-dir",
                    f"{d}/is", "--refs", f"{d}/refs", "--accessions", f"{d}/acc.txt",
                    "--skip", H37, "--out", f"{d}/iv.tsv"]
            r = subprocess.run(args, capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            row = next(csv.DictReader(open(f"{d}/iv.tsv"), delimiter="\t"))
            self.assertEqual((row["reference"], row["intervals"]), ("R1", "1"))
            write(os.path.join(d, "is", "R1.isclean.fasta.fai"), "NC_1\t950\t5\t60\t61\n")
            r = subprocess.run(args, capture_output=True, text=True)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("not made from this build's references", r.stderr)


# --------------------------------------------------------------------------
class NodeTables(unittest.TestCase):
    """P4P5-7: node tables are made from the build's graph and passed by P5."""

    GFA = ("H\tVN:Z:1.0\nS\t1\tAAAA\nS\t2\tCC\nS\t3\tGGG\nS\t4\tT\n"
           "P\tH#1#c\t1+,3+,4+\t*\n"
           "P\tA#1#c\t1+,2+,3+,2+\t*\n"
           "P\tB#1#c\t4-,2-,1-\t*\n")

    def test_membership_and_positions(self):
        with tempfile.TemporaryDirectory() as d:
            g = write(os.path.join(d, "g.gfa"), self.GFA)
            r = subprocess.run([py(), "accessory/bin/node_path_membership.py", "--gfa",
                                g, "--exclude-path", "H#1#c", "--out", f"{d}/np.tsv",
                                "--out-positions", f"{d}/pos.tsv"],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            np_ = {x["node"]: x["paths"] for x in
                   csv.DictReader(open(f"{d}/np.tsv"), delimiter="\t")}
            self.assertEqual(np_, {"2": "A,B"})      # nodes 1,3,4 are on H
            pos = {(x["node"], x["accession"]): x for x in
                   csv.DictReader(open(f"{d}/pos.tsv"), delimiter="\t")}
            self.assertEqual(pos[("2", "A")]["start"], "5")
            self.assertEqual(pos[("2", "A")]["n_occurrences"], "2")
            self.assertEqual(pos[("2", "B")]["start"], "2")
            self.assertEqual(pos[("2", "B")]["strand"], "-")

    def test_no_cx333_default_graph(self):
        r = subprocess.run([py(), "accessory/bin/node_path_membership.py",
                            "--out", "/dev/null"], capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)

    def test_p5_passes_build_assets(self):
        s = open("bin/p5_merge.sh").read()
        self.assertIn('--node-paths "$NODE_PATHS"', s)
        self.assertIn('NODE_PATHS="${BUILD}/assets/node_paths.tsv"', s)
        self.assertIn('--graph-vcf "$GRAPH_VCF"', s)
        self.assertFalse("ls graphs/CX333" in s, "p5_merge.sh globs the CX333 graph")
        self.assertNotIn("graphs/CX333", open("bin/p5_matrix.py").read().split(
            "def main")[1].split("a = ap.parse_args()")[0])


# --------------------------------------------------------------------------
@unittest.skipUnless(bcftools(), "bcftools not available")
class Stamp(unittest.TestCase):
    """PGB-14: --graph stamps the graph's own build, or nothing."""

    def test_graph_build(self):
        with tempfile.TemporaryDirectory() as d:
            gA = write(os.path.join(d, "A.og"), "graph A\n")
            gB = write(os.path.join(d, "B.og"), "graph B\n")
            idA = sha(gA)[:12]
            b = make_build(os.path.join(d, "build"), idA, sha(gA))
            v = write(os.path.join(d, "x.vcf"), vcf_text(["S"]))
            r = run(["bash", "bin/stamp_build_id.sh", "--graph", gB, v],
                    BUILD_ROOT=os.path.join(d, "build"))
            self.assertNotEqual(r.returncode, 0)
            self.assertNotIn("MTB_graph_build", open(v).read())
            r = run(["bash", "bin/stamp_build_id.sh", "--graph", gA, v],
                    BUILD_ROOT=os.path.join(d, "build"))
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn(f"##MTB_graph_build={idA}", open(v).read())
            # MTB_BUILD_DIR naming a directory without a record: an error
            v2 = write(os.path.join(d, "y.vcf"), vcf_text(["S"]))
            r = run(["bash", "bin/stamp_build_id.sh", v2],
                    MTB_BUILD_DIR=os.path.join(d, "nobuild"))
            self.assertNotEqual(r.returncode, 0)
            _ = b


# --------------------------------------------------------------------------
class PresenceGuard(unittest.TestCase):
    """P3IS-8: presence tables are kept only for this build's catalogue."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        d = self.d
        self.b = make_build(os.path.join(d, "build"), "newbuild")
        for ext, txt in (("tsv", "locus_id\nACC_1\n"), ("fasta", ">ACC_1\nAC\n")):
            write(os.path.join(d, "accessory", "assets", f"accessory_catalogue.{ext}"), txt)
            write(os.path.join(self.b, "assets", f"accessory_catalogue.{ext}"), txt)
        self.refmap = write(os.path.join(d, "refmap.tsv"), "sample\treference\nS1\tR1\n")
        self.out = os.path.join(d, "accessory", "coh", "S1.presence.tsv")

    def tearDown(self):
        shutil.rmtree(self.d)

    def arr(self, **kw):
        return run(["bash", os.path.join(ROOT, "accessory/bin/locus_presence_array.sh")],
                   SLURM_SUBMIT_DIR=self.d, SLURM_ARRAY_TASK_ID=1, COHORT_TAG="coh",
                   REFMAP=self.refmap, MTB_BUILD_DIR=self.b, **kw)

    def test_table_older_than_catalogue_refused(self):
        write(self.out, "locus_id\tcall\n")
        os.utime(self.out, (1, 1))
        r = self.arr()
        self.assertNotEqual(r.returncode, 0, r.stdout)
        self.assertIn("was made against build", r.stderr)

    def test_table_of_this_build_kept(self):
        write(self.out, "locus_id\tcall\n")
        cs = hashlib.sha256(b"locus_id\nACC_1\n>ACC_1\nAC\n").hexdigest()[:16]
        write(self.out + ".build", f"build_id\tnewbuild\ncatalogue_sha\t{cs}\n")
        r = self.arr()
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("already done", r.stdout)

    def test_catalogue_not_the_builds_refused(self):
        write(os.path.join(self.d, "accessory", "assets", "accessory_catalogue.tsv"),
              "locus_id\nACC_2\n")
        # the step reads the build's copy now (audit cleanup); a catalogue
        # named explicitly must still be the build's
        r = self.arr(ACC_CATALOGUE=os.path.join(self.d, "accessory", "assets",
                                                "accessory_catalogue"))
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("differs from build", r.stderr)


# --------------------------------------------------------------------------
class AssocTailGuard(unittest.TestCase):
    """TP-4: the association tail refuses products not made from this VCF."""

    def test_unrecorded_alignment_refused(self):
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "assoc"))
            os.symlink(os.path.join(ROOT, "assoc", "bin"), os.path.join(d, "assoc", "bin"))
            os.symlink(os.path.join(ROOT, "bin"), os.path.join(d, "bin"))
            b = make_build(os.path.join(d, "refbias", "build"), "newbuild")
            write(os.path.join(b, "assets", "panel_polarity.tsv"), "chrom\n")
            write(os.path.join(b, "assets", "graph_collapsed.vcf.gz"), "x")
            write_vcf_gz(os.path.join(d, "refbias", "coh", "p5", "merged.vcf.gz"),
                         ["S1"], "newbuild")
            write(os.path.join(d, "data", "trees", "coh.snps.fasta"), ">S1\nA\n")
            r = run(["bash", os.path.join(ROOT, "assoc/bin/cohort_assoc_tail.sh"), "coh"],
                    cwd=d, BUILD_ROOT="refbias/build")
            self.assertNotEqual(r.returncode, 0, r.stdout)
            self.assertIn("was not made from", r.stderr)


# --------------------------------------------------------------------------
class RunnerGuard(unittest.TestCase):
    """P0P2-1 / P3IS-8 at the runner: one output root holds one build, the
    build is verified at use, and the frame table is the build's own."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        d = self.d
        self.b = make_build(os.path.join(d, "build"), "newbuild")
        write(os.path.join(self.b, "manifest.tsv"), "asset\tpath\tbytes\tsha256\n")
        write(os.path.join(self.b, "assets", "graph_frame_offsets.tsv"), "accession\n")
        og = write(os.path.join(d, "g.og"), "g\n")
        write(os.path.join(self.b, "build_info.tsv"),
              f"build_id\tnewbuild\ngraph\t{og}\n")
        write(os.path.join(d, "cohort.tsv"), "sample\tlineage\nS1\tl4\n")
        write(os.path.join(d, "crams.tsv"), "S1\tx.cram\n")
        self.outroot = os.path.join(d, "out")
        self.name = f"rsafety{os.getpid()}"
        write(os.path.join(d, "cohorts.tsv"),
              "cohort\ttable\tcrams\toutroot\tpasses\tnote\tworkprefix\n"
              f"{self.name}\t{d}/cohort.tsv\t{d}/crams.tsv\t{self.outroot}\tp1,p2\t-\t{d}/w_\n")

    def tearDown(self):
        shutil.rmtree(self.d)

    def runner(self):
        return run(["bash", "bin/refbias_run.sh", self.name, "--dry-run"],
                   REGISTRY=f"{self.d}/cohorts.tsv", MTB_BUILD_DIR=self.b,
                   IO_CHECK=0)

    def test_clean_outroot_runs(self):
        r = self.runner()
        self.assertEqual(r.returncode, 0, r.stderr)

    def test_p3acc_waits_for_p1(self):
        """Review 2, R2-INT-5: p3acc reads P1's refmap.tsv (its sample list and
        each sample's matched reference), so it is submitted after P1's
        summary, not with no dependency."""
        write(os.path.join(self.d, "cohorts.tsv"),
              "cohort\ttable\tcrams\toutroot\tpasses\tnote\tworkprefix\n"
              f"{self.name}\t{self.d}/cohort.tsv\t{self.d}/crams.tsv\t"
              f"{self.outroot}\tp1,p3\t-\t{self.d}/w_\n")
        r = self.runner()
        self.assertEqual(r.returncode, 0, r.stderr)
        line = [l for l in r.stderr.splitlines() if "locus_presence_array.sh" in l]
        self.assertTrue(line, r.stderr)
        self.assertIn("--dependency=afterok:DRYRUN_p1sum", line[0])

    def test_outroot_of_other_build_refused(self):
        write(os.path.join(self.outroot, ".mtb_build"), "build_id\toldbuild\n")
        r = self.runner()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("holds outputs of build oldbuild", r.stderr)

    def test_unrecorded_outroot_with_old_stamps_refused(self):
        write_vcf_gz(os.path.join(self.outroot, "p2", "S1.vcf.gz"), ["S1"], "oldbuild")
        r = self.runner()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("stamped 'oldbuild'", r.stderr)

    def test_no_frame_fallback_outside_the_build(self):
        os.remove(os.path.join(self.b, "assets", "graph_frame_offsets.tsv"))
        r = self.runner()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("no frame table", r.stderr)


if __name__ == "__main__":
    unittest.main(warnings="ignore")
