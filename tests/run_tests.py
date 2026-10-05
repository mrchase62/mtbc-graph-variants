#!/usr/bin/env python3
"""Regression tests for defects found in the 2026-09-29 review (CODE_REVIEW.md).

Every test here pins one bug that produced plausible-looking wrong output
rather than an error. They run on small synthetic inputs, need nothing beyond
the standard library and the pipeline's own scripts, and take a few seconds:

    $MTB_PY tests/run_tests.py            # from the repository root
    $MTB_PY tests/run_tests.py -v

Run them before any cohort run that follows a code change.
"""
import csv
import shutil
import importlib.util
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
import warnings

# the scripts under test open files without context managers; that is not
# what these tests are about
warnings.simplefilter("ignore", ResourceWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)                      # several scripts load siblings by relative path
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


class RepeatMaskLookup(unittest.TestCase):
    """CODE_REVIEW 5.1: nested intervals missed, and a 0/1-based slip."""

    def setUp(self):
        self.p4 = load("p4_place", "bin/p4_place.py")

    def test_nested_interval_outer_tail_is_masked(self):
        # PE_PGRS4 336359-339273 contains 336559-339142 (BED, 0-based half-open)
        iv = self.p4.merge_intervals([(336359, 339273), (336559, 339142)])
        hit = self.p4.make_hit(iv)
        self.assertTrue(hit(339200))          # past the inner end: was "core"

    def test_bed_edges_against_1_based_positions(self):
        hit = self.p4.make_hit(self.p4.merge_intervals([(100, 200)]))
        # BED [100, 200) covers 1-based 101..200
        self.assertFalse(hit(100))
        self.assertTrue(hit(101))
        self.assertTrue(hit(200))
        self.assertFalse(hit(201))

    def test_class_prefix_not_substring(self):
        with tempfile.TemporaryDirectory() as d:
            bed = write(os.path.join(d, "m.bed"), """
                c\t10\t20\tpe_ppe|PPE1|named
                c\t30\t40\ttandem|pepA|tandem100
                """)
            pe, other = self.p4.load_intervals(bed)
        self.assertEqual(pe, [(10, 20)])
        self.assertEqual(other, [(30, 40)])   # "pepA" contains PE; not PE/PPE


class GeneLookup(unittest.TestCase):
    """CODE_REVIEW 3.8: a gene containing a shorter nested gene was missed."""

    def test_nested_gene(self):
        p6 = load("p6_annotate", "bin/p6_annotate.py")
        rows = [(100, 1000, "outer", "", ""), (200, 300, "inner", "", "")]
        hit = p6.make_lookup(rows)
        self.assertEqual(hit(250)[2], "inner")
        self.assertEqual(hit(500)[2], "outer")   # was None: scan stopped at inner
        self.assertIsNone(hit(1001))


class Is6110Stage2Keys(unittest.TestCase):
    """CODE_REVIEW 4.1: node keys lost their offset and never found a carrier."""

    def test_node_key_keeps_offset(self):
        st2 = load("is6110_p5_stage2", "is6110/bin/is6110_p5_stage2.py")
        self.assertEqual(st2.short("node:46502:0:N><INS>"), "node:46502:0:")
        self.assertEqual(st2.short("h37rv:932204:C><INS>"), "h37rv:932204:")
        # the carrier table builds f'node:{node}:' from node = "<id>:<off>"
        self.assertEqual(st2.short("node:46502:0:N><INS>"), f'node:{"46502:0"}:')


class GvcfEvidence(unittest.TestCase):
    """CODE_REVIEW 3.3: any gVCF line counted as REF evidence."""

    def test_variant_lines_are_not_reference_coverage(self):
        p5 = load("p5_states", "bin/p5_states.py")
        with tempfile.TemporaryDirectory() as d:
            g = write(os.path.join(d, "s.g.vcf"), """
                ##fileformat=VCFv4.2
                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts
                c\t1\t.\tA\t<NON_REF>\t.\t.\tEND=99\tGT:DP\t0:30
                c\t100\t.\tC\tA,<NON_REF>\t.\t.\t.\tGT:DP\t1:30
                c\t101\t.\tGTT\tG,<NON_REF>\t.\t.\t.\tGT:DP\t1:30
                c\t102\t.\tT\t<NON_REF>\t.\t.\tEND=200\tGT:DP\t0:30
                """)
            blocks, spans, alt_at = p5.load_gvcf(g)
        cov = p5.make_cov(blocks, 5)
        inside = p5.make_in_spans(spans)
        self.assertTrue(cov(50))
        self.assertFalse(cov(100))            # a SNP call is not a reference block
        self.assertTrue(inside(100))
        self.assertTrue(inside(103))          # deleted base, though a block covers it
        self.assertFalse(inside(104))
        self.assertEqual(alt_at[100][:2], ("C", "A"))

    def test_spans_with_long_early_span(self):
        p5 = load("p5_states", "bin/p5_states.py")
        inside = p5.make_in_spans([(1, 100), (5, 6), (50, 55), (200, 201)])
        self.assertTrue(inside(99))           # reached only through the first span
        self.assertFalse(inside(150))


class SvReverseStrand(unittest.TestCase):
    """CODE_REVIEW 5.2: reverse-strand SVs came out start > end, shifted."""

    def run_p4b(self, d, proj_lines, svtype="DEL", pos=1000, end=3000):
        vcf = write(os.path.join(d, "s.delly.vcf"), f"""
            ##fileformat=VCFv4.2
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts
            R\t{pos}\t.\tN\t<{svtype}>\t60\tPASS\tSVTYPE={svtype};END={end};SVLEN={end - pos}\tGT\t1
            """)
        posf = write(os.path.join(d, "bp.pos"),
                     "#source.path.pos\ttarget.path.pos\tdist.to.ref\tstrand.vs.ref\n"
                     + "".join(proj_lines))
        mask = write(os.path.join(d, "m.bed"), "c\t1\t2\ttandem|x|y\n")
        out = os.path.join(d, "out.tsv")
        subprocess.run([sys.executable, "bin/p4b_place_sv.py", "--sample", "s",
                        "--reference", "R", "--build-id", "b", "--delly", vcf,
                        "--positions", posf, "--mask", mask, "--no-inherited",
                        "--out", out], check=True, capture_output=True)
        with open(out) as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    def test_reverse_deletion(self):
        # R bases 1001..3000 deleted; R runs reverse to H37Rv, R x <-> H 5001-x
        # (0-based in the file: R 999 -> H 4001, R 2999 -> H 2001)
        with tempfile.TemporaryDirectory() as d:
            rows = self.run_p4b(d, [
                "R#1#c,999,+\tH#1#h,4001,+\t0\t-\t+\n",
                "R#1#c,2999,+\tH#1#h,2001,+\t0\t-\t+\n"])
        r = rows[0]
        self.assertEqual(r["frame"], "h37rv")
        # deleted in H: f(3000)=2002 .. f(1001)=4001; padding base 2001
        self.assertEqual((int(r["h37rv_pos"]), int(r["h37rv_end"])), (2001, 4001))
        self.assertEqual(r["key"], "sv:DEL:2001:4001")

    def test_forward_unchanged(self):
        with tempfile.TemporaryDirectory() as d:
            rows = self.run_p4b(d, [
                "R#1#c,999,+\tH#1#h,1099,+\t0\t+\t+\n",
                "R#1#c,2999,+\tH#1#h,3099,+\t0\t+\t+\n"])
        self.assertEqual(rows[0]["key"], "sv:DEL:1100:3100")

    def test_span_mismatch_is_breakend(self):
        with tempfile.TemporaryDirectory() as d:
            rows = self.run_p4b(d, [
                "R#1#c,999,+\tH#1#h,1099,+\t0\t+\t+\n",
                "R#1#c,2999,+\tH#1#h,90099,+\t0\t+\t+\n"])
        self.assertEqual(rows[0]["frame"], "bnd")

    def test_no_calls_still_writes_header(self):
        with tempfile.TemporaryDirectory() as d:
            empty = write(os.path.join(d, "e.vcf"), "##fileformat=VCFv4.2\n")
            mask = write(os.path.join(d, "m.bed"), "c\t1\t2\ttandem|x|y\n")
            posf = write(os.path.join(d, "bp.pos"), "")
            out = os.path.join(d, "out.tsv")
            subprocess.run([sys.executable, "bin/p4b_place_sv.py", "--sample", "s",
                            "--reference", "R", "--build-id", "b", "--delly", empty,
                            "--positions", posf, "--mask", mask, "--no-inherited",
                            "--out", out], check=True, capture_output=True)
            with open(out) as fh:
                self.assertTrue(fh.readline().startswith("sample\treference"))

    def test_missing_graph_vcf_is_fatal(self):
        with tempfile.TemporaryDirectory() as d:
            mask = write(os.path.join(d, "m.bed"), "c\t1\t2\ttandem|x|y\n")
            posf = write(os.path.join(d, "bp.pos"), "")
            r = subprocess.run([sys.executable, "bin/p4b_place_sv.py", "--sample", "s",
                                "--reference", "R", "--build-id", "b",
                                "--positions", posf, "--mask", mask,
                                "--graph-vcf", os.path.join(d, "nope.vcf.gz"),
                                "--out", os.path.join(d, "o.tsv")],
                               capture_output=True, text=True)
            self.assertNotEqual(r.returncode, 0)


class Importers(unittest.TestCase):
    """Found regenerating scale200: renaming a p5_states function broke the
    SV genotyper, which imports it. Every script that imports a sibling must
    still load."""

    def test_sibling_imports_resolve(self):
        for rel in ("bin/p5_sv_genotype.py", "bin/p5_matrix.py",
                    "bin/merge_cohort_vcf.py", "bin/p5_keys.py", "bin/p5_states.py"):
            r = subprocess.run([sys.executable, "-c",
                                f"import runpy,sys; sys.argv=['x','--help']; "
                                f"runpy.run_path('{rel}', run_name='__main__')"],
                               capture_output=True, text=True)
            self.assertNotIn("ImportError", r.stderr, rel)
            self.assertNotIn("ModuleNotFoundError", r.stderr, rel)


class ShardedMerge(unittest.TestCase):
    """Item 2 of the scaling plan: the merged VCF built in shards and
    assembled must be the single-process VCF, byte for byte once decompressed,
    with shard boundaries falling inside both the H37Rv range and the node
    contigs. Also checks the states array round-trips the sparse files."""

    def setUp(self):
        self.bgzip = os.environ.get("MTB_BGZIP", "")
        self.tabix = os.environ.get("MTB_TABIX", "")
        if not (self.bgzip and self.tabix and os.path.exists(self.bgzip)):
            self.skipTest("MTB_BGZIP/MTB_TABIX not set; source config/project_env.sh")

    def build(self, d):
        import hashlib
        samples = ["S1", "S2", "S3", "S4", "S5"]
        keys = []
        for pos in (100, 900_000, 1_800_000, 2_700_000, 3_600_000, 4_411_000):
            keys.append(dict(key=f"h37rv:{pos}:A>G", frame="h37rv", h37rv_pos=pos,
                             node="", node_offset="", canonical_ref="A",
                             canonical_alt="G", region="core", kind="SNP",
                             acc_locus=""))
        for nid in ("5", "12", "300", "4001"):
            keys.append(dict(key=f"node:{nid}:0:C>T", frame="node", h37rv_pos="",
                             node=nid, node_offset="0", canonical_ref="C",
                             canonical_alt="T", region="off_path_near",
                             kind="SNP", acc_locus=""))
        cols = list(keys[0])
        with open(os.path.join(d, "keys.tsv"), "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for k in keys:
                fh.write("\t".join(str(k[c]) for c in cols) + "\n")
        sha = hashlib.sha1("\n".join(k["key"] for k in keys).encode()).hexdigest()
        states = {"S1": {0: "ALT", 6: "ALT", 3: "ABSENT"}, "S2": {1: "ALT", 7: "ALT"},
                  "S3": {2: "ALT", 8: "ALT", 4: "NOCALL"}, "S4": {3: "ALT", 9: "ALT"},
                  "S5": {4: "ALT", 5: "ALT", 0: "NOCALL"}}
        for sm in samples:
            with open(os.path.join(d, f"{sm}.states.tsv"), "w") as fh:
                fh.write(f"#format\tsparse-v1\n#sample\t{sm}\n#default\tREF\n"
                         f"#n_keys\t{len(keys)}\n#keys_sha1\t{sha}\n"
                         "idx\tstate\tallele\n")
                for i, st in sorted(states[sm].items()):
                    fh.write(f"{i}\t{st}\t\n")
        with open(os.path.join(d, "refmap.tsv"), "w") as fh:
            fh.write("sample\treference\n" + "".join(f"{sm}\tR\n" for sm in samples))
        subprocess.run([sys.executable, "bin/p5_matrix.py", "--refmap",
                        os.path.join(d, "refmap.tsv"), "--keys", os.path.join(d, "keys.tsv"),
                        "--dir", d, "--out", os.devnull, "--no-dense",
                        "--graph-vcf", os.path.join(d, "none.vcf.gz")],
                       check=True, capture_output=True)
        return samples, keys

    def merge(self, d, *extra):
        r = subprocess.run([sys.executable, "bin/merge_cohort_vcf.py",
                            "--states-array", d, "--keys", os.path.join(d, "keys.tsv"),
                            "--cohort-name", "t", "--ancestral", "",
                            "--bgzip", self.bgzip, "--tabix", self.tabix, *extra],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)

    def body(self, path):
        import gzip
        with gzip.open(path, "rt") as fh:
            return [l for l in fh if not l.startswith("##fileDate")]

    def test_array_round_trip(self):
        p5io = load("p5_states_io", "bin/p5_states_io.py")
        with tempfile.TemporaryDirectory() as d:
            samples, keys = self.build(d)
            M, got = p5io.open_array(d, rows(os.path.join(d, "keys.tsv")))
            self.assertEqual(got, samples)
            ref = p5io.load_states(d, rows(os.path.join(d, "keys.tsv")), samples, quiet=True)
            self.assertTrue((M == ref).all())

    def test_sharded_equals_unsharded(self):
        with tempfile.TemporaryDirectory() as d:
            self.build(d)
            one = os.path.join(d, "one.vcf.gz")
            self.merge(d, "--out", one)
            for n in (2, 3, 5):
                pre = os.path.join(d, f"p{n}")
                for i in range(n):
                    self.merge(d, "--n-shards", str(n), "--shard", str(i),
                               "--part-prefix", pre, "--out", os.devnull)
                out = os.path.join(d, f"sh{n}.vcf.gz")
                self.merge(d, "--assemble", "--n-shards", str(n),
                           "--part-prefix", pre, "--out", out)
                self.assertEqual(self.body(one), self.body(out), f"{n} shards")
            recs = [l for l in self.body(one) if not l.startswith("#")]
            self.assertEqual(len(recs), 10)            # every key has an ALT carrier

    def test_uncatalogued_caller_deletion_is_kept(self):
        """With a catalogue, only the caller deletions it covers are dropped;
        one it does not cover stays, presence-only and flagged UNCATALOGUED."""
        import gzip
        with tempfile.TemporaryDirectory() as d:
            samples, _ = self.build(d)
            hdr = "key\tsvtype\th37rv_pos\tsvlen\tcomponent\tqual_band\t" + "\t".join(samples)
            alt = "\t".join(["ALT"] + ["NOCALL"] * (len(samples) - 1))
            write(f"{d}/sv.tsv", hdr + "\n"
                  f"sv:DEL:5000:400\tDEL\t5000\t400\tcalled\tlow\t{alt}\n"
                  f"sv:DEL:800000:30\tDEL\t800000\t30\tcalled\tlow\t{alt}\n")
            write(f"{d}/iv.tsv", "interval\tsource\tsvtype\tstart\tend\tsvlen\t"
                  "support_tier\tis6110_prox\tn_ref_carriers\tref_carriers\n"
                  "svi:DEL:5001:400\tcaller\tDEL\t5001\t5400\t400\tC_caller_only\t0\t0\t\n")
            write(f"{d}/ivst.tsv", "sample\tinterval\tsvtype\tstate\n" + "".join(
                f"{sm}\tsvi:DEL:5001:400\tDEL\t{'ALT' if sm == 'S1' else 'REF'}\n"
                for sm in samples))
            out = f"{d}/o.vcf.gz"
            self.merge(d, "--sv-matrix", f"{d}/sv.tsv", "--sv-intervals", f"{d}/iv.tsv",
                       "--sv-interval-states", f"{d}/ivst.tsv", "--out", out)
            with gzip.open(out, "rt") as fh:
                ids = {l.split("\t")[2]: l for l in fh if not l.startswith("#")}
            self.assertIn("svi:DEL:5001:400", ids)
            self.assertNotIn("sv:DEL:5000:400", ids)          # superseded
            self.assertIn("sv:DEL:800000:30", ids)            # not covered: kept
            self.assertIn("UNCATALOGUED", ids["sv:DEL:800000:30"])

    def test_assemble_refuses_missing_shard(self):
        with tempfile.TemporaryDirectory() as d:
            self.build(d)
            pre = os.path.join(d, "p")
            self.merge(d, "--n-shards", "3", "--shard", "0", "--part-prefix", pre,
                       "--out", os.devnull)
            r = subprocess.run([sys.executable, "bin/merge_cohort_vcf.py", "--assemble",
                                "--n-shards", "3", "--part-prefix", pre, "--cohort-name",
                                "t", "--out", os.path.join(d, "x.vcf.gz")],
                               capture_output=True, text=True)
            self.assertNotEqual(r.returncode, 0)


def rows(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class Is6110Stage2Split(unittest.TestCase):
    """Stage 2 run per reference and per sample: the array edges."""

    def base(self, d):
        write(os.path.join(d, "refmap.tsv"), "sample\treference\nA\tR1\nB\tR2\n")
        write(os.path.join(d, "keys.tsv"),
              "sample\treference\tr_pos\tframe\th37rv_pos\tnode\n")
        return [sys.executable, "is6110/bin/is6110_p5_stage2.py",
                "--refmap", os.path.join(d, "refmap.tsv"),
                "--cohort-keys", os.path.join(d, "keys.tsv"),
                "--stage1-dir", d, "--paths", os.path.join(d, "paths.txt"),
                "--graph", os.path.join(d, "g.og"),
                "--workdir", os.path.join(d, "w"), "--out", os.path.join(d, "o.tsv")]

    def test_project_index_past_last_reference_exits_cleanly(self):
        with tempfile.TemporaryDirectory() as d:
            write(os.path.join(d, "paths.txt"), "R1#1#c\nR2#1#c\n")
            r = subprocess.run(self.base(d) + ["--mode", "project", "--index", "3"],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("nothing to do", r.stdout)

    def test_merge_refuses_a_missing_sample(self):
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "w", "samples"))
            write(os.path.join(d, "w", "samples", "A.tsv"),
                  "sample\tkey\tstate\tallele\tevidence\nA\tk\tREF\tN\tx\n")
            r = subprocess.run(self.base(d) + ["--mode", "merge"],
                               capture_output=True, text=True)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("no stage-2 table for B", r.stderr)

    def test_projection_store_reuse_and_seed(self):
        """Projections are reused by (carrier, position): a rerun sends odgi
        only new positions, results are identical, and seeding from an
        earlier run's projection files reproduces the store."""
        with tempfile.TemporaryDirectory() as d:
            write(f"{d}/frames.tsv", "accession\tpanel_len\trefs_len\tstrand\toffset\tagree\n"
                  + "".join(f"{x}\t1000\t1000\t+\t0\tok\n" for x in ("R1", "R2")))
            write(f"{d}/paths.txt", "R1#1#c\nR2#1#c\n")
            log = f"{d}/odgi.log"
            fake = write(f"{d}/odgi", "#!/usr/bin/env python3\n"
                "import sys\n"
                "a = sys.argv; q = a[a.index('-F') + 1]; t = a[a.index('-r') + 1]\n"
                f"log = open({log!r}, 'a')\n"
                "for l in open(q):\n"
                "    p, o, s = l.strip().rsplit(',', 2); log.write(l)\n"
                "    print(f'{p},{o},{s}\\t{t},{int(o) + 7},+\\t0')\n")
            os.chmod(fake, 0o755)
            write(f"{d}/refmap.tsv", "sample\treference\nA\tR1\nB\tR2\n")
            def keys(pos):
                write(f"{d}/keys.tsv", "sample\treference\tr_pos\tframe\th37rv_pos\tnode\n"
                      + "".join(f"A\tR1\t{p}\th37rv\t{p}\t\n" for p in pos))
                os.makedirs(f"{d}/s1", exist_ok=True)
                write(f"{d}/s1/A.tsv", "sample\tkey\tstate\n"
                      + "".join(f"A\th37rv:{p}:A>T\tALT\n" for p in pos))
                write(f"{d}/s1/B.tsv", "sample\tkey\tstate\n"
                      + "".join(f"B\th37rv:{p}:A>T\tNOCALL\n" for p in pos))
            env = dict(os.environ, MTB_GRAPH_FRAMES=f"{d}/frames.tsv")
            cmd = [sys.executable, "is6110/bin/is6110_p5_stage2.py",
                   "--refmap", f"{d}/refmap.tsv", "--cohort-keys", f"{d}/keys.tsv",
                   "--stage1-dir", f"{d}/s1", "--paths", f"{d}/paths.txt",
                   "--graph", f"{d}/g.og", "--odgi", fake, "--workdir", f"{d}/w",
                   "--store", f"{d}/store", "--out", f"{d}/o.tsv"]
            def project():
                r = subprocess.run(cmd + ["--mode", "project", "--ref", "R2"],
                                   capture_output=True, text=True, env=env)
                self.assertEqual(r.returncode, 0, r.stderr)
                return open(f"{d}/w/proj/R2.tsv").read()
            sent = lambda: len(open(log).read().splitlines()) if os.path.exists(log) else 0
            keys([100, 200])
            first = project()
            self.assertEqual(sent(), 2)
            self.assertIn("h37rv:100:\t107\t0", first)
            keys([100, 200, 300])                      # one new key
            second = project()
            self.assertEqual(sent(), 3)                # only the new position
            self.assertTrue(set(first.splitlines()) < set(second.splitlines()))
            shutil.rmtree(f"{d}/store")                # an earlier run, no store
            r = subprocess.run(cmd + ["--mode", "seed"], capture_output=True,
                               text=True, env=env)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("seeded 3 positions", r.stdout)
            self.assertEqual(project(), second)
            self.assertEqual(sent(), 3)                # served from the seed


class SvMatrixReader(unittest.TestCase):
    """SV genotyping keeps one sample's column, not the whole matrix."""

    def test_keeps_only_its_own_sample(self):
        g = load("p5_sv_genotype", "bin/p5_sv_genotype.py")
        with tempfile.TemporaryDirectory() as d:
            m = write(os.path.join(d, "m.tsv"),
                      "key\tsvtype\th37rv_pos\tsvlen\tS1\tS2\tS3\n"
                      "sv:DEL:10:100\tDEL\t10\t100\tALT\tNOCALL\tABSENT\n")
            rows, samples = g.read_sv_matrix(m, keep_sample="S2")
            self.assertEqual(samples, ["S1", "S2", "S3"])
            self.assertEqual(rows[0], {"key": "sv:DEL:10:100", "svtype": "DEL",
                                       "h37rv_pos": "10", "svlen": "100", "S2": "NOCALL"})
            rows, _ = g.read_sv_matrix(m)                 # the probe step
            self.assertNotIn("S1", rows[0])
            rows, _ = g.read_sv_matrix(m, keep_sample="S9")   # not in the matrix
            self.assertEqual(rows[0].get("S9", "NOCALL"), "NOCALL")


class AlignmentArchive(unittest.TestCase):
    """bin/archive_alignments.sh: CRAM archive, verified, and back to BAM."""

    REF = ("ACGTTGCAAGGCTTACCGATCGATTACGGATCCATGCAAGTCGATCGTAGCTAGCTTAGG"
           "CATCGATCGGATCGAATTCGGCTAGCTAGGATCCGATAGC")          # 100 bp

    def setup_sample(self, d, header_len=100):
        st = os.environ.get("MTB_SAMTOOLS", "")
        if not (st and os.path.exists(st)):
            self.skipTest("MTB_SAMTOOLS not set; source config/project_env.sh")
        ref = self.REF
        self.assertEqual(len(ref), 100)
        for sub in ("build/refs", "clean", "p2w", "p1i", "p1w", "p1g"):
            os.makedirs(os.path.join(d, sub))
        for fa in (f"{d}/build/refs/R.fasta", f"{d}/clean/R.isclean.fasta",
                   f"{d}/clean/H37Rv.isclean.fasta"):
            write(fa, f">c\n{ref}\n")
        sam = (f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:{header_len}\n"
               "@RG\tID:S\tSM:S\n"
               + "".join(f"r{i}\t0\tc\t{1 + 10 * i}\t60\t20M\t*\t0\t0\t"
                         f"{ref[10 * i:10 * i + 20]}\t{'F' * 20}\tRG:Z:S\tNM:i:0\tMD:Z:20\n"
                         for i in range(5)))
        write(f"{d}/x.sam", sam)
        for bam in (f"{d}/p2w/S.bam", f"{d}/p1i/S.isclean.bam", f"{d}/p1w/S.h37rv.bam",
                    f"{d}/p1g/S.isclean.bam"):
            subprocess.run([st, "view", "-b", "-o", bam, f"{d}/x.sam"], check=True,
                           capture_output=True)
            subprocess.run([st, "index", bam], check=True)
        write(f"{d}/refmap.tsv", "sample\ta\tb\tc\treference\nS\t.\t.\t.\tR\n")
        env = dict(os.environ, MTB_BUILD_DIR=f"{d}/build", REFMAP=f"{d}/refmap.tsv",
                   P1WORK=f"{d}/p1w", P2WORK=f"{d}/p2w", P1IDIR=f"{d}/p1i",
                   P1GDIR=f"{d}/p1g", CLEANDIR=f"{d}/clean",
                   H37CLEAN=f"{d}/clean/H37Rv.isclean.fasta",
                   ARCHIVE_DELETE_BAM="1", ARCHIVE_DROP_UNUSED="1")
        return st, env

    def run_step(self, step, env):
        return subprocess.run(["bash", "bin/archive_alignments.sh", step, "S"],
                              capture_output=True, text=True, env=env)

    def test_archive_then_restore_round_trips(self):
        with tempfile.TemporaryDirectory() as d:
            st, env = self.setup_sample(d)
            r = self.run_step("--archive", env)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            self.assertIn("verified", r.stdout)
            for gone in ("p2w/S.bam", "p1i/S.isclean.bam", "p1w/S.h37rv.bam",
                         "p1g/S.isclean.bam"):
                self.assertFalse(os.path.exists(f"{d}/{gone}"), gone)
            self.assertTrue(os.path.exists(f"{d}/p2w/S.archive.cram"))
            # the p1g alignment is read by the SV two-frame step: archived, not dropped
            self.assertTrue(os.path.exists(f"{d}/p1g/S.isclean.archive.cram"))
            r = self.run_step("--restore", env)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            view = lambda f: subprocess.run([st, "view", f], capture_output=True,
                                            text=True, check=True).stdout
            norm = lambda t: sorted("\t".join(l.split("\t")[:11] + sorted(l.split("\t")[11:]))
                                    for l in t.splitlines())
            self.assertEqual(norm(view(f"{d}/p2w/S.bam")), norm(view(f"{d}/x.sam")))

    def test_mismatch_keeps_every_bam(self):
        # a header that disagrees with the reference: CRAM cannot reproduce it
        with tempfile.TemporaryDirectory() as d:
            st, env = self.setup_sample(d, header_len=101)
            r = self.run_step("--archive", env)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("does not reproduce", r.stderr)
            for kept in ("p2w/S.bam", "p1i/S.isclean.bam", "p1w/S.h37rv.bam",
                         "p1g/S.isclean.bam"):
                self.assertTrue(os.path.exists(f"{d}/{kept}"), kept)


class DrArrayRescue(unittest.TestCase):
    """IS6110 copies in the DR array: junction reads scatter over identical
    repeats with MAPQ 0, so no stack passes. The rescue merges them into one
    locus-level site and leaves every already-called sample unchanged."""

    DR = "GTCGTCAGACCCAAAACCCCGAGAGGGGACGGAAAC"
    HDR = ["sample", "clean_pos", "orig_pos", "reads", "reads_q", "positions", "span",
           "sa_mapq_max", "sa_mapq_mean", "el_start", "el_end", "el_internal",
           "chr_start", "chr_end", "both_el_termini", "both_chr_sides", "fwd", "rev"]

    def stack(self, s, pos, reads, rq, es, ee, mq=0):
        return [s, pos, pos, reads, rq, 1, 0, mq, float(mq), es, ee, 0, 1, 0,
                int(es > 0 and ee > 0), 0, reads, 0]

    def test_rescue(self):
        import random
        rng = random.Random(3)
        uniq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
        left = uniq(5000)
        array = "".join(self.DR + uniq(36) for _ in range(20))
        seq = left + array + uniq(5000)
        lo = len(left) + 1                      # 1-based start of the array
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(f"{d}/clean"); os.makedirs(f"{d}/p1i")
            write(f"{d}/clean/R.isclean.fasta", ">c\n" + seq + "\n")
            write(f"{d}/refmap.tsv", "sample\treference\nA\tR\nB\tR\nC\tR\n")
            tables = {
                # A: scattered over repeats, none passes -> rescued, two-sided
                "A": [self.stack("A", lo + 72 * i, 4, 0, 1 if i % 2 else 0, 0 if i % 2 else 1)
                      for i in range(6)] + [self.stack("A", 900, 30, 30, 15, 15, 60)],
                # B: a DR stack already passes -> unchanged
                "B": [self.stack("B", lo + 100, 25, 25, 12, 13, 60),
                      self.stack("B", lo + 300, 3, 0, 1, 0)],
                # C: nothing in the array -> unchanged
                "C": [self.stack("C", 900, 30, 30, 15, 15, 60)],
            }
            for smp, rows in tables.items():
                with open(f"{d}/p1i/{smp}.elstacks.tsv", "w") as fh:
                    fh.write("\t".join(self.HDR) + "\r\n")
                    for r in rows:
                        fh.write("\t".join(map(str, r)) + "\r\n")
            cmd = [sys.executable, "is6110/bin/is6110_repeat_rescue.py", "--dir", f"{d}/p1i",
                   "--refmap", f"{d}/refmap.tsv", "--clean-dir", f"{d}/clean"]
            for _ in range(2):                  # the second run must be a no-op
                r = subprocess.run(cmd, capture_output=True, text=True)
                self.assertEqual(r.returncode, 0, r.stderr)
            out = {smp: rows_of(f"{d}/p1i/{smp}.elstacks.tsv") for smp in tables}
            dr = [x for x in out["A"] if x["repeat_locus"] == "DR"]
            self.assertEqual(len(dr), 1)
            self.assertEqual(int(dr[0]["clean_pos"]), lo)
            self.assertEqual(int(dr[0]["reads_q"]), 24)
            self.assertEqual(dr[0]["both_el_termini"], "1")
            self.assertEqual(len(out["A"]), 2)  # the unique-region stack is kept
            self.assertEqual(len(out["B"]), 2)
            self.assertFalse(any(x["repeat_locus"] for x in out["B"] + out["C"]))
            self.assertTrue(os.path.exists(f"{d}/p1i/A.elstacks.raw.tsv"))

    def test_two_clusters(self):
        """Two DR clusters: rescue one only when the other holds no weak
        evidence, so one copy's ambiguous reads never become two calls."""
        import random
        rng = random.Random(5)
        uniq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
        arr = lambda: "".join(self.DR + uniq(36) for _ in range(8))
        u1, a1, u2, a2 = uniq(5000), arr(), uniq(20000), arr()
        seq = u1 + a1 + u2 + a2 + uniq(5000)
        c1 = len(u1) + 1
        c2 = c1 + len(a1) + len(u2)
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(f"{d}/clean"); os.makedirs(f"{d}/p1i")
            write(f"{d}/clean/R.isclean.fasta", ">c\n" + seq + "\n")
            write(f"{d}/refmap.tsv", "sample\treference\nA\tR\nB\tR\n")
            tables = {
                "A": [self.stack("A", c1 + 72 * i, 4, 0, 1, 1) for i in range(3)],
                "B": [self.stack("B", c1 + 72 * i, 4, 0, 1, 1) for i in range(3)]
                     + [self.stack("B", c2 + 72 * i, 4, 0, 1, 1) for i in range(3)],
            }
            for smp, rows in tables.items():
                with open(f"{d}/p1i/{smp}.elstacks.tsv", "w") as fh:
                    fh.write("\t".join(self.HDR) + "\n")
                    for r in rows:
                        fh.write("\t".join(map(str, r)) + "\n")
            r = subprocess.run([sys.executable, "is6110/bin/is6110_repeat_rescue.py",
                                "--dir", f"{d}/p1i", "--refmap", f"{d}/refmap.tsv",
                                "--clean-dir", f"{d}/clean"], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            a = rows_of(f"{d}/p1i/A.elstacks.tsv")
            b = rows_of(f"{d}/p1i/B.elstacks.tsv")
            self.assertEqual([int(x["clean_pos"]) for x in a if x["repeat_locus"]], [c1])
            self.assertFalse(any(x["repeat_locus"] for x in b))
            self.assertEqual(len(b), 6)


class Is6110SeamAndKeys(unittest.TestCase):
    """One ref_shared rule everywhere (review 4.4); one key per insertion (4.5)."""

    def test_seam_slop_and_side(self):
        sm = load("is6110_seam", "is6110/bin/is6110_seam.py")
        with tempfile.TemporaryDirectory() as d:
            cm = write(f"{d}/R.crossmap.tsv", "orig_start\torig_end\tdeleted_len\t"
                       "clean_junction\tcum_deleted\n1001\t2355\t1355\t1000\t1355\n")
            s = sm.Seams(cm)
            self.assertEqual(s.nearest(1000), ("L", 1001, 2355, 0))
            self.assertEqual(s.nearest(2358), ("R", 1001, 2355, 2))
            self.assertIsNone(s.nearest(996))           # 4 bp off: not shared
            self.assertTrue(s.shared(1003))             # inside the duplication
            self.assertFalse(sm.Seams(f"{d}/missing.tsv").known())
        # promotion and the writer agree, by construction
        pr = load("is6110_promote_sites", "is6110/bin/is6110_promote_sites.py")
        self.assertEqual(pr.classify_site(2357, s, sm.SEAM_SLOP), pr.SHARED)
        self.assertEqual(pr.classify_site(2360, s, sm.SEAM_SLOP), pr.LACKING)

    def test_cluster_keys(self):
        w = load("is6110_write_vcf", "is6110/bin/is6110_write_vcf.py")
        obs = [("h", 100, "A", 1), ("h", 101, "B", 1), ("h", 101, "C", 1),
               ("h", 106, "D", 1),          # within 6 of the cluster start
               ("h", 107, "E", 1),          # 7 from the start: a new cluster
               ("h", 103, "A", 2),          # A already in the cluster: separate
               ("n", 101, "F", 1)]          # another axis is never merged in
        c = w.cluster_keys(obs, 6)
        self.assertEqual(c[("A", 1)], 101)      # most carriers wins
        self.assertEqual(c[("D", 1)], 101)
        self.assertEqual(c[("A", 2)], 103)
        self.assertEqual(c[("E", 1)], 103)      # joins the cluster 103 opened
        self.assertEqual(c[("F", 1)], 101)
        self.assertEqual(w.cluster_keys(obs, 0)[("D", 1)], 106)


class AncestralAtIngroup(unittest.TestCase):
    """Review 3.7: AA is the MTBC ancestor's state, the outgroup breaking a
    tie there; a site the root already resolved keeps its allele."""

    def test_ingroup_state(self):
        with tempfile.TemporaryDirectory() as d:
            write(f"{d}/t.nwk", "(O,(A,(B,C)));\n")
            # site 1: root ties {C,T}, ingroup is C. site 2: root resolves A.
            # site 3: ingroup ties {A,G} and the outgroup (T) cannot break it.
            write(f"{d}/a.fa", ">O\nTAT\n>A\nCAA\n>B\nCGG\n>C\nCGG\n")
            write(f"{d}/s.tsv", "chrom\tpos\tref\talt\nc\t1\tC\tT\nc\t2\tA\tG\nc\t3\tA\tG\n")
            r = subprocess.run([sys.executable, "bin/ancestral_alleles.py", "--tree",
                                f"{d}/t.nwk", "--alignment", f"{d}/a.fa", "--sites",
                                f"{d}/s.tsv", "--out", f"{d}/o.tsv", "--outgroup", "O"],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            got = [(x["AA"], x["flag"]) for x in rows_of(f"{d}/o.tsv")]
            self.assertEqual(got, [("C", ""), ("A", ""), (".", "TIED:AG")])


def rows_of(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class ProjectionStore(unittest.TestCase):
    """Found regenerating scale200: `add` failed for any store path with /../"""

    def test_add_with_dotdot_in_store_path(self):
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "a"))
            store = os.path.join(d, "a", "..", "store")        # contains /../
            res = write(os.path.join(d, "r.pos"),
                        "#h\nR#1#c,9,+\tH#1#h,19,+\t0\t+\t+\n")
            r = subprocess.run([sys.executable, "bin/proj_store.py", "add",
                                "--store", store, "--ref", "R", "--result", res],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            got = os.listdir(os.path.join(d, "store", "R"))
            self.assertEqual(len([f for f in got if f.endswith(".pos")]), 1)
            self.assertFalse([f for f in got if f.endswith(".tmp")])


class RunnerCohortTable(unittest.TestCase):
    """CODE_REVIEW 6.7 and 6.3: blank lines inflated N; --cohort was rejected."""

    def runner(self, table_text, *extra):
        d = tempfile.mkdtemp()
        tab = write(os.path.join(d, "cohort.tsv"), table_text)
        crams = write(os.path.join(d, "crams.tsv"), "sample\tcram_relpath\n")
        reg = write(os.path.join(d, "reg.tsv"),
                    f"t\t{tab}\t{crams}\t{d}/out\tp1\tnote\t{d}/work/\n")
        # a minimal completed build, so the runner gets past build selection
        b = os.path.join(d, "build", "b0")
        os.makedirs(os.path.join(b, "logs"))
        og = write(os.path.join(d, "g.og"), "x\n")
        write(os.path.join(b, "logs", "manifest.done"), "done\n")
        write(os.path.join(b, "build_info.tsv"), f"build_id\tb0\ngraph\t{og}\n")
        frames = write(os.path.join(d, "frames.tsv"), "accession\n")
        env = dict(os.environ, REGISTRY=reg, BUILD_ROOT=os.path.join(d, "build"),
                   MTB_GRAPH_FRAMES=frames)
        return subprocess.run(["bash", "bin/refbias_run.sh", "--cohort", "t",
                               "--dry-run", *extra], capture_output=True,
                              text=True, env=env)

    def test_cohort_flag_accepted(self):
        r = self.runner("sample\tlineage\nA\t4\nB\t4\n")
        self.assertNotIn("unknown option", r.stderr)
        self.assertEqual(r.returncode, 0, r.stderr)

    def test_trailing_blank_line_not_counted(self):
        r = self.runner("sample\tlineage\nA\t4\nB\t4\n\n")
        self.assertRegex(r.stdout + r.stderr, r"cohort\s*:\s*t\s+2 isolates")

    def test_blank_line_between_rows_refused(self):
        r = self.runner("sample\tlineage\nA\t4\n\nB\t4\n")
        self.assertIn("between sample rows", r.stderr)


class ConfigSiteFile(unittest.TestCase):
    """CODE_REVIEW 2.1 and 2.2: site file order and working directory."""

    def test_site_file_moves_derived_paths_from_any_directory(self):
        d = tempfile.mkdtemp()
        cfg = os.path.join(d, "config")
        os.makedirs(cfg)
        with open(os.path.join(ROOT, "config", "project_env.sh")) as src, \
                open(os.path.join(cfg, "project_env.sh"), "w") as dst:
            dst.write(src.read())
        write(os.path.join(cfg, "site.local.sh"), "export MTB_WORK=/elsewhere\n")
        r = subprocess.run(
            ["env", "-i", f"HOME={os.environ.get('HOME', '/tmp')}", "bash", "-c",
             f"cd / && source {cfg}/project_env.sh && echo $MTB_DATA"],
            capture_output=True, text=True)
        self.assertEqual(r.stdout.strip(), "/elsewhere/data")


def load_tests(loader, standard_tests, pattern):
    """Also run the audit regression tests (tests/test_audit_*.py), so one
    command covers every pinned defect."""
    here = os.path.dirname(os.path.abspath(__file__))
    for fn in sorted(os.listdir(here)):
        if fn.startswith("test_audit_") and fn.endswith(".py"):
            mod = load(fn[:-3], os.path.join("tests", fn))
            standard_tests.addTests(loader.loadTestsFromModule(mod))
    return standard_tests


if __name__ == "__main__":
    unittest.main(warnings="ignore")
