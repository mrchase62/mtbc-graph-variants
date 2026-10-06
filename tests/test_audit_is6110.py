#!/usr/bin/env python3
"""Regression tests for the IS6110 findings of the 2026-10-05 audit
(analysis/audit/p3_is6110_accessory.md, P3IS-1 and P3IS-3 to P3IS-7).

Small synthetic inputs, standard library plus the pipeline's own scripts:

    $MTB_PY tests/test_audit_is6110.py            # from the repository root
"""
import csv
import gzip
import importlib.util
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
import warnings

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


def rows_of(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def tsv(path, cols, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    return path


class WriterFixture:
    """Inputs for is6110_write_vcf.py: two isolates on two references.

    R1 has one excised copy, orig 1001-2355 (clean junction 1000), so a site
    at orig 1000 is ref_shared for A. R2 has none."""

    RECON = ["sample", "orig_pos", "geometry", "chrom_side", "reads_q", "span"]
    FLANK = ["sample", "r_pos", "verdict", "h37rv_pos", "node", "node_offset",
             "ismapper", "node_occ"]

    def build(self, d, sites, flank_cols=None):
        os.makedirs(f"{d}/refs"); os.makedirs(f"{d}/cm"); os.makedirs(f"{d}/gff")
        write(f"{d}/refmap.tsv", "sample\treference\nA\tR1\nB\tR2\n")
        for r in ("R1", "R2"):
            write(f"{d}/refs/{r}.fasta", f">{r}c\n" + "A" * 9000 + "\n")
        write(f"{d}/cm/R1.crossmap.tsv", "orig_start\torig_end\tdeleted_len\t"
              "clean_junction\tcum_deleted\n1001\t2355\t1355\t1000\t1355\n")
        write(f"{d}/cm/R2.crossmap.tsv", "orig_start\torig_end\tdeleted_len\t"
              "clean_junction\tcum_deleted\n")
        write(f"{d}/h37.fasta", ">NC_000962.3\n" + "C" * 3000 + "\n")
        tsv(f"{d}/recon.tsv", self.RECON,
            [dict(sample=s, orig_pos=p, geometry="tsd", chrom_side="agreed",
                  reads_q=20, span=3) for s, p, _ in sites])
        tsv(f"{d}/flank.tsv", flank_cols or self.FLANK,
            [dict(sample=s, r_pos=p, ismapper="", **f) for s, p, f in sites])
        return [sys.executable, "is6110/bin/is6110_write_vcf.py",
                "--reconcile", f"{d}/recon.tsv", "--flank", f"{d}/flank.tsv",
                "--refmap", f"{d}/refmap.tsv", "--refs", f"{d}/refs",
                "--crossmap-dir", f"{d}/cm", "--gff-dir", f"{d}/gff",
                "--h37rv", f"{d}/h37.fasta", "--outdir", f"{d}/vcf",
                "--keys-out", f"{d}/keys.tsv", "--build-id", "b"]


EMPTY = lambda pos: dict(verdict="placed_h37rv_empty", h37rv_pos=pos)
OCC = lambda pos: dict(verdict="placed_h37rv_occupied", h37rv_pos=pos)
NODE = lambda nid, occ: dict(verdict="flanks_disagree_gap", node=nid,
                             node_offset=0, node_occ=occ)


class Is6110CohortFrameState(unittest.TestCase, WriterFixture):
    """P3IS-1: a carrier's state in the cohort key table is stated in the
    key's frame, not in its matched reference's."""

    def run_writer(self, d, sites):
        r = subprocess.run(self.build(d, sites), capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        return {(x["sample"], int(x["r_pos"])): x for x in rows_of(f"{d}/keys.tsv")}

    def test_carrier_states_in_the_key_frame(self):
        with tempfile.TemporaryDirectory() as d:
            k = self.run_writer(d, [
                ("A", 1000, EMPTY(500)),     # shared with R1, H37Rv empty
                ("B", 1000, EMPTY(500)),     # R2 lacks it, H37Rv empty
                ("A", 4000, OCC(900)),       # H37Rv holds the element
                ("A", 1003, NODE(77, 1)),    # shared with R1, node frame
            ])
            self.assertEqual(k[("A", 1000)]["site_class"], "ref_shared")
            self.assertEqual(k[("A", 1000)]["state"], "ALT")     # was REF
            self.assertEqual(k[("B", 1000)]["state"], "ALT")
            self.assertEqual(k[("A", 1000)]["key"], k[("B", 1000)]["key"])
            self.assertEqual(k[("A", 4000)]["state"], "REF")     # was ALT
            self.assertEqual(k[("A", 1003)]["frame"], "node")
            self.assertEqual(k[("A", 1003)]["state"], "ALT")     # was REF
            # the derived H37Rv-frame file: the shared carrier at the empty
            # locus has its record; nothing at the occupied one
            h = [l.split("\t") for l in open(f"{d}/vcf/A.is6110.h37rv.vcf")
                 if not l.startswith("#")]
            self.assertEqual([x[1] for x in h], ["500"])
            # matched frame unchanged: the shared copy is REF there, no record
            m = [l.split("\t")[1] for l in open(f"{d}/vcf/A.is6110.vcf")
                 if not l.startswith("#")]
            self.assertEqual(m, ["4000"])

    def test_meinfo_polarity_unknown(self):
        """P3IS-6: the element strand is not determined, so it is not '+'."""
        with tempfile.TemporaryDirectory() as d:
            self.run_writer(d, [("B", 1000, EMPTY(500))])
            rec = [l for l in open(f"{d}/vcf/B.is6110.vcf") if not l.startswith("#")]
            self.assertIn("MEINFO=IS6110,1,1355,.;", rec[0])


class Is6110RepeatedNode(unittest.TestCase, WriterFixture):
    """P3IS-3 and P3IS-5: a node its carrier's path visits more than once is
    no key, so unrelated insertions are not merged into one record."""

    def test_repeated_node_is_not_a_key(self):
        with tempfile.TemporaryDirectory() as d:
            r = subprocess.run(self.build(d, [
                ("B", 3000, NODE(46966, 289)),
                ("B", 7000, NODE(46966, 289)),   # 4 kb away, same 1 bp node
                ("A", 4000, NODE(46966, 276)),
                ("A", 5000, NODE(12, 1))]), capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            k = rows_of(f"{d}/keys.tsv")
            rep = [x for x in k if x["frame"] == "repeat_node"]
            self.assertEqual(len(rep), 3)
            self.assertTrue(all(x["key"] == "" and x["node"] == "" for x in rep))
            self.assertTrue(all(x["node_placed"] == "46966:0" for x in rep))
            self.assertEqual([x["key"] for x in k if x["frame"] == "node"],
                             ["node:12:0"])
            # no (sample, key) pair twice: the old table had B twice on
            # node:46966:0, settled downstream by row order
            keyed = [(x["sample"], x["key"]) for x in k if x["key"]]
            self.assertEqual(len(keyed), len(set(keyed)))

    def test_table_without_occurrence_count_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            cmd = self.build(d, [("A", 5000, NODE(12, 1))],
                             flank_cols=self.FLANK[:-1])
            r = subprocess.run(cmd, capture_output=True, text=True)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("node_occ", r.stderr)

    def test_two_sites_of_one_isolate_on_one_key_stop_the_writer(self):
        with tempfile.TemporaryDirectory() as d:
            r = subprocess.run(self.build(d, [("B", 3000, EMPTY(500)),
                                              ("B", 3100, EMPTY(500))]),
                               capture_output=True, text=True)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("(sample, key) pairs", r.stderr)

    def test_occurrence_count_from_odgi(self):
        with tempfile.TemporaryDirectory() as d:
            # the module builds its frame table at import
            old = os.environ.get("MTB_GRAPH_FRAMES")
            os.environ["MTB_GRAPH_FRAMES"] = write(
                f"{d}/frames.tsv",
                "accession\tpanel_len\trefs_len\tstrand\toffset\tagree\n")
            try:
                ps = load("is6110_project_sites_t", "is6110/bin/is6110_project_sites.py")
            finally:
                if old is None:
                    os.environ.pop("MTB_GRAPH_FRAMES")
                else:
                    os.environ["MTB_GRAPH_FRAMES"] = old
            fake = write(f"{d}/odgi", "#!/usr/bin/env python3\n"
                         "print('path.name\\tpath.length\\tpath.step.count\\t"
                         "node.1\\tnode.2\\tnode.3')\n"
                         "print('R1#1#c\\t10\\t5\\t1\\t3\\t0')\n"
                         "print('R2#1#c\\t10\\t5\\t2\\t1\\t1')\n")
            os.chmod(fake, 0o755)
            got = ps.node_occurrences(fake, "g.og", {("R1#1#c", 2), ("R2#1#c", 2),
                                                     ("R1#1#c", 1)}, 1)
        self.assertEqual(got, {("R1#1#c", 2): 3, ("R2#1#c", 2): 1,
                               ("R1#1#c", 1): 1})


class Is6110P5MergeDuplicates(unittest.TestCase):
    """P3IS-5 in stage 1, and repeat_node rows (P3IS-3) not keyed there."""

    COLS = ["sample", "reference", "frame", "key", "h37rv_pos", "h37rv_state",
            "node", "state", "r_pos"]

    def run_merge(self, d, rows):
        write(f"{d}/refmap.tsv", "sample\treference\nA\tR\nB\tR\n")
        write(f"{d}/h37.fasta", ">NC_000962.3\n" + "C" * 3000 + "\n")
        tsv(f"{d}/keys.tsv", self.COLS, rows)
        return subprocess.run(
            [sys.executable, "is6110/bin/is6110_p5_merge.py", "--cohort-keys",
             f"{d}/keys.tsv", "--refmap", f"{d}/refmap.tsv", "--h37rv",
             f"{d}/h37.fasta", "--out-keys", f"{d}/k.tsv", "--out-states",
             f"{d}/s.tsv"], capture_output=True, text=True)

    def test_conflicting_duplicate_stops(self):
        with tempfile.TemporaryDirectory() as d:
            r = self.run_merge(d, [
                dict(sample="A", frame="node", key="node:5:0", node="5:0", state="ALT", r_pos=1),
                dict(sample="A", frame="node", key="node:5:0", node="5:0", state="REF", r_pos=2)])
            self.assertNotEqual(r.returncode, 0)      # was: last row won (REF)

    def test_repeat_node_row_not_keyed(self):
        with tempfile.TemporaryDirectory() as d:
            r = self.run_merge(d, [
                dict(sample="A", frame="h37rv", key="h37rv:10", h37rv_pos=10,
                     h37rv_state="empty", state="ALT", r_pos=1),
                dict(sample="B", frame="repeat_node", key="", node="", state="ALT", r_pos=2)])
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual([x["key"] for x in rows_of(f"{d}/k.tsv")], ["h37rv:10:C><INS>"])
            self.assertIn("repeated graph node", r.stdout)


class Is6110Stage2Occupied(unittest.TestCase):
    """P3IS-1, stage 2: at a locus H37Rv holds, a non-carrier with depth and no
    junction lacks the element; it is not written REF (matching H37Rv)."""

    def test_non_carrier_at_occupied_locus_is_nocall(self):
        with tempfile.TemporaryDirectory() as d:
            write(f"{d}/refmap.tsv", "sample\treference\nA\tR\nB\tR\n")
            tsv(f"{d}/keys.tsv", ["sample", "reference", "r_pos", "frame",
                                  "h37rv_pos", "h37rv_state", "node"],
                [dict(sample="A", reference="R", r_pos=100, frame="h37rv",
                      h37rv_pos=500, h37rv_state="empty"),
                 dict(sample="A", reference="R", r_pos=300, frame="h37rv",
                      h37rv_pos=900, h37rv_state="occupied")])
            os.makedirs(f"{d}/s1"); os.makedirs(f"{d}/w/proj"); os.makedirs(f"{d}/cm")
            write(f"{d}/s1/B.tsv", "sample\tkey\tstate\tallele\n"
                  "B\th37rv:500:C><INS>\tNOCALL\t\nB\th37rv:900:C><INS>\tNOCALL\t\n")
            write(f"{d}/w/proj/R.tsv", "short_key\tpos\tdist\n"
                  "h37rv:500:\t100\t0\nh37rv:900:\t300\t0\n")
            write(f"{d}/cm/R.crossmap.tsv", "orig_end\tcum_deleted\n")
            write(f"{d}/B.isclean.bam", "")
            fake = write(f"{d}/samtools", "#!/usr/bin/env python3\n"
                         "import sys\n"
                         "if sys.argv[1] == 'idxstats':\n"
                         "    print('R_isclean\\t1000\\t50\\t0')\n"
                         "else:\n"
                         "    for l in open(sys.argv[sys.argv.index('-b') + 1]):\n"
                         "        c, s, e = l.split()\n"
                         "        print(f'{c}\\t{e}\\t30')\n")
            os.chmod(fake, 0o755)
            r = subprocess.run(
                [sys.executable, "is6110/bin/is6110_p5_stage2.py", "--mode", "sample",
                 "--sample", "B", "--refmap", f"{d}/refmap.tsv",
                 "--cohort-keys", f"{d}/keys.tsv", "--stage1-dir", f"{d}/s1",
                 "--isclean-dir", d, "--crossmap-dir", f"{d}/cm",
                 "--samtools", fake, "--workdir", f"{d}/w", "--out", f"{d}/o.tsv"],
                capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            st = {x["key"]: x["state"] for x in rows_of(f"{d}/w/samples/B.tsv")}
        self.assertEqual(st["h37rv:500:C><INS>"], "REF")
        self.assertEqual(st["h37rv:900:C><INS>"], "NOCALL")      # was REF


class Is6110MergedVcf(unittest.TestCase):
    """P3IS-5 and P3IS-6 in bin/merge_cohort_vcf.py's IS6110 block."""

    COLS = ["sample", "key", "frame", "h37rv_pos", "h37rv_state", "node",
            "state", "evidence", "site_class"]

    def setUp(self):
        self.bgzip = os.environ.get("MTB_BGZIP", "")
        _t = os.path.join(os.path.dirname(self.bgzip), "tabix") if self.bgzip else ""
        self.tabix = os.environ.get("MTB_TABIX", _t)
        if not (self.bgzip and os.path.exists(self.bgzip) and os.path.exists(self.tabix)):
            self.skipTest("MTB_BGZIP not set; source config/project_env.sh")

    def merge(self, d, rows):
        write(f"{d}/m.tsv", "key\tframe\tA\tB\tC\n")
        tsv(f"{d}/is.tsv", self.COLS, rows)
        return subprocess.run(
            [sys.executable, "bin/merge_cohort_vcf.py", "--matrix", f"{d}/m.tsv",
             "--is6110-keys", f"{d}/is.tsv", "--is6110-states", f"{d}/none.tsv",
             "--cohort-name", "t", "--ancestral", "", "--bgzip", self.bgzip,
             "--tabix", self.tabix, "--out", f"{d}/o.vcf.gz"],
            capture_output=True, text=True)

    def recs(self, d):
        with gzip.open(f"{d}/o.vcf.gz", "rt") as fh:
            return {l.split("\t")[2]: l for l in fh if not l.startswith("#")}

    def test_key_level_info_and_duplicates(self):
        h = dict(frame="h37rv", h37rv_pos=500, h37rv_state="empty", key="h37rv:500")
        with tempfile.TemporaryDirectory() as d:
            r = self.merge(d, [
                dict(h, sample="A", state="ALT", evidence="one_sided", site_class="ref_shared"),
                dict(h, sample="B", state="ALT", evidence="two_sided", site_class="ref_shared"),
                dict(h, sample="B", state="ALT", evidence="two_sided", site_class="ref_shared"),
                dict(sample="C", frame="repeat_node", key="", state="ALT",
                     evidence="two_sided", site_class="ref_lacking")])
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            recs = self.recs(d)
        self.assertEqual(list(recs), ["h37rv:500"])          # repeat_node: no record
        rec = recs["h37rv:500"]
        self.assertIn("EVIDENCE=mixed", rec)                 # was the first row's
        self.assertIn("SITECLASS=ref_shared", rec)
        self.assertIn("MEINFO=IS6110,1,1355,.", rec)
        self.assertTrue(rec.rstrip("\n").endswith("1:ALT\t1:ALT\t.:NOCALL"))

    def test_conflicting_duplicate_stops(self):
        h = dict(frame="h37rv", h37rv_pos=500, h37rv_state="empty", key="h37rv:500",
                 evidence="two_sided", site_class="ref_lacking")
        with tempfile.TemporaryDirectory() as d:
            r = self.merge(d, [dict(h, sample="A", state="ALT"),
                               dict(h, sample="A", state="REF")])
        self.assertNotEqual(r.returncode, 0)                 # was: last row won


class Is6110RescueOffset(unittest.TestCase):
    """P3IS-4: the DR-array start converts through the reference's crossmap,
    not with the nearest stack's offset."""

    DR = "GTCGTCAGACCCAAAACCCCGAGAGGGGACGGAAAC"
    HDR = ["sample", "clean_pos", "orig_pos", "reads", "reads_q", "positions", "span",
           "sa_mapq_max", "sa_mapq_mean", "el_start", "el_end", "el_internal",
           "chr_start", "chr_end", "both_el_termini", "both_chr_sides", "fwd", "rev"]

    def test_excision_between_stack_and_array_start(self):
        import random
        rng = random.Random(7)
        uniq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
        left = uniq(5000)
        seq = left + "".join(self.DR + uniq(36) for _ in range(20)) + uniq(5000)
        lo = len(left) + 1
        cj = lo + 10          # the reference's own copy was excised after this base
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(f"{d}/clean"); os.makedirs(f"{d}/p1i")
            write(f"{d}/clean/R.isclean.fasta", ">c\n" + seq + "\n")
            write(f"{d}/clean/R.crossmap.tsv", "orig_start\torig_end\tdeleted_len\t"
                  f"clean_junction\tcum_deleted\n{cj + 1}\t{cj + 1355}\t1355\t{cj}\t1355\n")
            write(f"{d}/refmap.tsv", "sample\treference\nA\tR\n")
            with open(f"{d}/p1i/A.elstacks.tsv", "w") as fh:
                fh.write("\t".join(self.HDR) + "\n")
                for i in range(1, 6):            # all past the excision
                    p = lo + 72 * i
                    fh.write("\t".join(map(str, [
                        "A", p, p + 1355, 4, 0, 1, 0, 0, 0.0, 1, 1, 0, 1, 0,
                        1, 0, 4, 0])) + "\n")
            r = subprocess.run([sys.executable, "is6110/bin/is6110_repeat_rescue.py",
                                "--dir", f"{d}/p1i", "--refmap", f"{d}/refmap.tsv",
                                "--clean-dir", f"{d}/clean"], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            dr = [x for x in rows_of(f"{d}/p1i/A.elstacks.tsv") if x["repeat_locus"]]
        self.assertEqual(int(dr[0]["clean_pos"]), lo)
        self.assertEqual(int(dr[0]["orig_pos"]), lo)        # was lo + 1355


class P1ivReconcileCrossmaps(unittest.TestCase):
    """P3IS-7: p1iv's reconcile uses each isolate's own crossmap."""

    def test_reconcile_gets_crossmap_dir(self):
        text = open("bin/p1i_vcf.sh").read()
        call = text[text.index("is6110/bin/is6110_reconcile.py"):]
        call = call[:call.index("\n\n")]
        self.assertIn("--crossmap-dir", call)
        self.assertIn('--ismapper-dir ""', call)


if __name__ == "__main__":
    unittest.main(warnings="ignore")
