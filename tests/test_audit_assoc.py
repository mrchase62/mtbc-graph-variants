#!/usr/bin/env python3
"""Regression tests for the association-chain findings of the 2026-10-05 audit
(analysis/audit/assoc_tests.md, ASSOC-2 to ASSOC-14; ASSOC-1 is deferred).

Same pattern as tests/run_tests.py: unittest, small synthetic inputs, the
standard library plus the pipeline's own scripts, seconds to run.

    $MTB_PY tests/test_audit_assoc.py -v

The scan and the burden read their event branches through
mtbvartools.CallBytestream, which $MTB_PY does not carry; a stand-in that
serves the same `.calls.col` / `.calls.loc[key]` interface from a JSON file is
installed in sys.modules for these tests.
"""
import contextlib
import csv
import importlib.util
import io
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import types
import unittest
import warnings

import numpy as np

warnings.simplefilter("ignore", ResourceWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
sys.path.insert(0, os.path.join(ROOT, "assoc", "bin"))
sys.path.insert(0, os.path.join(ROOT, "bin"))


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, rel))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---- a stand-in for mtbvartools.CallBytestream ------------------------------
class _Loc:
    def __init__(self, rows):
        self.rows = rows

    def __getitem__(self, key):
        return np.asarray(self.rows[repr(tuple(key))])


class _Calls:
    def __init__(self, d):
        self.col = d["cols"]
        self.loc = _Loc(d["rows"])


class FakeCallBytestream:
    def __init__(self, path):
        with open(os.path.join(path, "fake.json")) as fh:
            self.calls = _Calls(json.load(fh))

    def close(self):
        pass


sys.modules["mtbvartools"] = types.SimpleNamespace(
    CallBytestream=FakeCallBytestream)


# ---- a synthetic event matrix ------------------------------------------------
def balanced(leaves, counter):
    if len(leaves) == 1:
        return f"{leaves[0]}:1"
    counter[0] += 1
    lab = f"n{counter[0]:05d}"
    h = len(leaves) // 2
    return (f"({balanced(leaves[:h], counter)},{balanced(leaves[h:], counter)})"
            f"{lab}:1")


LEAVES = [f"L{i:02d}" for i in range(1, 17)]
CARRIERS = LEAVES[:4]
VCOLS = ["row_key", "id", "chrom", "pos", "ref", "alt", "class", "frame",
         "region", "acc_locus", "n_applicable", "n_gain", "n_loss", "n_undet",
         "n_derived_leaves"]


def make_events(d, variants):
    """`variants`: list of (row_key, id, pos, class, region, gain leaves/nodes,
    n_undet). Writes labelled.nwk, variants.tsv and the stand-in bytestream."""
    ev = os.path.join(d, "events")
    os.makedirs(os.path.join(ev, "event"), exist_ok=True)
    nwk = balanced(LEAVES, [0])
    nwk = nwk[:nwk.rindex(":")] + ";"
    with open(os.path.join(ev, "labelled.nwk"), "w") as fh:
        fh.write(nwk + "\n")
    from sv_scatter import tree_frame
    T = tree_frame(os.path.join(ev, "labelled.nwk"))
    cols = T["labels"]
    rows = {}
    with open(os.path.join(ev, "variants.tsv"), "w") as fh:
        fh.write("\t".join(VCOLS) + "\n")
        for rk, vid, pos, cls, reg, gains, und in variants:
            k = rk.split("|")
            key = (int(k[0]) if k[0].isdigit() else k[0],) + tuple(k[1:])
            v = [1 if c in gains else 0 for c in cols]
            rows[repr(key)] = v
            fh.write("\t".join(str(x) for x in (
                rk, vid, "NC_000962.3", pos, k[1], k[2], cls,
                "h37rv" if cls != "sv" else "", reg, "", "", len(gains), 0,
                und, len(gains))) + "\n")
    with open(os.path.join(ev, "event", "fake.json"), "w") as fh:
        json.dump({"cols": cols, "rows": rows}, fh)
    with open(os.path.join(d, "pheno.txt"), "w") as fh:
        fh.write("\n".join(CARRIERS) + "\n")
    with open(os.path.join(d, "lin.tsv"), "w") as fh:
        fh.write("sample\tlineage\n" + "".join(
            f"{s}\tlineage{1 + i // 8}\n" for i, s in enumerate(LEAVES)))
    return ev


def run_main(mod, argv):
    old = sys.argv
    sys.argv = ["x"] + argv
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            mod.main()
    finally:
        sys.argv = old
    return buf.getvalue()


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# background variants: two gains each on non-carrier leaves
BACKGROUND = [(f"{500 + i}|A|G", f"v{i}", 500 + i, "small", "core",
               [LEAVES[4 + (2 * i) % 12], LEAVES[5 + (2 * i) % 12]], 0)
              for i in range(6)]


class ScanBase(unittest.TestCase):
    def setUp(self):
        self.scan = load("assoc_scan", "assoc/bin/assoc_scan.py")
        self.d = tempfile.mkdtemp()

    def scan_run(self, variants, *extra, lineages=True):
        ev = make_events(self.d, variants)
        out = os.path.join(self.d, "scan.tsv")
        args = ["--events", ev, "--phenotype", os.path.join(self.d, "pheno.txt"),
                "--genes", "", "--permutations", "400", "--min-pool", "1",
                "--out", out, *extra]
        if lineages:
            args += ["--lineages", os.path.join(self.d, "lin.tsv")]
        log = run_main(self.scan, args)
        return {r["row_key"]: r for r in read_tsv(out)}, log


class LeaveOneOut(ScanBase):
    """ASSOC-2: the region null's leave-one-out removed nothing."""

    def test_helper_is_a_multiset_subtraction(self):
        self.assertEqual(self.scan.leave_one_out([1, 1, 2], [1]), [1, 2])
        self.assertEqual(self.scan.leave_one_out([3, 1, 1, 2], [1, 2, 9]),
                         [3, 1])

    def test_scan_null_pool_excludes_own_branches(self):
        rows, _ = self.scan_run(BACKGROUND)
        # every variant has 2 gains; the core pool holds 6 x 2 = 12 branches
        for r in rows.values():
            self.assertEqual(int(r["null_pool"]), 10)

    def test_burden_null_pool_excludes_own_branches(self):
        bur = load("is6110_gene_burden", "assoc/bin/is6110_gene_burden.py")
        genes = write_genes(self.d)
        v = [("1050|A|G", "a", 1050, "small", "core", ["L05", "L06"], 0),
             ("1150|A|G", "b", 1150, "small", "core", ["L07", "L08"], 0),
             ("1160|A|G", "c", 1160, "small", "core", ["L09"], 0)]
        ev = make_events(self.d, v)
        out = os.path.join(self.d, "b.tsv")
        run_main(bur, ["--events", ev, "--phenotype",
                       os.path.join(self.d, "pheno.txt"), "--genes", genes,
                       "--cls", "small", "--min-pool", "1",
                       "--permutations", "200", "--out", out])
        got = {r["gene"]: int(r["null_pool"]) for r in read_tsv(out)}
        # G1 owns {L05, L06}, G2 owns {L07, L08, L09}; pool 5 branches
        self.assertEqual(got, {"G1": 3, "G2": 2})


class DedupeSuffixKey(ScanBase):
    """ASSOC-3: a `#2` record was looked up as the first record."""

    V = [("200|C|<DEL>,*", "svi:DEL:201:50", 200, "sv", "",
          ["L05", "L06", "L07"], 0),
         ("200|C|<DEL>,*|#2", "svi:DEL:201:80", 200, "sv", "",
          ["L01", "L02"], 0)]

    def test_event_key_keeps_suffix(self):
        self.assertEqual(self.scan.event_key("200|C|<DEL>,*|#2"),
                         (200, "C", "<DEL>,*", "#2"))
        self.assertEqual(self.scan.event_key("n12|A|G"), ("n12", "A", "G"))

    def test_scan_uses_own_branches(self):
        rows, _ = self.scan_run(BACKGROUND + self.V, "--sv-intervals", "none")
        self.assertEqual(rows["200|C|<DEL>,*|#2"]["gains"], "2")
        self.assertEqual(rows["200|C|<DEL>,*|#2"]["obs"], "1.0000")
        self.assertEqual(rows["200|C|<DEL>,*"]["gains"], "3")

    def test_burden_uses_own_branches(self):
        bur = load("is6110_gene_burden", "assoc/bin/is6110_gene_burden.py")
        genes = write_genes(self.d)
        # the first deletion covers 1001-1010 (G1), the second 2001-2010 (G3);
        # both share the anchor key 1000|C|<DEL>,*
        v = [("1000|C|<DEL>,*", "svi:DEL:1001:10", 1000, "sv", "",
              ["L05", "L06"], 0),
             ("1000|C|<DEL>,*|#2", "svi:DEL:2001:10", 1000, "sv", "",
              ["L01", "L02"], 0)]
        ev = make_events(self.d, v)
        out = os.path.join(self.d, "b.tsv")
        run_main(bur, ["--events", ev, "--phenotype",
                       os.path.join(self.d, "pheno.txt"), "--genes", genes,
                       "--cls", "sv", "--min-origins", "1",
                       "--permutations", "100", "--out", out])
        got = {r["gene"]: r["obs"] for r in read_tsv(out)}
        self.assertEqual(got.get("G3"), "1.0000")   # its own carrier branches


class SvTiers(ScanBase):
    """ASSOC-4: svi: deletions scanned without (or with a stale) tier table."""

    V = DedupeSuffixKey.V

    def test_refuses_svi_without_tier_table(self):
        with self.assertRaises(SystemExit):
            self.scan_run(BACKGROUND + self.V)

    def test_refuses_table_missing_an_id(self):
        t = os.path.join(self.d, "tiers.tsv")
        with open(t, "w") as fh:
            fh.write("interval\tevidence_tier\nsvi:DEL:201:50\tE1_coherent\n")
        with self.assertRaises(SystemExit):
            self.scan_run(BACKGROUND + self.V, "--sv-intervals", t)

    def test_tiered_region(self):
        t = os.path.join(self.d, "tiers.tsv")
        with open(t, "w") as fh:
            fh.write("interval\tevidence_tier\nsvi:DEL:201:50\tE1_coherent\n"
                     "svi:DEL:201:80\tE3_scattered\n")
        rows, _ = self.scan_run(BACKGROUND + self.V, "--sv-intervals", t)
        self.assertEqual(rows["200|C|<DEL>,*"]["region"], "svi:E1_coherent")
        self.assertEqual(rows["200|C|<DEL>,*|#2"]["region"], "svi:E3_scattered")


def write_genes(d):
    """snpEff-dump genes; starts are 0-based. G1 1001-1100 (+), G2 1101-1200
    (+), G3 2001-2100 (-)."""
    p = os.path.join(d, "genes.txt")
    with open(p, "w") as fh:
        fh.write("chr\tstart\tend\tstrand\ttype\tid\tgeneName\tgeneId\tn\n")
        fh.write("NC_000962\t0\t5000\t+1\tChromosome\tNC\t\t\t\n")
        for s, e, st, n in ((1000, 1100, "+1", "G1"), (1100, 1200, "+1", "G2"),
                            (2000, 2100, "-1", "G3")):
            fh.write(f"NC_000962\t{s}\t{e}\t{st}\tGene\t{n}\t{n}\t{n}\t1\n")
    return p


class BurdenBase(unittest.TestCase):
    def setUp(self):
        self.bur = load("is6110_gene_burden", "assoc/bin/is6110_gene_burden.py")
        self.d = tempfile.mkdtemp()
        self.genes = write_genes(self.d)

    def units(self, variants, cls, *extra):
        ev = make_events(self.d, variants)
        out = os.path.join(self.d, "b.tsv")
        run_main(self.bur, ["--events", ev, "--phenotype",
                            os.path.join(self.d, "pheno.txt"),
                            "--genes", self.genes, "--cls", cls,
                            "--min-origins", "1", "--permutations", "100",
                            "--out", out, *extra])
        return {r["gene"]: r for r in read_tsv(out)}


class SvSpanCredit(BurdenBase):
    """ASSOC-5: a deletion was credited only to the unit at its anchor base."""

    def test_deletion_spanning_two_genes_credits_both(self):
        # anchor 1049; deletes 1050-1149, half of G1 and half of G2
        u = self.units([("1049|A|<DEL:-100>", "svi:DEL:1050:100", 1049, "sv",
                         "", ["L01", "L05"], 0)], "sv")
        self.assertIn("G1", u)
        self.assertIn("G2", u)

    def test_anchor_at_gene_end_is_not_credited_to_that_gene(self):
        # G3 (minus strand) ends at 2100; the deletion removes 2101-2130,
        # upstream of G3, and none of its body
        u = self.units([("2100|A|<DEL:-30>", "svi:DEL:2101:30", 2100, "sv",
                         "", ["L01", "L05"], 0)], "sv")
        self.assertNotIn("G3", u)
        self.assertIn("up:G3", u)

    def test_min_overlap_option(self):
        # removes 1095-1104: 6 bp of G1, 4 bp of G2
        u = self.units([("1094|A|<DEL:-10>", "svi:DEL:1095:10", 1094, "sv",
                         "", ["L01", "L05"], 0)], "sv", "--sv-min-overlap", "5")
        self.assertIn("G1", u)
        self.assertNotIn("G2", u)


class GeneCoordinates(BurdenBase):
    """ASSOC-7: snpEff-dump starts are 0-based."""

    def test_base_before_gene_is_promoter(self):
        # dump start 1000 -> G1 begins at 1001; 1000 is G1's promoter
        u = self.units([("1000|A|G", "s", 1000, "small", "core",
                         ["L01", "L05"], 0)], "small")
        self.assertIn("up:G1", u)
        self.assertNotIn("G1", u)

    def test_scan_genic_split(self):
        scan = load("assoc_scan", "assoc/bin/assoc_scan.py")
        d = self.d
        v = BACKGROUND + [("1000|A|G", "s", 1000, "small", "core",
                           ["L01", "L02"], 0)]
        ev = make_events(d, v)
        out = os.path.join(d, "scan.tsv")
        run_main(scan, ["--events", ev, "--phenotype",
                        os.path.join(d, "pheno.txt"), "--genes", self.genes,
                        "--lineages", os.path.join(d, "lin.tsv"),
                        "--permutations", "100", "--min-pool", "1",
                        "--out", out])
        r = {x["row_key"]: x for x in read_tsv(out)}
        self.assertEqual(r["1000|A|G"]["region"], "core:intergenic")


class BurdenFloor(BurdenBase):
    """ASSOC-8: the burdens admitted records below the scan's callability floor."""

    def test_uncallable_record_contributes_nothing(self):
        # 30 branches; 10 undetermined is 67% determined, below 80%
        u = self.units([("1050|A|G", "s", 1050, "small", "core",
                         ["L01", "L05"], 10),
                        ("1150|A|G", "t", 1150, "small", "core",
                         ["L02", "L06"], 1)], "small")
        self.assertNotIn("G1", u)
        self.assertIn("G2", u)


class SurvivorRule(ScanBase):
    """ASSOC-9: survival must mean q < 0.05 under all three nulls."""

    STRONG = ("100|A|G", "strong", 100, "small", "core", CARRIERS, 0)

    def survivors(self, log):
        m = re.search(r"the (\d+) that survive", log)
        return int(m.group(1)) if m else 0

    def test_missing_lineage_null_is_not_a_pass(self):
        rows, log = self.scan_run(BACKGROUND + [self.STRONG], lineages=False)
        r = rows["100|A|G"]
        self.assertLess(float(r["q_region"]), 0.05)    # passes region alone
        self.assertTrue(math.isnan(float(r["q_lineage"])))
        self.assertEqual(self.survivors(log), 0)


class SilentDefaults(ScanBase):
    """ASSOC-10: inputs that silently changed the test now stop the run."""

    def test_missing_lineage_file(self):
        with self.assertRaises(SystemExit):
            self.scan_run(BACKGROUND, "--lineages",
                          os.path.join(self.d, "nope.tsv"), lineages=False)

    def test_lineage_table_without_columns(self):
        bad = os.path.join(self.d, "bad.tsv")
        with open(bad, "w") as fh:
            fh.write("isolate\tclade\n" + "".join(f"{s}\tx\n" for s in LEAVES))
        with self.assertRaises(SystemExit):
            self.scan_run(BACKGROUND, "--lineages", bad, lineages=False)

    def test_lineage_table_for_another_tree(self):
        bad = os.path.join(self.d, "bad.tsv")
        with open(bad, "w") as fh:
            fh.write("sample\tlineage\n" + "".join(
                f"X{s}\tlineage1\n" for s in LEAVES))
        with self.assertRaises(SystemExit):
            self.scan_run(BACKGROUND, "--lineages", bad, lineages=False)

    def test_phenotype_name_not_a_tip(self):
        make_events(self.d, BACKGROUND)
        with open(os.path.join(self.d, "pheno.txt"), "a") as fh:
            fh.write("NOT_A_TIP\n")
        out = os.path.join(self.d, "scan.tsv")
        with self.assertRaises(SystemExit):
            run_main(self.scan, ["--events", os.path.join(self.d, "events"),
                                 "--phenotype",
                                 os.path.join(self.d, "pheno.txt"),
                                 "--genes", "", "--out", out])

    def test_missing_genes_file(self):
        with self.assertRaises(SystemExit):
            self.scan_run(BACKGROUND, "--genes",
                          os.path.join(self.d, "nope.txt"))

    def test_burden_refuses_missing_lineage_file(self):
        bur = load("is6110_gene_burden", "assoc/bin/is6110_gene_burden.py")
        ev = make_events(self.d, BACKGROUND)
        with self.assertRaises(SystemExit):
            run_main(bur, ["--events", ev, "--phenotype",
                           os.path.join(self.d, "pheno.txt"),
                           "--genes", write_genes(self.d), "--cls", "small",
                           "--lineages", os.path.join(self.d, "nope.tsv"),
                           "--out", os.path.join(self.d, "b.tsv")])


class AuditReadsGainsAndBurdens(unittest.TestCase):
    """ASSOC-3 and ASSOC-10: the audit compares scan gains with n_gain and
    reads the three burdens."""

    def setUp(self):
        self.au = load("audit_chain", "bin/audit_chain.py")
        self.d = tempfile.mkdtemp()

    def test_gains_mismatch_fails(self):
        scan = os.path.join(self.d, "scan.tsv")
        with open(scan, "w") as fh:
            fh.write("row_key\tgains\n1|A|G\t3\n1|A|G|#2\t3\n")
        fails = []
        self.au.check_gains(scan, [{"row_key": "1|A|G", "n_gain": "3"},
                                   {"row_key": "1|A|G|#2", "n_gain": "2"}],
                            fails)
        self.assertEqual(len(fails), 1)

    def test_burden_without_lineage_null_fails(self):
        p = os.path.join(self.d, "small_gene.tsv")
        with open(p, "w") as fh:
            fh.write("gene\tp_lineage\nrpoB\tnan\n")
        fails, warns = [], []
        self.au.check_burdens([p, os.path.join(self.d, "sv_gene.tsv")],
                              fails, warns)
        self.assertEqual(len(fails), 2)            # blank null, absent file


class ChainScript(unittest.TestCase):
    """ASSOC-4 and ASSOC-6: the chain runs only from runroot, with this
    script's own siblings, and passes a tier table rebuilt from the cohort's
    current catalogue."""

    CHAIN = os.path.join(ROOT, "assoc", "bin", "cohort_assoc_tail.sh")

    def test_refuses_a_directory_whose_assoc_bin_is_not_its_own(self):
        d = tempfile.mkdtemp()
        os.makedirs(os.path.join(d, "assoc", "bin"))     # a stale copy's place
        os.symlink(os.path.join(ROOT, "bin"), os.path.join(d, "bin"))
        r = subprocess.run(["bash", self.CHAIN, "t", "p.txt"], cwd=d,
                           capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("run this from runroot", r.stderr)

    def test_retiers_and_passes_sv_intervals(self):
        d = tempfile.mkdtemp()
        os.makedirs(os.path.join(d, "assoc"))
        os.symlink(os.path.join(ROOT, "assoc", "bin"),
                   os.path.join(d, "assoc", "bin"))
        os.symlink(os.path.join(ROOT, "bin"), os.path.join(d, "bin"))
        for f in ("refbias/t/p5/merged.vcf.gz", "refbias/t/p5/sv_intervals.tsv",
                  "refbias/t/p5/svgt_iv_states.tsv", "refbias/t/p1/refmap.tsv",
                  "refbias/t.phenotype.tsv", "data/trees/t.snps.fasta",
                  "data/trees/t.combined.fasta",
                  "data/trees/t.combined.rooted.nwk", "data/trees/t.rooted.nwk",
                  "assoc/t/events/summary.txt", "p.txt"):
            os.makedirs(os.path.dirname(os.path.join(d, f)) or d,
                        exist_ok=True)
            with open(os.path.join(d, f), "w") as fh:
                fh.write("x\n")
        # rerun safety (audit TP-4): the merged VCF names its build, the build
        # holds the panel assets, and each existing product records the VCF
        # and build it was made from, so steps 1-4 are skipped as done
        import gzip, hashlib
        with gzip.open(os.path.join(d, "refbias/t/p5/merged.vcf.gz"), "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n##MTB_graph_build=b0\n"
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        b = os.path.join(d, "refbias", "build", "b0")
        os.makedirs(os.path.join(b, "assets"))
        for f, txt in (("build_info.tsv", "build_id\tb0\n"),
                       ("assets/panel_polarity.tsv", "chrom\n"),
                       ("assets/graph_collapsed.vcf.gz", "x\n")):
            with open(os.path.join(b, f), "w") as fh:
                fh.write(txt)

        def sha16(p):
            return hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]
        prov = (f"build_id\tb0\nvcf_sha\t"
                f"{sha16(os.path.join(d, 'refbias/t/p5/merged.vcf.gz'))}\n")
        for f in ("data/trees/t.snps.fasta", "data/trees/t.combined.fasta",
                  "data/trees/t.combined.rooted.nwk", "data/trees/t.rooted.nwk"):
            with open(os.path.join(d, f + ".prov"), "w") as fh:
                fh.write(prov)
        with open(os.path.join(d, "assoc/t/events/summary.txt.prov"), "w") as fh:
            fh.write(prov + "polarity_sha\t" + sha16(
                os.path.join(b, "assets/panel_polarity.tsv")) + "\n")
        log = os.path.join(d, "calls.log")
        stub = os.path.join(d, "stub.sh")
        with open(stub, "w") as fh:
            fh.write(f'#!/bin/bash\necho "$@" >> {log}\n')
        os.chmod(stub, 0o755)
        env = dict(os.environ, MTB_PY=stub, MTB_PY_VT=stub)
        r = subprocess.run(["bash", self.CHAIN, "t", "p.txt"], cwd=d,
                           capture_output=True, text=True, env=env)
        self.assertEqual(r.returncode, 0, r.stderr)
        calls = open(log).read().splitlines()
        retier = [i for i, c in enumerate(calls) if "retier_intervals.py" in c]
        scan = [i for i, c in enumerate(calls) if "assoc_scan.py" in c]
        self.assertTrue(retier and scan and retier[0] < scan[0], calls)
        self.assertIn("--intervals refbias/t/p5/sv_intervals.tsv",
                      calls[retier[0]])
        self.assertIn("--sv-intervals assoc/t/sv_intervals.retiered.tsv",
                      calls[scan[0]])
        # every sibling is this repository's copy
        for c in calls:
            m = re.match(r"(\S+\.py)", c)
            if m:
                self.assertTrue(os.path.realpath(m.group(1)).startswith(
                    os.path.realpath(ROOT) + os.sep), c)


class TreeLabels(unittest.TestCase):
    """ASSOC-14: tree_frame collapsed duplicate node labels silently."""

    def test_duplicate_internal_labels_refused(self):
        from sv_scatter import tree_frame
        d = tempfile.mkdtemp()
        p = os.path.join(d, "t.nwk")
        with open(p, "w") as fh:
            fh.write("((A:1,B:1)x:1,(C:1,D:1)x:1)r;\n")
        with self.assertRaises(SystemExit):
            tree_frame(p)

    def test_unlabelled_internal_nodes_refused(self):
        from sv_scatter import tree_frame
        d = tempfile.mkdtemp()
        p = os.path.join(d, "t.nwk")
        with open(p, "w") as fh:
            fh.write("((A:1,B:1):1,(C:1,D:1):1);\n")
        with self.assertRaises(SystemExit):
            tree_frame(p)


if __name__ == "__main__":
    unittest.main(warnings="ignore")
