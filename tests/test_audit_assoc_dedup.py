#!/usr/bin/env python3
"""Regression tests for the two association-step options of 2026-10-16 (user:
"add both options to the association step"):

  1  the identical-pattern collapse (assoc_scan.py, default on,
     --no-collapse-identical): records whose state rows are identical at every
     node, with the same applicability, are one test unit -- tested once and
     counted once in BH -- and every member keeps its row, carrying the
     unit's result and its members.
  2  the IS6110 sensitivity scan (--exclude-acc-is6110; the tail's
     ASSOC_IS6110_SENSITIVITY, default on): the scan again without the
     accessory_presence records flagged ACCIS6110, and scan_compare.py's table
     of associations found by one and not the other.

Same stand-in bytestream and synthetic tree as tests/test_audit_assoc.py.

    $MTB_PY tests/test_audit_assoc_dedup.py -v
"""
import gzip
import math
import os
import subprocess
import sys
import tempfile
import unittest
import warnings

warnings.simplefilter("ignore", ResourceWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tests"))
import test_audit_assoc as ta          # noqa: E402  (installs the stand-in)

BACKGROUND, CARRIERS, LEAVES = ta.BACKGROUND, ta.CARRIERS, ta.LEAVES
# one event written twice: an IS6110 insertion and the accessory record of the
# catalogue locus that is the copy -- the same carriers, the same branches
INS = ("300|A|<INS>", "is:300", 300, "is6110", "is6110", CARRIERS, 0)
COPY = ("300|A|<INS:acc:ACC_0000300>", "acc:ACC_0000300", 300,
        "accessory_presence", "accessory_presence", CARRIERS, 0)
# a hierarchical case: an insertion inside that locus, in half its carriers
NESTED = ("301|A|<INS>", "is:301", 301, "is6110", "is6110", CARRIERS[:2], 0)


def read_text(path):
    with open(path) as fh:
        return fh.read()


def read_kv(path):
    return dict(l.split("\t", 1) for l in read_text(path).splitlines()
                if not l.startswith("key\t"))


class CollapseUnits(unittest.TestCase):
    def setUp(self):
        self.scan = ta.load("assoc_scan", "assoc/bin/assoc_scan.py")

    def test_representative_by_class_rank_then_order(self):
        rows = [{"class": "accessory_presence", "k": "a"},
                {"class": "small", "k": "b"},
                {"class": "is6110", "k": "a"},
                {"class": "sv", "k": "c"}]
        u = self.scan.collapse_units(rows, lambda r: r["k"])
        # units in first-member order; the IS6110 record represents its copy
        self.assertEqual(u, [[2, 0], [1], [3]])

    def test_none_key_stays_alone(self):
        rows = [{"class": "small"}, {"class": "small"}]
        self.assertEqual(self.scan.collapse_units(rows, lambda r: None),
                         [[0], [1]])


class IdenticalPatternCollapse(ta.ScanBase):
    def run_both(self, variants, *extra):
        on, log = self.scan_run(variants, *extra)
        off, _ = self.scan_run(variants, "--no-collapse-identical", *extra)
        return on, off, log

    def test_duplicate_is_one_unit_with_both_rows(self):
        on, off, log = self.run_both(BACKGROUND + [INS, COPY])
        a, b = on[INS[0]], on[COPY[0]]
        self.assertEqual(a["unit_id"], b["unit_id"])
        self.assertEqual(a["unit_id"], f"u:{INS[0]}")
        self.assertEqual((a["unit_size"], a["unit_rep"], b["unit_rep"]),
                         ("2", "1", "0"))
        for f in ("p_branch", "q_branch", "p_region", "q_region",
                  "p_lineage", "q_lineage"):
            self.assertEqual(a[f], b[f], f)
        # each member keeps its own class and region; the unit's BH family
        # is the representative's, and the unit is flagged as spanning two
        self.assertEqual(b["cls"], "accessory_presence")
        self.assertEqual(b["region"], "accessory_presence")
        self.assertEqual(b["unit_region"], "is6110")
        self.assertEqual(a["unit_regions"], "is6110,accessory_presence")
        self.assertEqual(a["unit_members"], "is:300;acc:ACC_0000300")
        self.assertIn("8 records -> 7 test units (1 tests removed", log)
        # off: one test per record, as before
        self.assertNotEqual(off[INS[0]]["unit_id"], off[COPY[0]]["unit_id"])
        self.assertEqual({r["unit_size"] for r in off.values()}, {"1"})

    def test_bh_family_counts_each_unit_once(self):
        dup = ("100|A|G", "s1", 100, "small", "core", CARRIERS, 0)
        dup2 = ("101|A|G", "s2", 101, "small", "core", CARRIERS, 0)
        on, off, _ = self.run_both(BACKGROUND + [dup, dup2])
        scan = self.scan
        for rows, n in ((on, 7), (off, 8)):
            reps = [r for r in rows.values() if r["unit_rep"] == "1"
                    and r["unit_region"] == "core"]
            self.assertEqual(len(reps), n)
            q = scan.bh([float(r["p_branch"]) for r in reps])
            for r, v in zip(reps, q):
                self.assertAlmostEqual(float(r["q_branch"]), v)

    def test_membership_and_run_record(self):
        out = os.path.join(self.d, "scan.tsv")
        self.scan_run(BACKGROUND + [INS, COPY])
        units = ta.read_tsv(os.path.join(self.d, "scan.units.tsv"))
        self.assertEqual(len(units), 8)                  # one line per record
        mem = [u for u in units if u["unit_id"] == f"u:{INS[0]}"]
        self.assertEqual(sorted((u["id"], u["class"], u["is_rep"]) for u in mem),
                         [("acc:ACC_0000300", "accessory_presence", "0"),
                          ("is:300", "is6110", "1")])
        kv = read_kv(os.path.join(self.d, "scan.run.tsv"))
        self.assertEqual(kv["collapse_identical"], "on")
        self.assertEqual((kv["records_tested"], kv["test_units"],
                          kv["tests_removed_by_collapse"]), ("8", "7", "1"))
        self.assertEqual(kv["exclude_acc_is6110"], "none")
        self.assertTrue(os.path.exists(out))

    def test_hierarchical_case_stays_separate(self):
        on, _, _ = self.run_both(BACKGROUND + [INS, COPY, NESTED])
        self.assertEqual(on[NESTED[0]]["unit_size"], "1")
        self.assertNotEqual(on[NESTED[0]]["unit_id"], on[INS[0]]["unit_id"])

    def test_no_duplicates_changes_nothing(self):
        on, off, _ = self.run_both(BACKGROUND + [INS, NESTED])
        for k in on:
            for f in ("p_branch", "q_branch", "p_region", "q_region",
                      "p_lineage", "q_lineage", "null_pool"):
                self.assertEqual(on[k][f], off[k][f], (k, f))

    def test_different_applicability_is_not_collapsed(self):
        # the same state row, but one record is conditional on a locus
        ev = ta.make_events(self.d, BACKGROUND + [INS, COPY])
        vt = os.path.join(ev, "variants.tsv")
        lines = open(vt).read().splitlines()
        hdr = lines[0].split("\t")
        ai = hdr.index("acc_locus")
        for j, l in enumerate(lines):
            f = l.split("\t")
            if f[0] == INS[0]:
                f[ai] = "ACC_0000999"
                lines[j] = "\t".join(f)
        open(vt, "w").write("\n".join(lines) + "\n")
        out = os.path.join(self.d, "s.tsv")
        ta.run_main(self.scan, ["--events", ev, "--phenotype",
                                os.path.join(self.d, "pheno.txt"),
                                "--genes", "", "--permutations", "200",
                                "--min-pool", "1", "--out", out])
        rows = {r["row_key"]: r for r in ta.read_tsv(out)}
        self.assertNotEqual(rows[INS[0]]["unit_id"], rows[COPY[0]]["unit_id"])

    def test_missing_state_stream_is_refused(self):
        ev = ta.make_events(self.d, BACKGROUND + [INS, COPY])
        import shutil
        shutil.rmtree(os.path.join(ev, "ancestor"))
        args = ["--events", ev, "--phenotype", os.path.join(self.d, "pheno.txt"),
                "--genes", "", "--permutations", "100", "--min-pool", "1",
                "--out", os.path.join(self.d, "s.tsv")]
        with self.assertRaises(SystemExit):
            ta.run_main(self.scan, args)
        ta.run_main(self.scan, args + ["--no-collapse-identical"])   # runs


class ExcludeAccIs6110(ta.ScanBase):
    def vcf(self, header=True, flag=("acc:ACC_0000300",)):
        p = os.path.join(self.d, "m.vcf.gz")
        with gzip.open(p, "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            if header:
                fh.write('##INFO=<ID=ACCIS6110,Number=0,Type=Flag,'
                         'Description="x">\n')
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            fh.write("NC_000962.3\t300\tis:300\tA\t<INS>\t.\tPASS\t"
                     "CLASS=is6110\n")
            for x in flag:
                fh.write(f"NC_000962.3\t300\t{x}\tA\t<INS>\t.\tPASS\t"
                         f"CLASS=accessory_presence;ACCIS6110\n")
        return p

    def test_flagged_records_are_excluded(self):
        rows, log = self.scan_run(BACKGROUND + [INS, COPY],
                                  "--exclude-acc-is6110", self.vcf())
        self.assertNotIn(COPY[0], rows)
        self.assertEqual(rows[INS[0]]["unit_size"], "1")
        self.assertIn("1 accessory_presence records flagged ACCIS6110 excluded",
                      log)
        kv = read_kv(os.path.join(self.d, "scan.run.tsv"))
        self.assertEqual(kv["acc_is6110_rows_excluded"], "1")

    def test_id_list_serves_too(self):
        p = os.path.join(self.d, "ids.txt")
        open(p, "w").write("# copies\nacc:ACC_0000300\n")
        rows, _ = self.scan_run(BACKGROUND + [INS, COPY],
                                "--exclude-acc-is6110", p)
        self.assertNotIn(COPY[0], rows)

    def test_vcf_without_the_flag_header_is_refused(self):
        with self.assertRaises(SystemExit):
            self.scan_run(BACKGROUND + [INS, COPY], "--exclude-acc-is6110",
                          self.vcf(header=False))

    def test_id_not_in_the_event_matrix_is_refused(self):
        with self.assertRaises(SystemExit):
            self.scan_run(BACKGROUND + [INS, COPY], "--exclude-acc-is6110",
                          self.vcf(flag=("acc:ACC_0000300", "acc:ACC_9")))

    def test_id_of_another_class_is_refused(self):
        p = os.path.join(self.d, "ids.txt")
        open(p, "w").write("is:300\n")
        with self.assertRaises(SystemExit):
            self.scan_run(BACKGROUND + [INS, COPY], "--exclude-acc-is6110", p)


class ScanCompare(unittest.TestCase):
    def setUp(self):
        self.cmp = ta.load("scan_compare", "assoc/bin/scan_compare.py")

    @staticmethod
    def row(k, q, rep="1", cls="small"):
        return {"row_key": k, "id": k, "cls": cls, "region": "core",
                "q_branch": q, "q_region": q, "q_lineage": q, "q_cond": "nan",
                "unit_rep": rep, "unit_id": f"u:{k}", "unit_members": k}

    def test_lists_only_the_differences(self):
        full = {"a": self.row("a", "0.01"), "b": self.row("b", "0.01"),
                "c": self.row("c", "0.2"),
                "x": self.row("x", "0.01", cls="accessory_presence"),
                "y": self.row("y", "0.5", cls="accessory_presence")}
        sens = {"a": self.row("a", "0.01"), "b": self.row("b", "0.06"),
                "c": self.row("c", "0.04")}
        listed, n = self.cmp.compare(full, sens)
        got = {o["row_key"]: o["change"] for o in listed}
        self.assertEqual(got, {
            "b": "significant_in_full_only",
            "c": "significant_in_sens_only",
            "x": "significant_in_full_only_excluded_from_sens"})
        self.assertEqual((n["surv_units_full"], n["surv_units_sens"],
                          n["rows_only_full"]), (3, 2, 2))

    def test_refuses_rows_only_in_the_sensitivity_scan(self):
        d = tempfile.mkdtemp()
        hdr = "row_key\tid\tcls\tregion\tq_branch\tq_region\tq_lineage\tq_cond\n"
        f = os.path.join(d, "f.tsv")
        s = os.path.join(d, "s.tsv")
        open(f, "w").write(hdr + "a\ta\tsmall\tcore\t1\t1\t1\tnan\n")
        open(s, "w").write(hdr + "b\tb\tsmall\tcore\t1\t1\t1\tnan\n")
        r = subprocess.run([sys.executable, "assoc/bin/scan_compare.py",
                            "--full", f, "--sens", s, "--out",
                            os.path.join(d, "c.tsv")], cwd=ROOT,
                           capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not of one event matrix", r.stderr)


class TailOptions(unittest.TestCase):
    """cohort_assoc_tail.sh: both options default on, are passed through,
    are recorded in assoc_options.tsv, and refuse a value that is not 0/1."""

    CHAIN = os.path.join(ROOT, "assoc", "bin", "cohort_assoc_tail.sh")

    def setup_run(self, vcf_extra_header="", vcf_records=""):
        d = tempfile.mkdtemp()
        os.makedirs(os.path.join(d, "assoc"))
        os.symlink(os.path.join(ROOT, "assoc", "bin"),
                   os.path.join(d, "assoc", "bin"))
        os.symlink(os.path.join(ROOT, "bin"), os.path.join(d, "bin"))
        for f in ("refbias/t/p5/sv_intervals.tsv",
                  "refbias/t/p5/svgt_iv_states.tsv", "refbias/t/p1/refmap.tsv",
                  "refbias/t.phenotype.tsv", "data/trees/t.snps.fasta",
                  "data/trees/t.combined.fasta",
                  "data/trees/t.combined.rooted.nwk", "data/trees/t.rooted.nwk",
                  "assoc/t/events/summary.txt", "p.txt"):
            os.makedirs(os.path.dirname(os.path.join(d, f)), exist_ok=True)
            with open(os.path.join(d, f), "w") as fh:
                fh.write("x\n")
        with gzip.open(os.path.join(d, "refbias/t/p5/merged.vcf.gz"), "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n##MTB_graph_build=b0\n"
                     + vcf_extra_header
                     + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                     + vcf_records)
        b = os.path.join(d, "refbias", "build", "b0")
        os.makedirs(os.path.join(b, "assets"))
        for f, txt in (("build_info.tsv", "build_id\tb0\n"),
                       ("assets/panel_polarity.tsv", "chrom\n"),
                       ("assets/graph_collapsed.vcf.gz", "x\n")):
            with open(os.path.join(b, f), "w") as fh:
                fh.write(txt)
        for f in ("data/trees/t.snps.fasta", "data/trees/t.combined.fasta",
                  "data/trees/t.combined.rooted.nwk", "data/trees/t.rooted.nwk",
                  "assoc/t/events/summary.txt"):
            r = subprocess.run(["bash", self.CHAIN, "t", "p.txt"], cwd=d,
                               capture_output=True, text=True,
                               env=dict(os.environ, MTB_CHAIN_PRINT_PROV=f))
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            with open(os.path.join(d, f + ".prov"), "w") as fh:
                fh.write("".join(l + "\n" for l in r.stdout.splitlines()
                                 if "\t" in l))
        log = os.path.join(d, "calls.log")
        stub = os.path.join(d, "stub.sh")
        with open(stub, "w") as fh:
            fh.write(f'#!/bin/bash\necho "$@" >> {log}\n')
        os.chmod(stub, 0o755)
        return d, log, dict(os.environ, MTB_PY=stub, MTB_PY_VT=stub)

    def run_chain(self, d, env, **kv):
        return subprocess.run(["bash", self.CHAIN, "t", "p.txt"], cwd=d,
                              capture_output=True, text=True,
                              env=dict(env, **kv))

    FLAGGED = ('##INFO=<ID=ACCIS6110,Number=0,Type=Flag,Description="x">\n',
               "NC_000962.3\t5\tacc:ACC_5\tA\t<INS>\t.\tPASS\t"
               "CLASS=accessory_presence;ACCIS6110\n")

    def test_defaults_run_both(self):
        d, log, env = self.setup_run(*self.FLAGGED)
        r = self.run_chain(d, env)
        self.assertEqual(r.returncode, 0, r.stderr)
        calls = read_text(log).splitlines()
        scans = [c for c in calls if "assoc_scan.py" in c]
        self.assertEqual(len(scans), 2, calls)
        self.assertNotIn("--no-collapse-identical", scans[0])
        self.assertIn("--out assoc/t/scan.tsv", scans[0])
        self.assertIn("--exclude-acc-is6110 refbias/t/p5/merged.vcf.gz",
                      scans[1])
        self.assertIn("--out assoc/t/scan.no_accis6110.tsv", scans[1])
        cmp_ = [c for c in calls if "scan_compare.py" in c]
        self.assertTrue(cmp_ and "--out assoc/t/scan.is6110_sensitivity.tsv"
                        in cmp_[0], calls)
        kv = read_kv(os.path.join(d, "assoc/t/assoc_options.tsv"))
        self.assertEqual((kv["collapse_identical"], kv["is6110_sensitivity"],
                          kv["accis6110_records"]), ("on", "on", "1"))

    def test_both_off(self):
        d, log, env = self.setup_run(*self.FLAGGED)
        stale = os.path.join(d, "assoc/t/scan.no_accis6110.tsv")
        open(stale, "w").write("old\n")
        r = self.run_chain(d, env, ASSOC_COLLAPSE_IDENTICAL="0",
                           ASSOC_IS6110_SENSITIVITY="0")
        self.assertEqual(r.returncode, 0, r.stderr)
        scans = [c for c in read_text(log).splitlines()
                 if "assoc_scan.py" in c]
        self.assertEqual(len(scans), 1)
        self.assertIn("--no-collapse-identical", scans[0])
        self.assertFalse(os.path.exists(stale))       # not read as this run's
        kv = read_kv(os.path.join(d, "assoc/t/assoc_options.tsv"))
        self.assertTrue(kv["collapse_identical"].startswith("off"))
        self.assertTrue(kv["is6110_sensitivity"].startswith("off"))

    def test_vcf_before_the_flag_skips_with_the_reason(self):
        d, log, env = self.setup_run()
        r = self.run_chain(d, env)
        self.assertEqual(r.returncode, 0, r.stderr)
        scans = [c for c in read_text(log).splitlines()
                 if "assoc_scan.py" in c]
        self.assertEqual(len(scans), 1)
        kv = read_kv(os.path.join(d, "assoc/t/assoc_options.tsv"))
        self.assertIn("no ACCIS6110 header", kv["is6110_sensitivity"])

    def test_bad_value_is_refused(self):
        d, log, env = self.setup_run(*self.FLAGGED)
        r = self.run_chain(d, env, ASSOC_IS6110_SENSITIVITY="yes")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("must be 0 or 1", r.stderr)


if __name__ == "__main__":
    unittest.main()
