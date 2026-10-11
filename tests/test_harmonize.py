#!/usr/bin/env python3
"""Tests for bin/harmonize_records.py (one event, one record; approved
2026-10-11). Same pattern as tests/run_tests.py: unittest, synthetic inputs.

    $MTB_PY tests/test_harmonize.py -v
"""
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))
import harmonize_records as h  # noqa: E402

C = "NC_000962.3"
GT = {"ALT": "1", "REF": "0", "ABSENT": "2", "NOCALL": "."}


def rec(pos, ident, ref, alts, info, states):
    return (C, pos, ident, ref, alts, info, "\t".join(f"{GT[s]}:{s}" for s in states))


def small(pos, ref, alt, states):
    return rec(pos, f"h37rv:{pos}:{ref}>{alt}", ref, [alt],
               ["CLASS=small", "KIND=indel", "REGION=core", "FRAME=h37rv"], states)


def sv(pos, typ, L, states, component="called", end=None):
    info = ["CLASS=sv", f"SVTYPE={typ}", f"SVLEN={'-' if typ == 'DEL' else ''}{L}",
            f"COMPONENT={component}", "QUALBAND=high"]
    if end:
        info.append(f"END={end}")
    return rec(pos, f"sv:{typ}:{pos}:{L}", "N", [f"<{typ}>"], info, states)


def isrec(pos, states):
    return rec(pos, f"h37rv:{pos}", "N", ["<INS:ME:IS6110>"],
               ["CLASS=is6110", "SVTYPE=INS", "SVLEN=1355"], states)


def by_id(out):
    return {r[2]: r for r in out}


def info(r):
    return dict(x.partition("=")[::2] for x in r[5])


class Rule1HaplotypeCallerOwnsSmall(unittest.TestCase):
    def test_sv_cell_handed_to_hc(self):
        recs = [small(1000, "A", "A" + "C" * 30, ["ALT", "REF", "REF"]),
                sv(1010, "INS", 35, ["ALT", "ALT", "NOCALL"])]
        out, n = h.harmonize(recs, C)
        r = by_id(out)["sv:INS:1010:35"]
        self.assertEqual(r[6].split("\t"), [h.DUP_CELL, "1:ALT", ".:NOCALL"])
        self.assertIn("sv_caller", info(by_id(out)["h37rv:1000:A>A" + "C" * 30])["SOURCES"])

    def test_large_sv_next_to_small_indel_kept(self):
        recs = [small(1000, "A", "A" + "C" * 20, ["ALT"]), sv(1010, "DEL", 400, ["ALT"])]
        out, _ = h.harmonize(recs, C)
        self.assertEqual(by_id(out)["sv:DEL:1010:400"][6], "1:ALT")

    def test_record_with_no_alt_left_dropped(self):
        recs = [small(1000, "A", "A" + "C" * 30, ["ALT"]), sv(1010, "INS", 35, ["ALT"])]
        out, n = h.harmonize(recs, C)
        self.assertNotIn("sv:INS:1010:35", by_id(out))


class Rule2SvBlockOwnsLarge(unittest.TestCase):
    def test_small_allele_merged_into_sv(self):
        recs = [small(2000, "G" + "T" * 80, "G", ["ALT", "REF", "NOCALL"]),
                sv(2000, "DEL", 80, ["ALT", "NOCALL", "NOCALL"], component="interval", end=2080)]
        out, n = h.harmonize(recs, C)
        ids = by_id(out)
        self.assertNotIn(recs[0][2], ids)
        r = ids["sv:DEL:2000:80"]
        self.assertEqual(info(r)["SMALLALLELE"], recs[0][2])
        self.assertEqual(r[6].split("\t"), ["1:ALT", "0:REF", ".:NOCALL"])
        self.assertEqual(info(r)["SOURCES"], "small,sv_catalogue")

    def test_dup_cell_not_refilled(self):
        # rule 3 hands sample 0's SV call to the IS6110 arm; rule 2's fill
        # from a merged small allele must not turn that DUP into REF
        big = "A" + "C" * 1382
        recs = [isrec(9000, ["ALT", "REF"]),
                small(9020, "A", big, ["REF", "ALT"]),
                sv(9020, "INS", 1382, ["ALT", "ALT"])]
        out, _ = h.harmonize(recs, C)
        self.assertEqual(by_id(out)["sv:INS:9020:1382"][6].split("\t")[0], h.DUP_CELL)

    def test_unmatched_large_small_allele_relabelled(self):
        recs = [small(3000, "G", "G" + "A" * 120, ["ALT"])]
        out, _ = h.harmonize(recs, C)
        i = info(out[0])
        self.assertEqual((i["CLASS"], i["SVTYPE"], i["SVLEN"], i["SVSOURCE"]),
                         ("sv", "INS", "120", "small_allele"))

    def test_different_carriers_not_merged(self):
        recs = [small(2000, "G" + "T" * 80, "G", ["ALT", "REF", "REF"]),
                sv(2000, "DEL", 80, ["REF", "ALT", "ALT"], component="interval", end=2080)]
        out, _ = h.harmonize(recs, C)
        self.assertEqual(len(out), 2)


class Rule3IS6110OwnsInsertions(unittest.TestCase):
    def test_truncated_ins_cell_handed_to_is6110(self):
        recs = [isrec(5000, ["ALT", "REF"]), sv(5040, "INS", 234, ["ALT", "ALT"])]
        out, n = h.harmonize(recs, C)
        r = by_id(out)["sv:INS:5040:234"]
        self.assertEqual(r[6].split("\t"), [h.DUP_CELL, "1:ALT"])
        self.assertEqual(n["rule3 SV INS ALT where the IS6110 arm says REF/ABSENT (kept)"], 1)
        self.assertIn("sv_caller", info(by_id(out)["h37rv:5000"])["SOURCES"])

    def test_relabelled_is_sized_allele_handed_over(self):
        recs = [isrec(6000, ["ALT"]), small(6010, "A", "A" + "C" * 1355, ["ALT"])]
        out, _ = h.harmonize(recs, C)
        self.assertNotIn(recs[1][2], by_id(out))     # its only ALT cell went to IS6110

    def test_relabelled_short_allele_near_is6110_kept(self):
        recs = [isrec(6000, ["ALT"]), small(6010, "A", "A" + "C" * 73, ["ALT"])]
        out, _ = h.harmonize(recs, C)
        self.assertEqual(by_id(out)[recs[1][2]][6], "1:ALT")

    def test_far_ins_untouched(self):
        recs = [isrec(5000, ["ALT"]), sv(5200, "INS", 234, ["ALT"])]
        out, _ = h.harmonize(recs, C)
        self.assertEqual(by_id(out)["sv:INS:5200:234"][6], "1:ALT")


class Rule4CopyLoss(unittest.TestCase):
    def test_copy_loss_labelled(self):
        recs = [sv(889019, "DEL", 1355, ["ALT", "REF"], component="interval", end=890380)]
        out, n = h.harmonize(recs, C, copies=[(889021, 890375)])
        i = info(out[0])
        self.assertEqual(i["CLASS"], "is6110")
        self.assertIn("IS6110LOSS", i)
        self.assertIn("LOWCONF", i)

    def test_relabelled_small_allele_copy_loss(self):
        ref = "T" + "G" * 1355
        recs = [small(889020, ref, "T", ["ALT", "REF"])]
        out, _ = h.harmonize(recs, C, copies=[(889021, 890375)])
        i = info(out[0])
        self.assertEqual((i["CLASS"], "IS6110LOSS" in i), ("is6110", True))

    def test_other_deletion_not_labelled(self):
        recs = [sv(880000, "DEL", 1355, ["ALT"], component="interval", end=881355)]
        out, _ = h.harmonize(recs, C, copies=[(889021, 890375)])
        self.assertEqual(info(out[0])["CLASS"], "sv")


class Sources(unittest.TestCase):
    def test_every_record_has_sources(self):
        recs = [small(10, "A", "G", ["ALT"]), isrec(100000, ["ALT"])]
        out, _ = h.harmonize(recs, C)
        self.assertTrue(all("SOURCES" in info(r) for r in out))

    def test_node_records_untouched(self):
        r = ("node_5", 1, "node:5:0:A>G", "A", ["G"], ["CLASS=small"], "1:ALT")
        out, _ = h.harmonize([r], C)
        self.assertEqual(out[0][6], "1:ALT")


if __name__ == "__main__":
    unittest.main()
