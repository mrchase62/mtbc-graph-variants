#!/usr/bin/env python3
"""Tests for two fixes from HANDOFF's fix list (2026-10-14):

AA-1  flag unreliable ancestral alleles. bin/panel_ancestral.py writes each
      key's Fitch change count over the panel's graph genomes (kit genomes
      left out) and a flag list: TIED / NODATA as before, OUTGROUP where only
      the outgroups decided, HOMOPLASTIC where changes >= --homoplastic-min
      (default 2, confirmed by the user 2026-10-14). bin/merge_cohort_vcf.py
      carries AA_CHANGES / AA_FLAG / AA_STATE and writes AA_SOURCE
      (reconstructed or assumed) on every record.

SPLIT-1 merge split SV records into one event, in the event writer
      (assoc/bin/write_event_matrix.py): the two P200 loci collapse to one
      event each; a VNTR copy-number series stays separate; the tatC-helY
      element keeps its polarity.

    $MTB_PY tests/test_audit_split_aa.py -v
"""
import csv
import hashlib
import importlib.util
import os
import random
import subprocess
import sys
import tempfile
import unittest
import warnings

import numpy as np

warnings.simplefilter("ignore", ResourceWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
sys.path.insert(0, os.path.join(ROOT, "bin"))
sys.path.insert(0, os.path.join(ROOT, "tests"))

import test_audit_parsimony as TP  # noqa: E402  (fixtures: VCF, tree, run)


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, rel))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


WEM = load("write_event_matrix", "assoc/bin/write_event_matrix.py")
MCV = load("merge_cohort_vcf", "bin/merge_cohort_vcf.py")
PA = load("panel_ancestral", "bin/panel_ancestral.py")


def rows_of(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def tsv(path, header, recs):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in recs:
            fh.write("\t".join(str(r.get(c, "")) for c in header) + "\n")
    return path


# ==========================================================================
# AA-1 in bin/panel_ancestral.py
# ==========================================================================
class ChangeCount(unittest.TestCase):
    """fitch_changes: the count over the graph genomes, kit leaves pruned."""

    def test_count_and_kit_skipped(self):
        children, root, leaves = PA.parse_newick("(K,((A,B),(C,D)));")
        R, X, N = PA.REF_BIT, PA.ALT_BIT, PA.ANY
        bits = {"A": np.array([X, X], np.uint8), "B": np.array([R, X], np.uint8),
                "C": np.array([X, R], np.uint8), "D": np.array([R, R], np.uint8),
                "K": np.array([X, X], np.uint8)}
        # key 0: A and C carry, B and D not -> 2 changes; key 1: (A,B) clade -> 1
        self.assertEqual(PA.fitch_changes(children, root, bits, 2, skip={"K"}).tolist(),
                         [2, 1])
        # K's ALT at key 1 would cost nothing either way (root takes ALT);
        # at key 0 a counted K adds nothing beyond the ingroup's 2
        bits["K"] = np.array([R, R], np.uint8)
        self.assertEqual(PA.fitch_changes(children, root, bits, 2, skip={"K"}).tolist(),
                         [2, 1])
        self.assertEqual(PA.fitch_changes(children, root, bits, 2).tolist(), [2, 1])
        # a missing leaf adds nothing
        bits["A"] = np.array([N, N], np.uint8)
        self.assertEqual(PA.fitch_changes(children, root, bits, 2, skip={"K"}).tolist(),
                         [1, 1])

    def test_flags(self):
        self.assertEqual(PA.aa_flags("ingroup", 1, 2), "")
        self.assertEqual(PA.aa_flags("ingroup", 2, 2), "HOMOPLASTIC")
        self.assertEqual(PA.aa_flags("outgroups", 1, 2), "OUTGROUP")
        self.assertEqual(PA.aa_flags("outgroups", 5, 2), "OUTGROUP,HOMOPLASTIC")
        self.assertEqual(PA.aa_flags("tied", 3, 2, "AG"), "TIED:AG,HOMOPLASTIC")
        self.assertEqual(PA.aa_flags("tied", 0, 2), "TIED")
        self.assertEqual(PA.aa_flags("nodata", 0, 2), "NODATA")


class TablesCarryFlags(TP._Base):
    """End to end on the P200 tree shape: kit genomes do not add changes,
    the threshold is a parameter, OUTGROUP marks method=outgroups."""

    TREE = TP.TieOnlyOutgroupsSolve.TREE
    OG = TP.TieOnlyOutgroupsSolve.OG

    def run_it(self, extra=()):
        seq = list("A" * 2000)
        for p, b in ((100, "C"), (200, "G"), (300, "T")):
            seq[p - 1] = b
        seq = "".join(seq)
        ref = TP.h37rv(f"{self.d}/h.fasta", seq)
        # samples A1 A2 B1 B2 L5.1 L6.1
        body = ("100 C T 1 0 1 0 0 0\n"    # homoplastic inside: 2 changes
                "200 G A 1 1 0 0 1 1\n"    # ingroup ties, outgroups say ALT: 1 change
                "300 T C 1 1 0 0 0 0\n")   # one change; the kit carries C (+0)
        vcf = TP.panel_vcf(f"{self.d}/p.vcf", ["A1", "A2", "B1", "B2", TP.L51, TP.L61],
                           body, 2000)
        kit = [(TP.L8, TP.states(f"{self.d}/l8.npz", seq, {300: "C"})),
               (TP.ET, TP.states(f"{self.d}/et.npz", seq, {300: "C"}))]
        r, anc, pol = TP.ancestral(self.d, self.TREE, vcf, ref, self.OG, kit, extra)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        return TP.by_pos(anc), TP.by_pos(pol)

    def test_default_threshold(self):
        a, p = self.run_it()
        self.assertEqual((a[(100, "C", "T")]["changes"], a[(100, "C", "T")]["flag"]),
                         ("2", "HOMOPLASTIC"))
        self.assertEqual((a[(200, "G", "A")]["AA"], a[(200, "G", "A")]["flag"],
                          a[(200, "G", "A")]["method"]), ("A", "OUTGROUP", "outgroups"))
        self.assertEqual(p[(200, "G", "A")]["flag"], "OUTGROUP")
        # the kit genomes carry C at 300: counted, they would make it 2
        # changes (HOMOPLASTIC); left out, it is 1. They still decide the
        # ingroup's tie (OUTGROUP)
        self.assertEqual((a[(300, "T", "C")]["changes"], a[(300, "T", "C")]["flag"]),
                         ("1", "OUTGROUP"))
        self.assertEqual(p[(100, "C", "T")]["changes"], "2")

    def test_threshold_parameter(self):
        a, _ = self.run_it(["--homoplastic-min", "3"])
        self.assertEqual(a[(100, "C", "T")]["flag"], "")

    def test_bad_threshold_fatal(self):
        r = TP.run([sys.executable, "bin/panel_ancestral.py", "--tree", "t",
                    "--panel-vcf", "v", "--h37rv", "h", "--ref-sample", "x",
                    "--outgroups", "o", "--out-ancestral", "a", "--out-polarity", "p",
                    "--homoplastic-min", "0"])
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("homoplastic-min", r.stderr)


# ==========================================================================
# AA-1 in the merged VCF
# ==========================================================================
class AaSourceInVcf(unittest.TestCase):

    def test_annotate(self):
        pol = {(300, "A", "AT"): ("REF", "3", "HOMOPLASTIC"),
               (500, "C", "C" + "G" * 60): ("ALT", "1", "OUTGROUP"),
               (700, "T", "T" + "A" * 70): ("REF", "0", ""),
               (701, "T", "T" + "C" * 70): ("ALT", "0", "")}
        recs = [
            ("NC_000962.3", 100, "h37rv:100:A>G", "A", ["G"], ["CLASS=small", "AA=G"], ""),
            ("NC_000962.3", 200, "h37rv:200:C>T", "C", ["T"],
             ["CLASS=small", "AA=.", "AA_FLAG=TIED:CT"], ""),
            ("NC_000962.3", 300, "h37rv:300:A>AT", "A", ["AT", "*"], ["CLASS=small"], ""),
            ("NC_000962.3", 500, "sv:INS:500:60", "C", ["<INS>"],
             ["CLASS=sv", "SMALLALLELE=h37rv:500:C>C" + "G" * 60], ""),
            ("NC_000962.3", 700, "sv:INS:700:70", "T", ["<INS>"],
             ["CLASS=sv", "SMALLALLELE=h37rv:700:T>T" + "A" * 70
              + ",h37rv:701:T>T" + "C" * 70], ""),
            ("NC_000962.3", 900, "sv:INS:900:80", "G", ["<INS>"], ["CLASS=sv"], ""),
            ("node_5", 1, "node:5:3:A>C", "A", ["C"], ["CLASS=small"], ""),
        ]
        out, n = MCV.annotate_aa_source(recs, "NC_000962.3", pol)
        info = {r[2]: r[5] for r in out}
        self.assertIn("AA_SOURCE=reconstructed", info["h37rv:100:A>G"])
        self.assertIn("AA_SOURCE=assumed", info["h37rv:200:C>T"])
        self.assertEqual([x for x in info["h37rv:300:A>AT"] if x.startswith("AA")],
                         ["AA_STATE=REF", "AA_CHANGES=3", "AA_FLAG=HOMOPLASTIC",
                          "AA_SOURCE=reconstructed"])
        self.assertIn("AA_STATE=ALT", info["sv:INS:500:60"])
        self.assertIn("AA_FLAG=OUTGROUP", info["sv:INS:500:60"])
        # two small alleles that disagree: no claim
        self.assertIn("AA_SOURCE=assumed", info["sv:INS:700:70"])
        self.assertIn("AA_SOURCE=assumed", info["sv:INS:900:80"])
        self.assertIn("AA_SOURCE=assumed", info["node:5:3:A>C"])
        self.assertEqual(n, {"reconstructed": 3, "assumed": 4})
        # every record has exactly one AA_SOURCE
        for r in out:
            self.assertEqual(sum(x.startswith("AA_SOURCE=") for x in r[5]), 1)

    def test_merge_end_to_end(self):
        """The AA join keeps AA_INVERTED on a flagged site (it used to be
        written only when AA_FLAG was empty), and the header declares the
        new tags."""
        with tempfile.TemporaryDirectory() as d:
            samples = ["S1", "S2", "S3"]
            keys = [dict(key="h37rv:100:A>G", frame="h37rv", h37rv_pos=100, node="",
                         node_offset="", canonical_ref="A", canonical_alt="G",
                         region="core", kind="SNP", acc_locus=""),
                    dict(key="h37rv:300:A>AT", frame="h37rv", h37rv_pos=300, node="",
                         node_offset="", canonical_ref="A", canonical_alt="AT",
                         region="core", kind="INDEL", acc_locus="")]
            tsv(f"{d}/keys.tsv", list(keys[0]), keys)
            sha = hashlib.sha1("\n".join(k["key"] for k in keys).encode()).hexdigest()
            for sm in samples:
                with open(f"{d}/{sm}.states.tsv", "w") as fh:
                    fh.write(f"#format\tsparse-v1\n#sample\t{sm}\n#default\tREF\n"
                             f"#n_keys\t2\n#keys_sha1\t{sha}\nidx\tstate\tallele\n"
                             + ("0\tALT\tG\n1\tALT\tAT\n" if sm == "S1" else ""))
            tsv(f"{d}/refmap.tsv", ["sample", "reference"],
                [dict(sample=s, reference="R") for s in samples])
            tsv(f"{d}/anc.tsv", ["chrom", "pos", "ref", "alt", "AA", "flag", "polarity",
                                 "changes", "method"],
                [dict(chrom="NC_000962.3", pos=100, ref="A", alt="G", AA="G",
                      flag="OUTGROUP,HOMOPLASTIC", polarity="alt_ancestral",
                      changes=4, method="outgroups")])
            tsv(f"{d}/pol.tsv", ["chrom", "pos", "ref", "alt", "outgroup_gt", "ancestral",
                                 "panel_af", "n_records", "method", "changes", "flag"],
                [dict(chrom="NC_000962.3", pos=300, ref="A", alt="AT", outgroup_gt=0,
                      ancestral="REF", panel_af="0.1", n_records=1, method="ingroup",
                      changes=1, flag="")])
            r = subprocess.run(
                [sys.executable, "bin/merge_cohort_vcf.py", "--states-dir", d,
                 "--refmap", f"{d}/refmap.tsv", "--keys", f"{d}/keys.tsv",
                 "--cohort-name", "t", "--build-id", "b", "--ancestral", f"{d}/anc.tsv",
                 "--panel-polarity", f"{d}/pol.tsv", "--out", f"{d}/o.vcf"],
                capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
            lines = open(f"{d}/o.vcf").read().splitlines()
            r2 = subprocess.run(
                [sys.executable, "bin/merge_cohort_vcf.py", "--states-dir", d,
                 "--refmap", f"{d}/refmap.tsv", "--keys", f"{d}/keys.tsv",
                 "--cohort-name", "t", "--build-id", "b", "--ancestral", f"{d}/anc.tsv",
                 "--panel-polarity", f"{d}/missing.tsv", "--out", f"{d}/o2.vcf"],
                capture_output=True, text=True)
        self.assertNotEqual(r2.returncode, 0)
        self.assertIn("FATAL", r2.stderr + r2.stdout)
        hdr = "\n".join(l for l in lines if l.startswith("##INFO"))
        for tag in ("AA_FLAG,Number=.", "AA_CHANGES", "AA_SOURCE", "AA_STATE"):
            self.assertIn(tag, hdr)
        info = {l.split("\t")[2]: l.split("\t")[7].split(";")
                for l in lines if not l.startswith("#")}
        snp = info["h37rv:100:A>G"]
        for x in ("AA=G", "AA_FLAG=OUTGROUP,HOMOPLASTIC", "AA_CHANGES=4",
                  "AA_INVERTED", "AA_SOURCE=reconstructed"):
            self.assertIn(x, snp)
        ind = info["h37rv:300:A>AT"]
        for x in ("AA_STATE=REF", "AA_CHANGES=1", "AA_SOURCE=reconstructed"):
            self.assertIn(x, ind)
        self.assertFalse(any(x.startswith("AA=") for x in ind))


# ==========================================================================
# SPLIT-1 in the event writer
# ==========================================================================
def rseq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(s, k, rng):
    s = list(s)
    for _ in range(k):
        i = rng.randrange(len(s))
        s[i] = "ACGT"[("ACGT".index(s[i]) + 1) % 4]
    return "".join(s)


def ins(pos, ins_seq, svlen=None, rid=None, cls="sv", anchor="C"):
    """A variant dict as read_vcf builds it, for a sequence-resolved INS."""
    alt = anchor + ins_seq
    return dict(key=(pos, anchor, alt), chrom="NC_000962.3", pos=pos,
                id=rid or f"h37rv:{pos}:{anchor}>{alt[:12]}", ref=anchor, alt=alt,
                cls=cls, frame="h37rv", aa="", panel_af="", aa_flag="",
                polarity="presence_is_derived", derived="ALT", region="core",
                evidence="", siteclass="", aa_changes="", aa_source="",
                svtype="INS", svlen=str(svlen if svlen is not None else len(ins_seq)),
                smallallele="", acclocus="")


def gts(n_samples, carriers):
    g = np.full(n_samples, ord("0"), dtype=np.uint8)
    for c in carriers:
        g[c] = ord("1")
    return g


class SplitSvMerge(unittest.TestCase):

    def setUp(self):
        self.rng = random.Random(20261014)

    def merged(self, vs, carriers, pol=None, n=200):
        rows = [gts(n, c) for c in carriers]
        return WEM.merge_split_svs(list(vs), rows, pol or {}, 0.05)

    def test_2268721_one_event(self):
        # P200: ~5 kb insertion at 2,268,721 as 57 records; the main family
        # (4,999-5,025 bp, one-base differences, starts spread over 50 bp)
        # is one event. Length outliers (6,358 = the 5 kb plus an IS6110
        # copy; 752) are other events by the 5% rule
        base = rseq(5000, self.rng)
        lens = [5000] * 40 + [5001, 5002, 4999, 5025]
        vs = []
        for k, L in enumerate(lens):
            s = mutate(base, 1 + k % 3, self.rng)[:L] if L <= 5000 else \
                mutate(base, 1, self.rng) + rseq(L - 5000, self.rng)
            vs.append(ins(2268721 + (k % 3), s, rid=f"m{k}"))
        vs.append(ins(2268721, base + rseq(1358, self.rng), rid="plusIS"))
        vs.append(ins(2268721, rseq(752, self.rng), rid="short"))
        carriers = [[k] for k in range(len(vs))]
        out, rows, n = self.merged(vs, carriers)
        self.assertEqual(len(out), 3)
        ev = [v for v in out if v.get("sv_members")]
        self.assertEqual(len(ev), 1)
        self.assertEqual(ev[0]["sv_members"], len(lens))
        self.assertEqual(sorted(ev[0]["sv_member_ids"].split(",")),
                         sorted(f"m{k}" for k in range(len(lens))))
        i = out.index(ev[0])
        self.assertEqual(int((rows[i] == ord("1")).sum()), len(lens))  # any member
        self.assertEqual(n["SV records merged away"], len(lens) - 1)

    def test_2219417_one_event(self):
        base = rseq(4496, self.rng)
        lens = [4496] * 9 + [4495] * 8 + [4494] * 5
        vs = [ins(2219417, mutate(base, 1, self.rng)[:L], rid=f"r{k}")
              for k, L in enumerate(lens)]
        vs.append(ins(2219417, rseq(3417, self.rng), rid="other"))
        out, _, n = self.merged(vs, [[k] for k in range(len(vs))])
        self.assertEqual(len(out), 2)
        self.assertEqual(max(v.get("sv_members", 1) for v in out), 22)

    def test_vntr_series_separate(self):
        # 802,427: unit 54, one to six copies; each count is its own event,
        # and two records of the SAME count (phases 1 bp apart) merge
        unit = rseq(54, self.rng)
        vs = [ins(802427, unit * k, rid=f"c{k}") for k in range(1, 7)]
        vs.append(ins(802428, unit[1:] + unit[:1], rid="c1b"))
        out, _, _ = self.merged(vs, [[k] for k in range(len(vs))])
        self.assertEqual(len(out), 6)
        one = [v for v in out if v.get("sv_members")]
        self.assertEqual(len(one), 1)
        self.assertEqual(sorted(one[0]["sv_member_ids"].split(",")), ["c1", "c1b"])

    def test_long_vntr_counts_within_5_percent_stay_separate(self):
        # 21 and 22 copies of a 56 bp unit differ by 4.7%: the copy-number
        # guard, not the length rule, keeps them apart
        unit = rseq(56, self.rng)
        vs = [ins(38982, unit * 21, rid="c21"), ins(38982, unit * 22, rid="c22")]
        out, _, _ = self.merged(vs, [[0], [1]])
        self.assertEqual(len(out), 2)

    def test_duplication_copies_separate(self):
        # 1,480,930: the 1,674 bp Rv1318c-Rv1320c duplication. One copy
        # (1,671-1,701 bp: one event) and two copies (3,347-3,348) stay apart
        unit = rseq(1674, self.rng)
        one = [ins(1480930, mutate(unit, 1, self.rng)[:L], rid=f"one{L}")
               for L in (1674, 1673, 1671)] + \
            [ins(1480935, unit + rseq(27, self.rng), rid="one1701")]
        two = [ins(1480930, unit * 2, rid="two3348"),
               ins(1480940, (unit * 2)[:3347], rid="two3347")]
        out, _, _ = self.merged(one + two, [[k] for k in range(6)])
        self.assertEqual(len(out), 2)
        self.assertEqual(sorted(v["sv_members"] for v in out), [2, 4])

    def test_tatc_hely_element_unchanged(self):
        # the main 68 bp allele is ALT-ancestral (outgroups), its rare
        # one-base variants REF-ancestral: they are derived variants of the
        # element, not a disagreement. The merged event keeps the main
        # allele's polarity and key; carriers are the union
        el = TP.ELEMENT[1:]
        main = ins(2352065, el, rid="main")
        var = [ins(2352065, mutate(el, 1, self.rng), rid=f"v{k}") for k in range(3)]
        var.append(ins(2352065, el[:67], rid="v67"))
        vs = [main] + var
        pol = {(v["chrom"], v["pos"], v["ref"], v["alt"]): ("REF", "0.0050")
               for v in var}
        pol[("NC_000962.3", 2352065, "C", TP.ELEMENT)] = ("ALT", "0.7538")
        carriers = [list(range(0, 150)), [150], [151], [152], [153]]
        out, rows, n = self.merged(vs, carriers, pol)
        self.assertEqual(len(out), 1)
        ev = out[0]
        self.assertEqual(ev["id"], "main")
        self.assertEqual(ev["key"], main["key"])
        self.assertEqual((ev["polarity"], ev["derived"], ev["split_flag"]),
                         ("alt_ancestral_outgroup", "REF", ""))
        self.assertEqual(int((rows[0] == ord("1")).sum()), 154)
        # alone, the main record gets the same polarity
        out1, _, _ = self.merged([ins(2352065, el, rid="main")], [carriers[0]], pol)
        self.assertEqual(out1[0]["polarity"], "alt_ancestral_outgroup")

    def test_polarity_rules(self):
        s = rseq(300, self.rng)
        a, b = ins(5000, s, rid="a"), ins(5010, mutate(s, 1, self.rng), rid="b")
        k = lambda v: (v["chrom"], v["pos"], v["ref"], v["alt"])  # noqa: E731
        # two members ALT-ancestral: conflict
        out, _, n = self.merged([dict(a), dict(b)], [[0], [1]],
                                {k(a): ("ALT", "0.5"), k(b): ("ALT", "0.4")})
        self.assertEqual((out[0]["polarity"], out[0]["split_flag"]),
                         ("unpolarised", "SPLIT_CONFLICT"))
        # all REF
        out, _, _ = self.merged([dict(a), dict(b)], [[0], [1]],
                                {k(a): ("REF", "0.1"), k(b): ("REF", "0.1")})
        self.assertEqual((out[0]["polarity"], out[0]["split_flag"]),
                         ("ref_ancestral_outgroup", ""))
        # REF plus a member without a row: REF, flagged partial
        out, _, _ = self.merged([dict(a), dict(b)], [[0], [1]], {k(a): ("REF", "0.1")})
        self.assertEqual((out[0]["polarity"], out[0]["split_flag"]),
                         ("ref_ancestral_outgroup", "SPLIT_PARTIAL"))
        # a rare ALT-ancestral row is not trusted (the 5% guard): no row
        out, _, _ = self.merged([dict(a)], [[0]], {k(a): ("ALT", "0.01")})
        self.assertEqual(out[0]["polarity"], "presence_is_derived")
        # no rows at all: the record's own polarity
        out, _, _ = self.merged([dict(a), dict(b)], [[0], [1]])
        self.assertEqual(out[0]["polarity"], "presence_is_derived")

    def test_symbolic_member_by_smallallele(self):
        s = rseq(100, self.rng)
        v = ins(7000, s, rid="sv:INS:7000:100")
        v.update(alt="<INS>", ref="C", key=(7000, "C", "<INS:100>"),
                 smallallele=f"h37rv:7000:C>C{s}")
        out, _, _ = self.merged([v], [[0]], {("NC_000962.3", 7000, "C", "C" + s):
                                             ("REF", "0.2")})
        self.assertEqual(out[0]["polarity"], "ref_ancestral_outgroup")

    def test_types_and_classes_kept_apart(self):
        s = rseq(200, self.rng)
        i1 = ins(9000, s, rid="ins")
        d1 = ins(9000, s, rid="del")
        d1.update(svtype="DEL", svlen="-200", ref="C" + s, alt="C",
                  key=(9000, "C" + s, "C"))
        small = ins(9005, s, rid="small", cls="small")
        out, _, _ = self.merged([i1, d1, small], [[0], [1], [2]])
        self.assertEqual(len(out), 3)

    def test_merged_genotype(self):
        o = ord
        rows = [np.array([o("1"), o("0"), o("."), o("."), o("2")], np.uint8),
                np.array([o("0"), o("."), o("2"), o("."), o(".")], np.uint8),
                np.array([o("."), o("1"), o("."), o("."), o("0")], np.uint8)]
        self.assertEqual(bytes(WEM.merged_gt(rows)).decode(), "112.0")


class SplitSvEndToEnd(unittest.TestCase):
    """The writer run on a VCF: one variants.tsv row for the split event,
    its members in sv_merged.tsv, and --no-merge-split-sv restoring them."""

    def setUp(self):
        p = os.environ.get("MTB_PY_VT", "")
        ok = p and os.access(p, os.X_OK) and \
            TP.run([p, "-c", "import mtbvartools, dendropy"]).returncode == 0
        if not ok:
            self.skipTest("MTB_PY_VT with mtbvartools not set")
        self.py = p

    def test_rows(self):
        rng = random.Random(1)
        s = rseq(120, rng)
        s2 = mutate(s, 1, rng)
        with tempfile.TemporaryDirectory() as d:
            hdr = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t"
                   "INFO\tFORMAT\tS1\tS2\tS3\tS4\n")
            body = (f"NC_000962.3\t500\th37rv:500:C>C{s}\tC\tC{s}\t.\tPASS\t"
                    "CLASS=sv;SVTYPE=INS;SVLEN=120;FRAME=h37rv;AA_SOURCE=reconstructed;"
                    "AA_STATE=REF;AA_CHANGES=2\tGT\t1\t0\t0\t0\n"
                    f"NC_000962.3\t510\th37rv:510:C>C{s2}\tC\tC{s2}\t.\tPASS\t"
                    "CLASS=sv;SVTYPE=INS;SVLEN=120;FRAME=h37rv;AA_SOURCE=assumed\t"
                    "GT\t0\t0\t1\t0\n"
                    "NC_000962.3\t600\th37rv:600:G>GA\tG\tGA\t.\tPASS\t"
                    "CLASS=small;FRAME=h37rv;AA_SOURCE=assumed\tGT\t1\t1\t0\t0\n")
            TP.write(f"{d}/m.vcf", hdr + body)
            TP.write(f"{d}/t.nwk", "(O:1,((S1:1,S2:1)x:1,(S3:1,S4:1)y:1)ing:1);\n")
            TP.write(f"{d}/og.fa", ">O\nT\n")
            TP.write(f"{d}/og.tsv", "column\tchrom\tpos\tref\talt\n0\tNC_000962.3\t1\tA\tT\n")
            pol = TP.write(f"{d}/pol.tsv",
                           "chrom\tpos\tref\talt\toutgroup_gt\tancestral\tpanel_af\t"
                           "n_records\tmethod\tchanges\tflag\n"
                           f"{TP.PANSN}\t500\tC\tC{s}\t0\tREF\t0.2\t1\tingroup\t2\tHOMOPLASTIC\n")
            base = [self.py, "assoc/bin/write_event_matrix.py", "--vcf", f"{d}/m.vcf",
                    "--tree", f"{d}/t.nwk", "--outgroup-name", "O",
                    "--outgroup-fasta", f"{d}/og.fa", "--outgroup-sites", f"{d}/og.tsv",
                    "--panel-polarity", pol]
            r = TP.run(base + ["--out", f"{d}/ev"])
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            v = rows_of(f"{d}/ev/variants.tsv")
            m = rows_of(f"{d}/ev/sv_merged.tsv")
            summ = open(f"{d}/ev/summary.txt").read()
            r = TP.run(base + ["--out", f"{d}/ev0", "--no-merge-split-sv"])
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            v0 = rows_of(f"{d}/ev0/variants.tsv")
        sv = [x for x in v if x["class"] == "sv"]
        self.assertEqual(len(sv), 1)
        self.assertEqual((sv[0]["sv_members"], sv[0]["polarity"], sv[0]["split_flag"]),
                         ("2", "ref_ancestral_outgroup", "SPLIT_PARTIAL"))
        self.assertEqual(sv[0]["n_derived_leaves"], "2")      # S1 and S3
        self.assertEqual(sv[0]["aa_changes"], "")             # merged: see members
        self.assertEqual([x["member_id"] for x in m],
                         [f"h37rv:500:C>C{s}", f"h37rv:510:C>C{s2}"])
        self.assertEqual([x["member_state"] for x in m], ["REF", "."])
        small = [x for x in v if x["class"] == "small"][0]
        self.assertEqual(small["aa_source"], "assumed")
        self.assertIn("split SV records  1 merged into 1 events", summ)
        self.assertIn("AA_SOURCE in the VCF", summ)
        self.assertEqual(len([x for x in v0 if x["class"] == "sv"]), 2)
        self.assertTrue(all(x["polarity"] == "presence_is_derived"
                            for x in v0 if x["class"] == "sv"))


# ==========================================================================
# SPLIT-1 extension: truncated small-allele insertions -> the graph allele
# ==========================================================================
def call(pos, seq, rid, svsource="small_allele"):
    v = ins(pos, seq, rid=rid)
    v["svsource"] = svsource
    return v


def graph(pos, seq):
    return dict(pos=pos, ref="C", alt="C" + seq, seq=seq, len=len(seq))


def other(cls, pos, rid, length):
    """An accessory_presence or IS6110 record as read_vcf builds it."""
    v = ins(pos, "", rid=rid, cls=cls)
    v.update(alt="<INS>", ref="C", key=(pos, "C", f"<INS:{cls}>"), svlen=str(length),
             svtype="INS", acclen=str(length) if cls == "accessory_presence" else "",
             acclocus=f"ACC_{pos}" if cls == "accessory_presence" else "",
             polarity="presence_is_derived")
    return v


class SplitTruncated(unittest.TestCase):
    """HaplotypeCaller's assembled insertion stops short of a long insertion
    (ctrl25_p200: 130 of 136 right-site/wrong-length calls; 533 bp at
    742,634 called as 60 and 74 bp)."""

    def setUp(self):
        self.rng = random.Random(742634)

    def run_merge(self, vs, carriers, gal, pol=None, n=30):
        rows = [gts(n, c) for c in carriers]
        return WEM.merge_split_svs(list(vs), rows, pol or {}, 0.05, gal)

    def test_742634_calls_join_the_accessory_event(self):
        x = rseq(533, self.rng)
        gal = [graph(742634, x)] + [graph(742634, mutate(x, 1, self.rng)) for _ in range(3)]
        acc = other("accessory_presence", 742634, "acc:ACC_0742634", 533)
        vs = [call(742634, x[:60], "c60"), call(742634, x[:74], "c74"), acc]
        out, rows, n = self.run_merge(vs, [[0, 1], [1, 2], list(range(25))], gal)
        self.assertEqual(len(out), 1)
        ev = out[0]
        self.assertEqual((ev["id"], ev["cls"], ev["polarity"]),
                         ("acc:ACC_0742634", "accessory_presence", "presence_is_derived"))
        self.assertEqual(int((rows[0] == ord("1")).sum()), 25)
        why = {m["id"]: (m["reason"], m["graph_allele"]) for m in ev["_members"]}
        self.assertEqual(why["c60"], ("truncated_to_graph_allele", "742634:533"))
        self.assertEqual(why["c74"][0], "truncated_to_graph_allele")
        self.assertEqual(why["acc:ACC_0742634"][0], "same_event")

    def test_calls_alone_form_the_graph_allele_event(self):
        # no cohort record of the 533 bp allele: the two calls become one
        # event, polarised by the graph allele's own panel row
        x = rseq(533, self.rng)
        gal = [graph(742634, x)]
        vs = [call(742634, x[:60], "c60"), call(742638, x[10:84], "c74")]
        pol = {("NC_000962.3", 742634, "C", "C" + x): ("REF", "0.9")}
        out, rows, n = self.run_merge(vs, [[0], [1]], gal, pol)
        self.assertEqual(len(out), 1)
        self.assertEqual((out[0]["sv_members"], out[0]["polarity"]),
                         (2, "ref_ancestral_outgroup"))
        self.assertEqual(n["truncated calls formed a graph-allele event"], 2)

    def test_calls_join_a_full_length_sv_record(self):
        x = rseq(400, self.rng)
        gal = [graph(5000, x)]
        full = call(5000, x, "full")
        vs = [call(5000, x[:80], "short"), full]
        out, rows, _ = self.run_merge(vs, [[0], [1, 2]], gal)
        self.assertEqual(len(out), 1)
        self.assertEqual(out[0]["id"], "full")             # full-length representative
        self.assertEqual(int((rows[0] == ord("1")).sum()), 3)

    def test_is6110_piece_joins_the_site(self):
        x = rseq(1355, self.rng)
        gal = [graph(34866, x)]
        site = other("is6110", 34869, "h37rv:34869", 1355)
        vs = [call(34869, x[:64], "piece"), site]
        out, rows, _ = self.run_merge(vs, [[0], [0, 1]], gal)
        self.assertEqual(len(out), 1)
        self.assertEqual((out[0]["id"], out[0]["cls"]), ("h37rv:34869", "is6110"))

    def test_802427_vntr_unit_call_left_and_flagged(self):
        # a one-unit call: itself the 54 bp allele, and inside every longer
        # copy number -- cannot be placed
        u = rseq(54, self.rng)
        gal = [graph(802427, u * k) for k in (1, 2, 3, 4)]
        vs = [call(802427, u, "u1")]
        out, _, n = self.run_merge(vs, [[0]], gal)
        self.assertEqual(len(out), 1)
        self.assertEqual(out[0]["split_flag"], "TRUNC_POSSIBLE")
        self.assertEqual(out[0]["_members"][0]["reason"], "truncated_possible")
        self.assertEqual(out[0]["sv_members"], 1)          # merged with nothing
        self.assertEqual(n["truncated calls possible (TRUNC_POSSIBLE)"], 1)

    def test_vntr_single_longer_candidate_joins(self):
        # the graph has only the three-copy allele: one candidate, so the
        # one-unit call is its truncated reading
        u = rseq(54, self.rng)
        gal = [graph(802427, u * 3)]
        vs = [call(802427, u, "u1"), call(802427, u * 3, "u3")]
        out, _, _ = self.run_merge(vs, [[0], [1]], gal)
        self.assertEqual(len(out), 1)
        self.assertEqual(out[0]["id"], "u3")

    def test_two_longer_candidates_ambiguous(self):
        s = rseq(60, self.rng)
        gal = [graph(9000, s + rseq(200, self.rng)), graph(9010, rseq(300, self.rng) + s)]
        vs = [call(9000, s, "c")]
        out, _, n = self.run_merge(vs, [[0]], gal)
        self.assertEqual(len(out), 1)
        self.assertEqual(out[0]["split_flag"], "TRUNC_AMBIGUOUS")
        self.assertEqual(n["truncated calls ambiguous (TRUNC_AMBIGUOUS)"], 1)

    def test_not_contained_left_alone(self):
        x = rseq(500, self.rng)
        gal = [graph(9000, x)]
        vs = [call(9000, rseq(70, self.rng), "novel")]
        out, _, n = self.run_merge(vs, [[0]], gal)
        self.assertEqual((len(out), out[0]["split_flag"]), (1, ""))
        self.assertNotIn("_members", out[0])
        # half contained (below 90%): also left alone
        vs = [call(9000, x[:35] + rseq(35, self.rng), "half")]
        out, _, _ = self.run_merge(vs, [[0]], gal)
        self.assertNotIn("_members", out[0])

    def test_far_or_caller_records_not_moved(self):
        x = rseq(500, self.rng)
        gal = [graph(9000, x)]
        far = call(9100, x[:70], "far")                       # 100 bp away
        hc = call(9000, x[:70], "caller", svsource="")        # not small_allele
        out, _, _ = self.run_merge([far, hc], [[0], [1]], gal)
        self.assertEqual(len(out), 2)
        self.assertTrue(all("_members" not in v for v in out))

    def test_containment(self):
        x = rseq(300, self.rng)
        self.assertEqual(WEM.contained_frac(x[100:180], x), 1.0)
        self.assertLess(WEM.contained_frac(rseq(80, self.rng), x), 0.1)
        self.assertEqual(WEM.contained_frac("ACGT", x), 0.0)    # shorter than k


class SplitTruncatedEndToEnd(SplitSvEndToEnd):
    """--graph-vcf on the writer: the truncated call joins and sv_merged.tsv
    says why."""

    def test_rows(self):
        rng = random.Random(7)
        x = rseq(300, rng)
        with tempfile.TemporaryDirectory() as d:
            TP.write(f"{d}/g.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\t"
                     "QUAL\tFILTER\tINFO\tFORMAT\tG1\n"
                     f"{TP.PANSN}\t500\t.\tC\tC{x}\t.\t.\t.\tGT\t1\n")
            TP.write(f"{d}/m.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\t"
                     "QUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
                     f"NC_000962.3\t500\th37rv:500:C>C{x}\tC\tC{x}\t.\tPASS\t"
                     "CLASS=sv;SVTYPE=INS;SVLEN=300;SVSOURCE=small_allele;FRAME=h37rv\t"
                     "GT\t1\t0\t0\t0\n"
                     f"NC_000962.3\t500\th37rv:500:C>C{x[:70]}\tC\tC{x[:70]}\t.\tPASS\t"
                     "CLASS=sv;SVTYPE=INS;SVLEN=70;SVSOURCE=small_allele;FRAME=h37rv\t"
                     "GT\t0\t1\t0\t0\n")
            TP.write(f"{d}/t.nwk", "(O:1,((S1:1,S2:1)x:1,(S3:1,S4:1)y:1)ing:1);\n")
            TP.write(f"{d}/og.fa", ">O\nT\n")
            TP.write(f"{d}/og.tsv", "column\tchrom\tpos\tref\talt\n0\tNC_000962.3\t1\tA\tT\n")
            base = [self.py, "assoc/bin/write_event_matrix.py", "--vcf", f"{d}/m.vcf",
                    "--tree", f"{d}/t.nwk", "--outgroup-name", "O",
                    "--outgroup-fasta", f"{d}/og.fa", "--outgroup-sites", f"{d}/og.tsv",
                    "--panel-polarity", ""]
            r = TP.run(base + ["--out", f"{d}/ev", "--graph-vcf", f"{d}/g.vcf"])
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            v = rows_of(f"{d}/ev/variants.tsv")
            m = rows_of(f"{d}/ev/sv_merged.tsv")
            r0 = TP.run(base + ["--out", f"{d}/ev0"])
            self.assertEqual(r0.returncode, 0, r0.stdout + r0.stderr)
            v0 = rows_of(f"{d}/ev0/variants.tsv")
            bad = TP.run(base + ["--out", f"{d}/ev1", "--graph-vcf", f"{d}/missing.vcf"])
        self.assertEqual(len(v), 1)
        self.assertEqual(v[0]["n_derived_leaves"], "2")
        self.assertEqual({x["reason"] for x in m}, {"same_event", "truncated_to_graph_allele"})
        self.assertEqual([x["graph_allele"] for x in m if x["reason"] != "same_event"],
                         ["500:300"])
        self.assertEqual(len(v0), 2)                 # without --graph-vcf: apart
        self.assertNotEqual(bad.returncode, 0)


if __name__ == "__main__":
    unittest.main()
