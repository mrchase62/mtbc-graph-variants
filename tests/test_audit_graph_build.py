#!/usr/bin/env python3
"""Regression tests for the 2026-10-05 audit, group "graph_build".

PGB-1 (pggb_build.sh could not start a clean build), PGB-2 (its defaults were
not CX333's settings; effective settings unrecorded), PGB-9 (production readers
took the decomposed graph VCF; sync_back.sh dropped it from the mirror while
they did), PGB-12 (snp_nonredundant.py on duplicate-key input), and the
deconstruct header (vg_deconstruct.sh recorded no settings).

pggb and vg are never run: a stub `singularity` on PATH records the command
line it was given. Same pattern as tests/run_tests.py: small synthetic inputs,
standard library plus the pipeline's own scripts, seconds to run.

    $MTB_PY tests/test_audit_graph_build.py
"""
import csv
import gzip
import os
import re
import stat
import subprocess
import sys
import tempfile
import textwrap
import unittest
import warnings

warnings.simplefilter("ignore", ResourceWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

COLLAPSED = "graphs/CX333.s10k.k23.K15/all_variants.collapsed.vcf.gz"


def write(path, text, mode=None):
    with open(path, "w") as fh:
        fh.write(textwrap.dedent(text).lstrip("\n"))
    if mode:
        os.chmod(path, mode)
    return path


def py():
    return sys.executable


def stub_env(d, extra=None):
    """A fake MTB_WORK under d, and a stub singularity that logs its argv and
    answers --version / vg version / vg deconstruct."""
    os.makedirs(f"{d}/work/data/fastas", exist_ok=True)
    os.makedirs(f"{d}/work/graphs", exist_ok=True)
    os.makedirs(f"{d}/bin", exist_ok=True)
    write(f"{d}/work/data/fastas/in.fa", ">A#1#c\nACGT\n>B#1#c\nACGT\n")
    write(f"{d}/work/data/fastas/in.fa.fai", "A#1#c\t4\t7\t4\t5\nB#1#c\t4\t18\t4\t5\n")
    write(f"{d}/pggb.sif", "stub\n")
    write(f"{d}/vg.sif", "stub\n")
    write(f"{d}/bin/singularity", r'''
        #!/usr/bin/env bash
        echo "$*" >> "$STUB_LOG"
        for x in "$@"; do
          [[ "$x" == --version ]] && { echo "pggb v0.7.4-stub"; exit 0; }
          [[ "$x" == version ]] && { echo 'vg version v1.69.0 "stub"'; exit 0; }
          if [[ "$x" == deconstruct ]]; then
            printf '##fileformat=VCFv4.2\n##contig=<ID=r,length=9>\n'
            printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
            [[ -n "${STUB_FAIL:-}" ]] && exit 1
            printf 'r\t3\t>1>4\tA\tG\t60\t.\tLV=0\tGT\t1\n'
            exit 0
          fi
        done
        exit 0
        ''', mode=0o755)
    env = {"HOME": os.environ.get("HOME", "/tmp"),
           "PATH": f"{d}/bin:/usr/bin:/bin",
           "STUB_LOG": f"{d}/calls.log",
           "MTB_ENV_FILE": os.path.join(ROOT, "config/project_env.sh"),
           "MTB_WORK": f"{d}/work",
           "MTB_PGGB_SIF": f"{d}/pggb.sif",
           "MTB_VG_SIF": f"{d}/vg.sif"}
    env.update(extra or {})
    return env


def run_pggb(d, *args, extra=None):
    r = subprocess.run(["bash", "bin/pggb_build.sh", "in.fa", *args],
                       env=stub_env(d, extra), capture_output=True, text=True)
    calls = open(f"{d}/calls.log").read().splitlines() \
        if os.path.exists(f"{d}/calls.log") else []
    pggb = [c for c in calls if "/usr/local/bin/pggb -i" in c]
    return r, pggb


def provenance(d, name):
    p = f"{d}/work/graphs/{name}/graph_provenance.tsv"
    return dict(l.rstrip("\n").split("\t", 1) for l in open(p))


# ---------------------------------------------------------------- PGB-1

class PggbCleanBuild(unittest.TestCase):
    """PGB-1: the sidecar was written into the output directory before the
    non-empty check, so every fresh build stopped unless resume was on."""

    def test_clean_build_runs_without_resume(self):
        with tempfile.TemporaryDirectory() as d:
            r, pggb = run_pggb(d, "CX333.s10k.k23.K15")
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(len(pggb), 1)
            self.assertNotRegex(pggb[0], r" -r( |$)")
            prov = provenance(d, "CX333.s10k.k23.K15")
            self.assertIn("finished", prov)
            self.assertEqual(prov["resume"], "0")

    def test_non_empty_directory_still_refused(self):
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(f"{d}/work/graphs/G")
            write(f"{d}/work/graphs/G/old.gfa", "x\n")
            r, pggb = run_pggb(d, "G")
            self.assertEqual(r.returncode, 1)
            self.assertIn("not empty", r.stderr)
            self.assertEqual(pggb, [])
            # and nothing was written into it
            self.assertEqual(os.listdir(f"{d}/work/graphs/G"), ["old.gfa"])

    def test_resume_is_opt_in(self):
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(f"{d}/work/graphs/G")
            write(f"{d}/work/graphs/G/old.gfa", "x\n")
            r, pggb = run_pggb(d, "G", extra={"MTB_PGGB_RESUME": "1"})
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertRegex(pggb[0], r" -r$")
            self.assertEqual(provenance(d, "G")["resume"], "1")


# ---------------------------------------------------------------- PGB-2

class PggbSettings(unittest.TestCase):
    """PGB-2: defaults were -k 51 -K 21 on 112 threads; CX333 was -k 23 -K 15
    on 48. Every effective setting goes into the provenance file, sparse
    mapping is an explicit option, and a directory named for settings must
    match them."""

    def test_defaults_are_cx333(self):
        with tempfile.TemporaryDirectory() as d:
            r, pggb = run_pggb(d, "CX333.s10k.k23.K15")
            self.assertEqual(r.returncode, 0, r.stderr)
            args = pggb[0].split("/usr/local/bin/pggb ", 1)[1]
            self.assertEqual(
                args, "-i /data/fastas/in.fa -o /graphs/CX333.s10k.k23.K15 -n 2 "
                      "-t 48 -s 10000 -l 30000 -p 95 -k 23 -K 15 --keep-temp-files")
            prov = provenance(d, "CX333.s10k.k23.K15")
            for k, v in (("segment_length", "10000"), ("block_length", "30000"),
                         ("map_pct_id", "95"), ("min_match_len", "23"),
                         ("mash_kmer", "15"), ("sparse_map", "off"),
                         ("threads", "48"), ("n_haplotypes", "2")):
                self.assertEqual(prov[k], v, k)
            self.assertEqual(prov["pggb_args"], args)
            self.assertRegex(prov["pggb_container_sha256"], r"^[0-9a-f]{64}$")
            self.assertRegex(prov["input_fasta_sha256"], r"^[0-9a-f]{64}$")
            self.assertEqual(prov["pggb_version"], "pggb v0.7.4-stub")
        self.assertRegex(open("bin/pggb_build.sh").read(), r"(?m)^#SBATCH -c 48$")

    def test_sparse_map_is_explicit_and_recorded(self):
        with tempfile.TemporaryDirectory() as d:
            r, pggb = run_pggb(d, "T.s10k.k23.K15", "-x", "auto")
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn(" -x auto ", pggb[0])
            self.assertEqual(provenance(d, "T.s10k.k23.K15")["sparse_map"], "auto")

    def test_overrides_fold_in_once(self):
        # appended overrides used to leave both values on the command line
        with tempfile.TemporaryDirectory() as d:
            r, pggb = run_pggb(d, "O.s10k.k51.K21", "-k", "51", "-K", "21",
                               "-c", "2", "-t", "8")
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(pggb[0].count(" -k "), 1)
            self.assertIn(" -k 51 -K 21 ", pggb[0])
            self.assertIn(" -t 8 ", pggb[0])
            self.assertTrue(pggb[0].endswith("--keep-temp-files -c 2"))
            prov = provenance(d, "O.s10k.k51.K21")
            self.assertEqual((prov["min_match_len"], prov["mash_kmer"],
                              prov["threads"], prov["extra_args"]),
                             ("51", "21", "8", "-c 2"))

    def test_name_must_match_settings(self):
        with tempfile.TemporaryDirectory() as d:
            r, pggb = run_pggb(d, "Q.s10k.k51.K21")
            self.assertEqual(r.returncode, 1)
            self.assertIn("k51 vs -k 23", r.stderr)
            self.assertEqual(pggb, [])
            self.assertFalse(os.path.exists(f"{d}/work/graphs/Q.s10k.k51.K21"))


# ---------------------------------------------------------------- PGB-9

class ReadersTakeCollapsed(unittest.TestCase):
    """PGB-9: every production reader defaulted to the decomposed VCF."""

    TEXT_DEFAULTS = ["assoc/bin/add_outgroup.py", "bin/panel_polarity.py",
                     "bin/sv_intervals.py", "insgt/bin/insertion_contigs.py",
                     "accessory/bin/merge_catalogues.py"]

    def test_no_reader_defaults_to_decomposed(self):
        for p in self.TEXT_DEFAULTS + ["bin/p5_svgt.sh", "bin/vcf_split_classes.sh"]:
            s = open(p).read()
            self.assertNotRegex(s, r'(default=|= ?"|:-)[^\n]*all_variants\.decomposed',
                                p)
            # the graph's collapsed file, or the build's copy of it
            # (cleanup2: production readers have no CX333 default and are
            # passed <build>/assets/graph_collapsed.vcf.gz)
            self.assertTrue("all_variants.collapsed.vcf.gz" in s
                            or "graph_collapsed.vcf.gz" in s, p)

    @staticmethod
    def graph(d, records, samples=("S1", "S2", "S3")):
        os.makedirs(os.path.dirname(os.path.join(d, COLLAPSED)), exist_ok=True)
        with gzip.open(os.path.join(d, COLLAPSED), "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\t"
                     "FILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
            for pos, ref, alt, gts in records:
                fh.write(f"c\t{pos}\t.\t{ref}\t{alt}\t60\t.\t.\tGT\t"
                         + "\t".join(gts) + "\n")

    def test_python_readers_run_on_the_collapsed_default(self):
        # Only the collapsed file exists, as after the next build's cleanup:
        # each reader must find it with no --graph-vcf.
        dele = "A" + "C" * 80
        ins = "G" + "".join("ACGT"[(i * 7) % 4] for i in range(400))
        recs = [(100, dele, "A", ["1", "1", "0"]),
                (5000, "G", ins, ["1", "0", "1"])]
        with tempfile.TemporaryDirectory() as d:
            self.graph(d, recs)
            r = subprocess.run([py(), os.path.join(ROOT, "bin/sv_intervals.py"),
                                "--graph-vcf", COLLAPSED,
                                "--is6110-gff", "", "--out", "iv.tsv"],
                               cwd=d, capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            iv = list(csv.DictReader(open(f"{d}/iv.tsv"), delimiter="\t"))
            self.assertEqual([(x["start"], x["end"], x["ref_carriers"]) for x in iv],
                             [("101", "180", "S1,S2")])
            r = subprocess.run([py(), os.path.join(ROOT, "insgt/bin/insertion_contigs.py"),
                                "--out-fasta", "i.fa", "--out-table", "i.tsv"],
                               cwd=d, capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            t = list(csv.DictReader(open(f"{d}/i.tsv"), delimiter="\t"))
            self.assertEqual([(x["h37rv_pos"], x["length"], x["panel_carriers"])
                              for x in t], [("5000", "400", "S1,S3")])
            write(f"{d}/loci.tsv", "locus_id\tpos\tklass\trep_len\tn_alleles\t"
                  "carriers_any\tcarrier_frac\th37rv_cov\tnovelty\n"
                  "ACC_5000\t5000\tpolymorphic\t400\t1\t2\t0.6\t0\tnovel\n")
            r = subprocess.run([py(), os.path.join(ROOT, "accessory/bin/merge_catalogues.py"),
                                "--loci", "loci.tsv", "--graph-vcf", COLLAPSED,
                                "--out", "c.tsv", "--out-fasta", "c.fa"],
                               cwd=d, capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            c = list(csv.DictReader(open(f"{d}/c.tsv"), delimiter="\t"))
            self.assertEqual((c[0]["graph_len"], c[0]["panel_carriers"]),
                             ("400", "S1,S3"))

    def test_merge_catalogues_on_collapsed_picks_the_most_carried_allele(self):
        # Why the input matters: decomposed, a 400 bp allele carried by three
        # genomes is three records of one carrier each, and loses to a 300 bp
        # allele with two carriers in one record. Collapsed, it is one record
        # with three carriers and wins.
        a400 = "G" + "".join("ACGT"[(i * 7) % 4] for i in range(400))
        a300 = "G" + "".join("ACGT"[(i * 3) % 4] for i in range(300))
        smp = ("S1", "S2", "S3", "S4", "S5")
        dec = [(5000, "G", a400, ["1", "0", "0", "0", "0"]),
               (5000, "G", a400, ["0", "1", "0", "0", "0"]),
               (5000, "G", a400, ["0", "0", "1", "0", "0"]),
               (5000, "G", a300, ["0", "0", "0", "1", "1"])]
        col = [(5000, "G", a300, ["0", "0", "0", "1", "1"]),
               (5000, "G", a400, ["1", "1", "1", "0", "0"])]
        got = {}
        for name, recs in (("dec", dec), ("col", col)):
            with tempfile.TemporaryDirectory() as d:
                self.graph(d, recs, smp)
                write(f"{d}/loci.tsv", "locus_id\tpos\tklass\trep_len\tn_alleles\t"
                      "carriers_any\tcarrier_frac\th37rv_cov\tnovelty\n"
                      "ACC_5000\t5000\tpolymorphic\t400\t2\t5\t1\t0\tnovel\n")
                subprocess.run([py(), os.path.join(ROOT, "accessory/bin/merge_catalogues.py"),
                                "--loci", "loci.tsv", "--graph-vcf", COLLAPSED,
                                "--out", "c.tsv", "--out-fasta", "c.fa"],
                               cwd=d, capture_output=True, text=True, check=True)
                got[name] = list(csv.DictReader(open(f"{d}/c.tsv"),
                                                delimiter="\t"))[0]["graph_len"]
        self.assertEqual(got, {"dec": "300", "col": "400"})


# ---------------------------------------------------------------- PGB-9 sync

class SyncKeepsGraphVcf(unittest.TestCase):
    """PGB-9: sync_back.sh excluded the decomposed file as superseded while
    production read it; the product and provenance must always be kept."""

    def test_mirror_keeps_variant_files(self):
        s = open("bin/sync_back.sh").read()
        m = re.search(r"(?ms)^    EXCLUDES=\(\n(.*?)^    \)\n", s)
        self.assertIsNotNone(m)
        names = ["all_variants.collapsed.vcf.gz", "all_variants.collapsed.vcf.gz.tbi",
                 "all_variants.nolab.vcf.gz", "all_variants.decomposed.vcf.gz",
                 "all_variants.decomposed.vcf.gz.tbi", "graph_provenance.tsv",
                 "x.seqwish.gfa", "variants.vcf"]
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(f"{d}/src/G")
            os.makedirs(f"{d}/dst")
            for n in names:
                write(f"{d}/src/G/{n}", "x\n")
            script = "EXCLUDES=(\n" + m.group(1) + ")\n" \
                     f'rsync -an --out-format="%n" "${{EXCLUDES[@]}}" {d}/src/ {d}/dst/\n'
            r = subprocess.run(["bash", "-c", script], capture_output=True,
                               text=True, check=True)
        sent = {os.path.basename(l) for l in r.stdout.split()}
        for n in names[:6]:
            self.assertIn(n, sent)
        self.assertNotIn("x.seqwish.gfa", sent)
        self.assertNotIn("variants.vcf", sent)


# ---------------------------------------------------------------- deconstruct

class DeconstructHeader(unittest.TestCase):
    """vg_deconstruct.sh: settings into the header, write under a temporary
    name, and never beside a variants.vcf.gz that vcf_decompose.sh prefers."""

    def setup(self, d):
        os.makedirs(f"{d}/work/graphs/G")
        write(f"{d}/work/graphs/G/g.smooth.final.gfa", "H\tVN:Z:1.0\n")

    def test_header_records_command_version_and_graph(self):
        with tempfile.TemporaryDirectory() as d:
            self.setup(d)
            r = subprocess.run(["bash", "bin/vg_deconstruct.sh", "G"],
                               env=stub_env(d), capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            out = f"{d}/work/graphs/G/variants.vcf"
            lines = open(out).read().splitlines()
            self.assertEqual(sorted(os.listdir(f"{d}/work/graphs/G")),
                             ["g.smooth.final.gfa", "variants.vcf"])
        hdr = [l for l in lines if l.startswith("##MTB_deconstruct")]
        self.assertEqual(len(hdr), 2)
        self.assertIn("vg deconstruct -a -t 8 -P GCF_000195955#1#NC_000962.3 "
                      "g.smooth.final.gfa", hdr[0])
        self.assertIn('vg version v1.69.0 "stub"', hdr[0])
        self.assertRegex(hdr[1], r"^##MTB_deconstruct_graph_sha256=[0-9a-f]{64}$")
        chrom = [i for i, l in enumerate(lines) if l.startswith("#CHROM")][0]
        self.assertEqual(lines[chrom - 1], hdr[1])
        self.assertEqual(lines[-1].split("\t")[:5], ["r", "3", ">1>4", "A", "G"])

    def test_failed_run_leaves_no_output(self):
        with tempfile.TemporaryDirectory() as d:
            self.setup(d)
            r = subprocess.run(["bash", "bin/vg_deconstruct.sh", "G"],
                               env=stub_env(d, {"STUB_FAIL": "1"}),
                               capture_output=True, text=True)
            self.assertNotEqual(r.returncode, 0)
            self.assertEqual(os.listdir(f"{d}/work/graphs/G"), ["g.smooth.final.gfa"])

    def test_refuses_beside_stale_bgzipped_vcf(self):
        with tempfile.TemporaryDirectory() as d:
            self.setup(d)
            write(f"{d}/work/graphs/G/variants.vcf.gz", "old\n")
            r = subprocess.run(["bash", "bin/vg_deconstruct.sh", "G"],
                               env=stub_env(d), capture_output=True, text=True)
            self.assertEqual(r.returncode, 1)
            self.assertIn("variants.vcf.gz exists", r.stderr)


# ---------------------------------------------------------------- PGB-12

FAKE_BCFTOOLS = r'''
#!/usr/bin/env python3
"""bcftools query -l / -f on a plain-text VCF, enough for snp_nonredundant.py."""
import sys
a = sys.argv[1:]
vcf = a[-1]
lines = [l.rstrip("\n").split("\t") for l in open(vcf) if not l.startswith("##")]
if "-l" in a:
    print("\n".join(lines[0][9:])); sys.exit(0)
fmt = a[a.index("-f") + 1]
for f in lines[1:]:
    if fmt.startswith("[%GT"):
        print("\t".join(f[9:]) + "\t")
    else:
        print("\t".join([f[0], f[1], f[3], f[4]] + f[9:]))
'''


class NonRedundantInput(unittest.TestCase):
    """PGB-12: duplicate-key input inflated distances; multiallelic GT 1 and 2
    were flattened; without --rank-by the representative was the largest
    accession."""

    def run_nr(self, d, records, rank=True):
        samples = ["G1", "G2", "G3"]
        with open(f"{d}/in.vcf", "w") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\t"
                     "FILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
            for pos, ref, alt, gts in records:
                fh.write(f"c\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t"
                         + "\t".join(gts) + "\n")
        write(f"{d}/acc.txt", "\n".join(samples) + "\n")
        write(f"{d}/bcftools", FAKE_BCFTOOLS, mode=0o755)
        write(f"{d}/rank.tsv", "G1\t3\nG2\t1\nG3\t2\n")
        cmd = [py(), "bin/snp_nonredundant.py", "--vcf", f"{d}/in.vcf",
               "--accessions", f"{d}/acc.txt", "--bcftools", f"{d}/bcftools",
               "--thresholds", "0", "1", "--pick-threshold", "0",
               "--out", f"{d}/out.tsv"]
        if rank:
            cmd += ["--rank-by", f"{d}/rank.tsv"]
        return subprocess.run(cmd, capture_output=True, text=True)

    def test_duplicate_keys_refused(self):
        with tempfile.TemporaryDirectory() as d:
            r = self.run_nr(d, [(10, "A", "C", ["1", "0", "0"]),
                                (10, "A", "C", ["0", "1", "0"])])
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("repeat a (CHROM, POS, REF, ALT) key", r.stderr)

    def test_multiallelic_alleles_differ(self):
        with tempfile.TemporaryDirectory() as d:
            r = self.run_nr(d, [(10, "A", "C,G", ["1", "2", "1"])])
            self.assertEqual(r.returncode, 0, r.stderr)
            out = list(csv.DictReader(open(f"{d}/out.tsv"), delimiter="\t"))
        rep = {x["accession"]: x["representative"] for x in out}
        # G1 and G3 carry C, G2 carries G: two clusters at distance 0
        self.assertEqual(rep, {"G1": "G1", "G3": "G1", "G2": "G2"})

    def test_rank_by_required(self):
        with tempfile.TemporaryDirectory() as d:
            r = self.run_nr(d, [(10, "A", "C", ["1", "0", "0"])], rank=False)
        self.assertEqual(r.returncode, 2)
        self.assertIn("--rank-by", r.stderr)


if __name__ == "__main__":
    unittest.main(warnings="ignore")
