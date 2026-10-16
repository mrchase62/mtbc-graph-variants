"""bin/cohort_gate.py: the cohort acceptance gate fails the problems that
reached new1k_p200 unseen (SV-EVAL-1..4) and passes a clean cohort."""
import gzip, os, subprocess, sys, tempfile, unittest

R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GATE = os.path.join(R, "bin", "cohort_gate.py")
HDR = ("##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">\n"
       "##FORMAT=<ID=ST,Number=1,Type=String,Description=\"\">\n")


def rec(pos, info, cells):
    return f"NC_000962.3\t{pos}\t.\tA\tG\t.\tPASS\t{info}\tGT:ST\t" + "\t".join(cells) + "\n"


def write(d, rows, n=20):
    p = os.path.join(d, "m.vcf.gz")
    with gzip.open(p, "wt") as fh:
        fh.write(HDR + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
                 "\t".join(f"s{i}" for i in range(n)) + "\n")
        for r in rows: fh.write(r)
    return p


def build(d, n_loci):
    os.makedirs(os.path.join(d, "b", "assets"))
    with open(os.path.join(d, "b", "assets", "accessory_catalogue.tsv"), "w") as fh:
        fh.write("locus_id\n" + "".join(f"ACC_{i}\n" for i in range(n_loci)))
    return os.path.join(d, "b")


def clean_rows(n=20, loci=4):
    alt = ["1:ALT"] * 5 + ["0:REF"] * (n - 5)
    rows = [rec(100 + i, "CLASS=small;KIND=SNP;FRAME=h37rv", alt) for i in range(50)]
    rows += [rec(10000, "CLASS=sv;SVTYPE=INS;SVLEN=80;SVSOURCE=small_allele", alt)]
    rows += [rec(20000, "CLASS=is6110;FRAME=h37rv", alt)]
    rows += [rec(30000 + i, f"CLASS=accessory_presence;ACCLOCUS=ACC_{i}", alt) for i in range(loci)]
    return rows


class Gate(unittest.TestCase):
    def run_gate(self, d, rows, loci=4):
        b = build(d, loci)
        r = subprocess.run([sys.executable, GATE, "--vcf", write(d, rows), "--build", b,
                            "--out", os.path.join(d, "g.tsv")], capture_output=True, text=True)
        return r

    def test_clean_cohort_passes(self):
        with tempfile.TemporaryDirectory() as d:
            r = self.run_gate(d, clean_rows())
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)

    def test_missing_accessory_loci_fail(self):
        with tempfile.TemporaryDirectory() as d:
            r = self.run_gate(d, clean_rows(loci=1), loci=4)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("FAIL  accessory:loci_written", r.stdout)

    def test_sitewide_is6110_absent_fails(self):
        with tempfile.TemporaryDirectory() as d:
            rows = clean_rows() + [rec(40000 + i, "CLASS=is6110;FRAME=h37rv", ["2:ABSENT"] * 20) for i in range(6)]
            r = self.run_gate(d, rows)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("FAIL  sitewide:is6110", r.stdout)

    def test_long_small_allele_fails(self):
        with tempfile.TemporaryDirectory() as d:
            rows = clean_rows() + [rec(50000, "CLASS=sv;SVTYPE=INS;SVLEN=900;SVSOURCE=small_allele",
                                       ["1:ALT"] * 3 + ["0:REF"] * 17)]
            r = self.run_gate(d, rows)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("FAIL  small_allele:too_long", r.stdout)

    def test_accessory_nocall_jump_fails(self):
        with tempfile.TemporaryDirectory() as d:
            rows = clean_rows(loci=0) + [rec(30000 + i, f"CLASS=accessory_presence;ACCLOCUS=ACC_{i}",
                                             ["1:ALT"] * 5 + [".:NOCALL"] * 15) for i in range(4)]
            r = self.run_gate(d, rows)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("FAIL  nocall:accessory_presence", r.stdout)


if __name__ == "__main__":
    unittest.main()
