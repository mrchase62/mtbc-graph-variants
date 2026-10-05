#!/usr/bin/env python3
"""Find insertions that are foreign to the MTBC — vector, plasmid, engineered DNA.

Metadata cannot detect an engineered strain. A ~5 kb pJEB integration in one panel
genome was invisible in its RefSeq record: strain field "Erdman = ATCC 35801",
isolate "na", assembly method Flye, nothing amiss. It was found only by comparing
two assemblies of the same strain, where it was the single difference.

That comparison generalises without needing a matched pair. Take every insertion
of >= --min-len that a genome carries relative to the reference, and ask whether it
aligns anywhere in ANY OTHER panel genome:

  native   IS6110, prophage, RD-region sequence and duplications all have
           homologues elsewhere in the panel, because other genomes carry the
           same mobile elements or the source locus.
  foreign  vector backbone, selection markers and synthetic barcodes align
           nowhere, because no MTBC genome contains them.

Self-hits are excluded using the PanSN target names, so an insert is never
credited to the genome it came from.

The screen is deliberately run at >= 1 kb: below that, short native duplications
and repeat fragments produce enough unmatched sequence to bury the signal.
"""
import argparse, collections, csv, os, subprocess, sys


def extract_inserts(acc, asm, ref, mm2, k8, paftools, min_len, min_block, out_fh):
    p1 = subprocess.Popen([mm2, "-cx", "asm5", "--cs", "-t", "2", ref, asm],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    srt = subprocess.Popen(["sort", "-k6,6", "-k8,8n"], stdin=p1.stdout,
                           stdout=subprocess.PIPE)
    p1.stdout.close()
    out = subprocess.run([k8, paftools, "call", "-L", str(min_block),
                          "-l", str(min_block), "-"],
                         stdin=srt.stdout, capture_output=True, text=True)
    srt.stdout.close()
    n = 0
    for line in out.stdout.splitlines():
        f = line.split("\t")
        if f[0] != "V":
            continue
        r = "" if f[6] == "-" else f[6]
        a = "" if f[7] == "-" else f[7]
        if len(a) - len(r) < min_len:
            continue
        n += 1
        out_fh.write(f">{acc}|{f[2]}|{len(a)-len(r)}\n{a}\n")
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accessions", required=True)
    ap.add_argument("--assembly-dir", required=True)
    ap.add_argument("--assembly-suffix", default=".dnaA_rotated.fasta")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--panel-fasta", required=True,
                    help="PanSN-named concatenation used as the homology background")
    ap.add_argument("--min-len", type=int, default=1000)
    ap.add_argument("--min-block", type=int, default=5000)
    ap.add_argument("--max-foreign-frac", type=float, default=0.10,
                    help="an insert with less than this fraction of its length "
                         "aligned to any other genome is called foreign")
    ap.add_argument("--minimap2", required=True)
    ap.add_argument("--k8", required=True)
    ap.add_argument("--paftools", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    os.makedirs(a.workdir, exist_ok=True)
    accs = [l.strip() for l in open(a.accessions) if l.strip()]
    ins_fa = os.path.join(a.workdir, "inserts.fa")

    if not os.path.exists(ins_fa) or os.path.getsize(ins_fa) == 0:
        print(f"  extracting insertions >= {a.min_len} bp from {len(accs)} genomes ...")
        with open(ins_fa, "w") as fh:
            for i, acc in enumerate(accs, 1):
                asm = os.path.join(a.assembly_dir, acc + a.assembly_suffix)
                if not os.path.exists(asm):
                    continue
                extract_inserts(acc, asm, a.ref, a.minimap2, a.k8, a.paftools,
                                a.min_len, a.min_block, fh)
                if i % 50 == 0:
                    print(f"    {i}/{len(accs)}", flush=True)
    n_ins = sum(1 for l in open(ins_fa) if l.startswith(">"))
    print(f"  {n_ins} insertions >= {a.min_len} bp across the panel")
    if n_ins == 0:
        sys.exit("no insertions extracted")

    paf = os.path.join(a.workdir, "inserts_vs_panel.paf")
    if not os.path.exists(paf) or os.path.getsize(paf) == 0:
        print("  aligning inserts against the panel background ...")
        with open(paf, "w") as fh:
            # Secondary alignments are REQUIRED. Each insert's best hit is
            # always to the genome it came from, so with --secondary=no the
            # self-hit is the only alignment reported and excluding it leaves
            # nothing -- which flagged 3,169 of 4,936 inserts as foreign across
            # 331 of 334 genomes. -N 50 -p 0.1 lets homologues in other genomes
            # be reported alongside the self-hit.
            subprocess.run([a.minimap2, "-cx", "asm10", "-t", str(a.threads),
                            "-I", "8G", "--secondary=yes", "-N", "50", "-p", "0.1",
                            a.panel_fasta, ins_fa],
                           stdout=fh, stderr=subprocess.DEVNULL, check=True)

    # sum aligned bases per insert, EXCLUDING hits to the insert's own genome
    hit = collections.Counter()
    qlen = {}
    for line in open(paf):
        f = line.split("\t")
        q, tgt = f[0], f[5]
        qlen[q] = int(f[1])
        src = q.split("|")[0]
        tsrc = tgt.split("#")[0]
        if tsrc == src:
            continue
        hit[q] += int(f[9])

    rows = []
    for line in open(ins_fa):
        if not line.startswith(">"):
            continue
        q = line[1:].strip()
        acc, pos, ln = q.split("|")
        ln = int(ln)
        aligned = hit.get(q, 0)
        frac = aligned / max(ln, 1)
        rows.append(dict(accession=acc, ref_pos=int(pos), insert_len=ln,
                         aligned_elsewhere=aligned, frac=round(frac, 4),
                         verdict="FOREIGN" if frac < a.max_foreign_frac else "native"))
    rows.sort(key=lambda r: (r["verdict"] != "FOREIGN", -r["insert_len"]))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    fo = [r for r in rows if r["verdict"] == "FOREIGN"]
    print(f"\n  {len(fo)} of {len(rows)} insertions are FOREIGN "
          f"(<{100*a.max_foreign_frac:.0f}% aligned to any other genome)")
    bg = collections.Counter(r["accession"] for r in fo)
    print(f"  affecting {len(bg)} genomes\n")
    print(f"  {'accession':<16}{'n foreign':>10}{'total bp':>10}{'largest':>9}  positions")
    for acc, n in bg.most_common():
        sub = [r for r in fo if r["accession"] == acc]
        tot = sum(r["insert_len"] for r in sub)
        big = max(r["insert_len"] for r in sub)
        pos = ", ".join(str(r["ref_pos"]) for r in sorted(sub, key=lambda r: -r["insert_len"])[:4])
        print(f"  {acc:<16}{n:>10}{tot:>10,}{big:>9,}  {pos}")
    print(f"\n  wrote {a.out}")


if __name__ == "__main__":
    main()
