#!/usr/bin/env python3
"""Screen assemblies for reference-guided provenance and short-read limitation.

This is the screen that removed 150 of 484 RefSeq "Complete Genome" assemblies
from the panel. It was developed interactively and is written up here so it is
reproducible; the finding is documented in the artifact linked from STATUS.md.

The problem it solves: a substantial share of public MTBC "complete genomes" are
consensus sequences built by mapping short reads to H37Rv, or short-read
assemblies scaffolded on it. They type correctly by SNP -- the SNPs are the
sample's own -- so lineage assignment succeeds and nothing looks wrong. But their
structural content is the reference's, which makes them worse than useless for a
pangenome: they contribute phylogenetically inert structure while counting as
taxa.

Crucially, the screens that already existed cannot see this:

  * collinearity against the reference passes by construction -- a
    reference-scaffolded assembly IS collinear with the reference;
  * ambiguous-base content is zero (89 of 94 carried no Ns at all);
  * contig count is 1, and assembly level reads "Complete Genome";
  * SNP-based lineage typing succeeds.

Two measurements separate them.

**Insertion content (sequence only, no metadata).** A reference-guided assembly
cannot contain sequence absent from its template. Counting insertions of >= 50 bp
relative to the reference gave, on the 484-genome panel: median 0 (max 1) for
reference-structured assemblies against median 56 (min 5) for long-read ones.
A threshold of < 5 excluded 118 genomes with zero false exclusions. This is the
primary screen and needs nothing but the assembly and a reference.

**Sequencing technology (NCBI metadata).** 93 of 94 reference-structured genomes
reported short-read platforms only, none long-read or hybrid; the retained set was
217 long-read, 100 hybrid, 56 short-read-only. Requiring long-read or hybrid data
for a "Complete Genome" claim catches all 150 including the intermediate class
that passes the insertion test but still under-reports insertions (median 3 novel
IS6110 against 11 for long-read assemblies).

Report both, and treat the sequence measurement as authoritative: a composition
metric should never be used to overturn an exclusion, only to support one. An
assembly can carry ample novel sequence and still be missing regions and carry
errors in genes -- that is a different failure mode this screen does not detect.
"""
import argparse, csv, json, os, subprocess, sys, time

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession"
LONG = ("pacbio", "nanopore", "oxford", "sequel", "rs ii", "promethion",
        "gridion", "minion", "hifi", "revio")
SHORT = ("illumina", "ion torrent", "454", "solid", "bgi",
         "complete genomics", "mgi")
ALIGNER = ("bwa", "bowtie", "samtools", "mpileup", "smalt", "novoalign",
           "minimap", "bcftools", "tanoti", "consensus")


def tech_class(t):
    tl = (t or "").lower()
    if not tl:
        return "unknown"
    lr = any(k in tl for k in LONG)
    sr = any(k in tl for k in SHORT)
    return ("hybrid" if lr and sr else "long_read" if lr
            else "short_read" if sr else "unknown")


def fetch_metadata(versioned, batch=20, pause=0.35):
    """assembly_method and sequencing_tech from NCBI Datasets"""
    out = {}
    items = list(versioned.items())
    for i in range(0, len(items), batch):
        chunk = items[i:i + batch]
        url = f"{API}/{','.join(v for _, v in chunk)}/dataset_report?page_size=100"
        raw = subprocess.run(["curl", "-s", "--max-time", "90", url],
                             capture_output=True, text=True).stdout
        try:
            d = json.loads(raw)
        except ValueError:
            continue
        for r in d.get("reports", []):
            acc = r.get("accession", "").split(".")[0]
            ai = r.get("assembly_info", {})
            out[acc] = (ai.get("assembly_method", "") or "",
                        ai.get("sequencing_tech", "") or "")
        time.sleep(pause)
    return out


def insertion_content(asm, ref, mm2, k8, paftools, min_block=5000, min_ins=50):
    """(inserted_bp, deleted_bp, n_insertions >= min_ins)"""
    p1 = subprocess.Popen([mm2, "-cx", "asm5", "--cs", "-t", "2", ref, asm],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    srt = subprocess.Popen(["sort", "-k6,6", "-k8,8n"], stdin=p1.stdout,
                           stdout=subprocess.PIPE)
    p1.stdout.close()
    res = subprocess.run([k8, paftools, "call", "-L", str(min_block),
                          "-l", str(min_block), "-"],
                         stdin=srt.stdout, capture_output=True, text=True)
    srt.stdout.close()
    ins = dele = nbig = 0
    for line in res.stdout.splitlines():
        f = line.split("\t")
        if f[0] != "V":
            continue
        r = 0 if f[6] == "-" else len(f[6])
        a = 0 if f[7] == "-" else len(f[7])
        if a > r:
            ins += a - r
            if a - r >= min_ins:
                nbig += 1
        else:
            dele += r - a
    return ins, dele, nbig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accessions", required=True)
    ap.add_argument("--assembly-dir", required=True)
    ap.add_argument("--assembly-suffix", default=".dnaA_rotated.fasta")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--assembly-summary",
                    help="NCBI assembly_summary_refseq.txt, to resolve versioned "
                         "accessions for the metadata query; without it only the "
                         "sequence screen runs")
    ap.add_argument("--min-big-ins", type=int, default=5,
                    help="minimum insertions >= 50 bp; below this the assembly "
                         "cannot contain novel sequence and is excluded")
    ap.add_argument("--require-long-read", action="store_true",
                    help="also exclude assemblies reporting no long-read or "
                         "hybrid data, which catches the intermediate class")
    ap.add_argument("--minimap2", required=True)
    ap.add_argument("--k8", required=True)
    ap.add_argument("--paftools", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    accs = [l.strip() for l in open(a.accessions) if l.strip()]
    meta = {}
    if a.assembly_summary:
        ver = {}
        for line in open(a.assembly_summary):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            base = f[0].split(".")[0]
            if base in set(accs):
                ver[base] = f[0]
        print(f"  querying NCBI for {len(ver)} accessions ...", flush=True)
        meta = fetch_metadata(ver)
        print(f"  metadata retrieved for {len(meta)}")

    rows = []
    for i, acc in enumerate(accs, 1):
        asm = os.path.join(a.assembly_dir, acc + a.assembly_suffix)
        if not os.path.exists(asm):
            print(f"  WARNING: no assembly for {acc}", file=sys.stderr)
            continue
        ins, dele, nbig = insertion_content(asm, a.ref, a.minimap2, a.k8, a.paftools)
        m, t = meta.get(acc, ("", ""))
        tc = tech_class(t)
        reasons = []
        if nbig < a.min_big_ins:
            reasons.append("no_novel_sequence")
        if a.require_long_read and tc == "short_read":
            reasons.append("short_read_only")
        if any(k in m.lower() for k in ALIGNER):
            reasons.append("method_names_aligner")
        rows.append(dict(accession=acc, assembly_method=m, sequencing_tech=t,
                         tech_class=tc, inserted_bp=ins, deleted_bp=dele,
                         insertions_ge50=nbig,
                         verdict="EXCLUDE" if reasons else "ok",
                         reasons=";".join(reasons)))
        if i % 50 == 0:
            print(f"    {i}/{len(accs)}", flush=True)

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    ex = [r for r in rows if r["verdict"] == "EXCLUDE"]
    import collections, statistics
    print(f"\n  {len(rows)} assemblies screened, {len(ex)} EXCLUDE")
    for k, v in collections.Counter(r["reasons"] for r in ex).most_common():
        print(f"    {k:<44}{v:>4}")
    keep = [r for r in rows if r["verdict"] == "ok"]
    for nm, s in (("excluded", ex), ("retained", keep)):
        if not s:
            continue
        print(f"  {nm:<10} median insertions >=50bp "
              f"{statistics.median(r['insertions_ge50'] for r in s):>5.0f}"
              f"   median inserted bp "
              f"{statistics.median(r['inserted_bp'] for r in s):>9,.0f}")
    print(f"\n  wrote {a.out}")


if __name__ == "__main__":
    main()
