#!/usr/bin/env python3
"""Screen assemblies for reference-guided provenance and short-read limitation.

This is the screen that removed 150 of 484 RefSeq "Complete Genome" assemblies
from the panel. It was developed interactively; this script now implements the
rule that production applied (data/qc/assembly_provenance.tsv), so the screen
is reproducible.

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

THE RULE (production's, audit PGB-4), in this order:

  REFERENCE_STRUCTURE  the genome's IS6110 profile is identical to H37Rv's:
                       the same number of copies, every H37Rv copy present at
                       its H37Rv site, and no copy anywhere else. A genome of
                       another strain cannot share H37Rv's exact IS6110
                       profile; a consensus built on H37Rv inherits it. 94 of
                       484 in production, all short-read or unknown
                       technology.
  SHORT_READ_SUSPECT   otherwise, NCBI sequencing technology is short-read only
                       (tech_class short_read). 56 of 484.
  ok                   everything else, INCLUDING unknown technology (17
                       genomes passed this way in production). Unknown is
                       written to the `review` column so it is not silent.

An earlier version of this script excluded on "< 5 insertions of >= 50 bp"
(--min-big-ins) instead. That is not the rule that built CX333: the
production table holds no insertion count. The count is still measured and
reported (`insertions_ge50`), and it feeds `tech_from_sequence`, but it does
not decide the flag.

IS6110 PROFILE. Copies are found by aligning the canonical IS6110
(--is6110) to the assembly with minimap2 (an alignment covering >= 50% of the
element; overlapping hits merged). Each copy's two 300 bp flanks are aligned to
the reference (asm5, primary alignment, MAPQ >= 1), and the reference coordinate
next to the element is the copy's site. A copy whose placed flanks all lie
within 50 bp of the same reference copy is that reference copy "kept" (36
of production's 94 reference-structured genomes have a flank 26 bp off, an
alignment-end effect, so 10 bp is too tight); any other placed copy is
"novel"; a copy with neither flank placed counts only in the total. These
counts are this script's own; production's is6110_ref_kept / is6110_novel
columns were made by an unrecorded method and do not match them copy for
copy. The identical-profile call, the one the rule uses, agrees with
production on 483 of 484 genomes; the exception, GCF_002886775, is called
REFERENCE_STRUCTURE here and SHORT_READ_SUSPECT in production (excluded
either way). The
reference's own copies come from the same procedure run on the reference, so
the two sides are measured identically.

TECHNOLOGY FROM SEQUENCE (a flag, not a rule). Metadata are unreliable (SY-1,
GCF_050259585, says PacBio and was assembled with SOAPdenovo). From the
sequence alone:
  reference_structured  identical IS6110 profile
  short_read_like       fewer than --min-big-ins (5) insertions >= 50 bp
                        against the reference: an assembly that cannot carry
                        novel sequence (median 0 for reference-structured,
                        min 5 for long-read assemblies in the 484)
  long_read_like        otherwise
`tech_conflict` records disagreements: metadata long-read or hybrid with
sequence short-read-like or reference-structured; a long-read label with a
short-read assembler named and no long-read assembler; metadata short-read
with a long-read-like sequence. Uncorrected long-read error (private
homopolymer indels) is the separate check in variant_counts.py ->
snp_outlier_screen.py.

Do not use the insertion count to overturn an exclusion. An assembly can carry
ample novel sequence and still be missing regions and carry errors in genes --
that is a different failure mode this screen does not detect.
"""
import argparse, csv, json, os, re, subprocess, sys, tempfile, time

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession"
LONG = ("pacbio", "nanopore", "oxford", "sequel", "rs ii", "rsii", "promethion",
        "gridion", "minion", "hifi", "revio", "smrt")
SHORT = ("illumina", "ion torrent", "454", "solid", "bgi",
         "complete genomics", "mgi", "miseq", "hiseq", "nextseq", "novaseq",
         "iseq", "dnbseq", "genome analyzer")
# whole-word only: "ont" must not match inside another word
LONG_WORDS = ("ont",)
ALIGNER = ("bwa", "bowtie", "samtools", "mpileup", "smalt", "novoalign",
           "minimap", "bcftools", "tanoti", "consensus")
SR_ASSEMBLER = ("soapdenovo", "spades", "velvet", "abyss", "skesa", "shovill",
                "newbler", "a5-miseq", "megahit")
LR_ASSEMBLER = ("flye", "canu", "hgap", "falcon", "unicycler", "miniasm",
                "raven", "wtdbg", "nextdenovo", "necat", "hifiasm", "trycycler",
                "shasta", "celera", "smrt", "ratatosk", "redbean")
REF_SITE_TOL = 50
FLANK = 300


def tech_class(t):
    tl = (t or "").lower()
    if not tl:
        return "unknown"
    lr = any(k in tl for k in LONG) or \
        any(re.search(rf"\b{k}\b", tl) for k in LONG_WORDS)
    sr = any(k in tl for k in SHORT)
    return ("hybrid" if lr and sr else "long_read" if lr
            else "short_read" if sr else "unknown")


def provenance_flag(identical_profile, tc):
    """production's rule; see the module docstring"""
    if identical_profile:
        return "REFERENCE_STRUCTURE"
    if tc == "short_read":
        return "SHORT_READ_SUSPECT"
    return "ok"


def tech_from_sequence(identical_profile, nbig, min_big_ins):
    if identical_profile:
        return "reference_structured"
    if nbig is None:
        return ""
    return "short_read_like" if nbig < min_big_ins else "long_read_like"


def tech_conflicts(tc, seq_tech, method):
    m = (method or "").lower()
    out = []
    if tc in ("long_read", "hybrid") and seq_tech in ("short_read_like",
                                                      "reference_structured"):
        out.append(f"metadata_{tc}_but_sequence_{seq_tech}")
    if tc == "short_read" and seq_tech == "long_read_like":
        out.append("metadata_short_read_but_sequence_long_read_like")
    if tc in ("long_read", "hybrid") and any(k in m for k in SR_ASSEMBLER) \
            and not any(k in m for k in LR_ASSEMBLER):
        out.append("short_read_assembler")
    return out


def read_fasta(path):
    out, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        out[name] = "".join(buf).upper()
    return out


def merge_hits(hits):
    """[(contig, start, end)] -> merged intervals per contig, sorted"""
    out = []
    for c, s, e in sorted(hits):
        if out and out[-1][0] == c and s <= out[-1][2]:
            out[-1][2] = max(out[-1][2], e)
        else:
            out.append([c, s, e])
    return [tuple(x) for x in out]


def is6110_copies(asm, is6110, mm2, min_frac=0.5):
    """IS6110 copies in an assembly: [(contig, start, end)], 0-based half-open"""
    res = subprocess.run([mm2, "-c", "--secondary=yes", "-N", "100", "-p", "0.1",
                          asm, is6110], capture_output=True, text=True)
    hits = []
    for line in res.stdout.splitlines():
        f = line.split("\t")
        qlen, qs, qe = int(f[1]), int(f[2]), int(f[3])
        if qe - qs >= min_frac * qlen:
            hits.append((f[5], int(f[7]), int(f[8])))
    return merge_hits(hits)


def flank_site(strand, side, tstart, tend):
    """reference coordinate next to the element, from one flank's alignment.
    side L = the flank left of the element in the assembly."""
    return tend if (side == "L") == (strand == "+") else tstart


def classify_sites(n_copies, sites_by_copy, ref_sites, tol=REF_SITE_TOL):
    """(total, ref_kept, novel). sites_by_copy: {copy index: [reference
    positions from its placed flanks]}; ref_sites: [(start, end)] of the
    reference's own copies."""
    kept, novel = set(), 0
    for i in range(n_copies):
        ps = sites_by_copy.get(i, [])
        hit = None
        # kept only when EVERY placed flank sits at the same reference copy: a
        # copy with one flank at an H37Rv site and the other displaced (a
        # deletion or rearrangement beside the element) is not H37Rv's copy
        for k, (a, b) in enumerate(ref_sites):
            if ps and all(a - tol <= p <= b + tol for p in ps):
                hit = k
                break
        if hit is not None:
            kept.add(hit)
        elif ps:
            novel += 1
    return n_copies, len(kept), novel


def is6110_profile(asm, ref, is6110, ref_sites, mm2, flank=FLANK):
    """(total, ref_kept, novel, identical_to_reference)"""
    seqs = read_fasta(asm)
    copies = is6110_copies(asm, is6110, mm2)
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        for i, (c, s, e) in enumerate(copies):
            sq = seqs[c]
            if s - flank >= 0:
                fh.write(f">{i}_L\n{sq[s - flank:s]}\n")
            if e + flank <= len(sq):
                fh.write(f">{i}_R\n{sq[e:e + flank]}\n")
        fl = fh.name
    try:
        res = subprocess.run([mm2, "-c", "-x", "asm5", ref, fl],
                             capture_output=True, text=True)
    finally:
        os.unlink(fl)
    sites = {}
    for line in res.stdout.splitlines():
        f = line.split("\t")
        if int(f[11]) < 1 or "tp:A:P" not in line:
            continue
        i, side = f[0].split("_")
        sites.setdefault(int(i), []).append(
            flank_site(f[4], side, int(f[7]), int(f[8])))
    total, kept, novel = classify_sites(len(copies), sites, ref_sites)
    identical = (total == len(ref_sites) and kept == len(ref_sites)
                 and novel == 0)
    return total, kept, novel, identical


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
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accessions", required=True)
    ap.add_argument("--assembly-dir", required=True)
    ap.add_argument("--assembly-suffix", default=".dnaA_rotated.fasta")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--is6110", required=True,
                    help="canonical IS6110 FASTA (data/annotation/IS6110.fasta)")
    ap.add_argument("--assembly-summary",
                    help="NCBI assembly_summary_refseq.txt, to resolve versioned "
                         "accessions for the metadata query")
    ap.add_argument("--metadata",
                    help="recorded metadata instead of an NCBI query: TSV with "
                         "accession, assembly_method, sequencing_tech (e.g. a "
                         "previous assembly_provenance.tsv)")
    ap.add_argument("--min-big-ins", type=int, default=5,
                    help="insertions >= 50 bp below which the sequence is "
                         "called short_read_like (a flag, not the rule)")
    ap.add_argument("--no-insertion-count", action="store_true",
                    help="skip the whole-genome alignment that counts insertions "
                         "(the rule does not need it; tech_from_sequence then "
                         "rests on the IS6110 profile only)")
    ap.add_argument("--minimap2", required=True)
    ap.add_argument("--k8", required=True)
    ap.add_argument("--paftools", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    accs = [l.strip() for l in open(a.accessions) if l.strip()]
    meta = {}
    if a.metadata:
        for r in csv.DictReader(open(a.metadata), delimiter="\t"):
            meta[r["accession"]] = (r.get("assembly_method", "") or "",
                                    r.get("sequencing_tech", "") or "")
        print(f"  metadata read from {a.metadata} for "
              f"{sum(1 for x in accs if x in meta)} of {len(accs)}")
    elif a.assembly_summary:
        ver = {}
        want = set(accs)
        for line in open(a.assembly_summary):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            base = f[0].split(".")[0]
            if base in want:
                ver[base] = f[0]
        print(f"  querying NCBI for {len(ver)} accessions ...", flush=True)
        meta = fetch_metadata(ver)
        print(f"  metadata retrieved for {len(meta)}")
    else:
        print("  WARNING: no --metadata or --assembly-summary; technology is "
              "unknown for every genome", file=sys.stderr)

    rc = is6110_copies(a.ref, a.is6110, a.minimap2)
    ref_sites = [(s, e) for _, s, e in rc]
    print(f"  reference IS6110 copies: {len(ref_sites)}")
    if not ref_sites:
        sys.exit("no IS6110 copy found in the reference; check --is6110")

    rows, missing = [], []
    for i, acc in enumerate(accs, 1):
        asm = os.path.join(a.assembly_dir, acc + a.assembly_suffix)
        if not os.path.exists(asm):
            missing.append(acc)
            continue
        tot, kept, novel, ident = is6110_profile(asm, a.ref, a.is6110,
                                                 ref_sites, a.minimap2)
        if a.no_insertion_count:
            ins = dele = nbig = None
        else:
            ins, dele, nbig = insertion_content(asm, a.ref, a.minimap2, a.k8,
                                                a.paftools)
        m, t = meta.get(acc, ("", ""))
        tc = tech_class(t)
        seq = tech_from_sequence(ident, nbig, a.min_big_ins)
        review = tech_conflicts(tc, seq, m)
        if tc == "unknown":
            review.append("unknown_technology")
        if any(k in m.lower() for k in ALIGNER):
            review.append("method_names_aligner")
        rows.append(dict(accession=acc, assembly_method=m, sequencing_tech=t,
                         tech_class=tc, is6110_total=tot, is6110_ref_kept=kept,
                         is6110_novel=novel,
                         identical_profile_ref="yes" if ident else "no",
                         inserted_bp="" if ins is None else ins,
                         deleted_bp="" if dele is None else dele,
                         insertions_ge50="" if nbig is None else nbig,
                         tech_from_sequence=seq,
                         flag=provenance_flag(ident, tc),
                         review=";".join(review)))
        if i % 50 == 0:
            print(f"    {i}/{len(accs)}", flush=True)
    if missing:
        # a genome that was not screened must not read as "ok" downstream
        sys.exit(f"no assembly for {len(missing)} accession(s): "
                 f"{' '.join(missing[:10])}")

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    import collections
    print(f"\n  {len(rows)} assemblies screened")
    for k, v in collections.Counter(r["flag"] for r in rows).most_common():
        print(f"    {k:<24}{v:>5}")
    print("  review notes (do not change the flag):")
    for k, v in collections.Counter(x for r in rows for x in r["review"].split(";")
                                    if x).most_common():
        print(f"    {k:<60}{v:>5}")
    print(f"\n  wrote {a.out}")


if __name__ == "__main__":
    main()
