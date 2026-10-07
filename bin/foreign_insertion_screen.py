#!/usr/bin/env python3
"""Find insertions that are foreign to the MTBC — vector, plasmid, engineered DNA.

Metadata cannot detect an engineered strain. A ~5 kb pJEB integration in one panel
genome was invisible in its RefSeq record: strain field "Erdman = ATCC 35801",
isolate "na", assembly method Flye, nothing amiss. It was found only by comparing
two assemblies of the same strain, where it was the single difference.

That comparison generalises without needing a matched pair. Take every insertion
of >= --min-len that a genome carries relative to the reference, and ask whether it
aligns in OTHER genomes of a small, diverse background:

  native   IS6110, prophage, RD-region sequence and duplications all have
           homologues elsewhere, because other genomes carry the same mobile
           elements or the source locus.
  foreign  vector backbone, selection markers and synthetic barcodes align
           nowhere, because no MTBC genome contains them.

WHAT IS EXTRACTED. Two kinds of insert, both >= --min-len:
  indel  insertions inside one alignment block (paftools call);
  gap    assembly sequence between two consecutive alignment blocks of
         >= --min-block that aligns to neither. An artifact that breaks the
         alignment is otherwise never seen: GCF_039770655's 1,144 bp
         (CCATT)n / poly-T artifact replaces 955 bp of PE_PGRS31, so its 40 kb
         window aligns as two blocks and paftools calls no indel (audit PGB-5e).

HOW AN INSERT IS SCORED (audit PGB-5a, b):
  * per genome, not summed. For each background genome other than the insert's
    own, the insert bases covered by that genome's alignments are merged; the
    best single genome gives best_homologue_bp. Summing over up to 50
    secondary hits let a short native piece that hits many genomes exceed its
    own length.
  * FOREIGN when best_homologue_bp < --max-foreign-frac (0.10) of the insert.
  * REVIEW_one_homologue when fewer than --min-homologue-genomes (2) other
    genomes carry >= that fraction. Two genomes carrying the same construct
    otherwise vouch for each other: GCF_044324775 (pJEB) and GCF_021535155
    (attB vector) both read "native" when both were in the background.
  * native otherwise.

THE RE-CHECK (ported from analysis/external_assemblies/bin/foreign_recheck.py).
minimap2 -x asm10's default frequency filter drops IS6110's minimizers (about
16 copies per background genome), so an IS6110 copy at a new site finds no
homologue: 952 false FOREIGN calls among the external assemblies, 27 among
CX333's. Every insert not called native is re-aligned (1) to the background
with the filter off (-f 1000000) and (2) to IS6110 (--is6110). Covered >= 50%
by IS6110 -> native_is6110; >= 50% by each of >= --min-homologue-genomes other
genomes -> native_recheck. `final_verdict` holds the outcome.

THE BACKGROUND (audit PGB-5c, d). It must be small and diverse: against the
whole 1.48 Gb panel minimap2 discards repetitive queries, and 2,427 of 4,936
inserts got no alignment (QC_PIPELINE.md 1.4; RUNBOOK.md's command used the
whole panel and called 1,167 inserts foreign in the external assemblies).
Pass --lineages and the script writes one genome per sublineage plus every
genome with no lineage call, PanSN-named, into the workdir. The genome kept
for a sublineage is the best by --rank-by, the same file and rule
snp_nonredundant.py uses to pick a cluster's representative (the user's
decisions D27 and D32, 2026-10-07): accession<TAB>score, higher is better,
a tie goes to the larger accession. It was the first accession in sorted
order, so the background's quality was an accident of numbering. A
sublineage member with no score stops the run rather than scoring 0. A --background FASTA is accepted only if every
name is PanSN (sample#hap#contig), because self-hits are excluded by sample,
and if it holds no more than --max-background genomes.

MINIMAP2 SETTINGS. Inserts against the background: -cx asm10 --secondary=yes
-N 100 -p 0.05, as QC_PIPELINE.md 1.4 documents. Secondary alignments are
REQUIRED: each insert's best hit is to the genome it came from, and with
--secondary=no the self-hit is the only alignment reported.

The workdir records its inputs (inputs.key); a workdir written for different
accessions, settings or background is refused rather than reused (PGB-5f).

THE VECTOR CHECK (the user's decision D29, 2026-10-07). The homologue rule
cannot see a construct that two genomes share; REVIEW_one_homologue catches
two, not three. So every insert, native or not, is also searched against
NCBI's UniVec_Core (--univec; bin/fetch_univec.sh) the way VecScreen does it:
blastn -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes
-soft_masking true -evalue 700 -searchsp 1750000000000, and a hit is STRONG at
score >= 30, or >= 24 within 25 bases of either end of the insert. Any strong
hit makes final_verdict VECTOR, whatever the homologue rule said;
univec_strong_bp is the insert bases those hits cover. Measured on the
inserts of CX333 (449) and of the 156 external assemblies (2,254): strong
hits on exactly the two known constructs, GCF_044324775's pJEB integration
(2,295 bp) and GCF_021535155's attB vector (2,107 bp), and on nothing else;
no moderate hit either.

The outgroup is not exempted here: M. canettii's divergent sequence reads
FOREIGN and is reviewed by hand (QC_PIPELINE.md 1.4).
"""
import argparse, collections, csv, hashlib, os, subprocess, sys

ASM10 = ["-cx", "asm10", "--secondary=yes", "-N", "100", "-p", "0.05"]
VECSCREEN = ["-task", "blastn", "-reward", "1", "-penalty", "-5", "-gapopen", "3",
             "-gapextend", "3", "-dust", "yes", "-soft_masking", "true",
             "-evalue", "700", "-searchsp", "1750000000000"]
BLAST_FMT = "6 qseqid sseqid qstart qend score qlen stitle"


def read_fasta(path):
    import gzip
    op = gzip.open if path.endswith(".gz") else open
    out, name, buf = {}, None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        out[name] = "".join(buf)
    return out


def alignment_gaps(paf_lines, min_len, min_block):
    """Unaligned query stretches >= min_len between consecutive primary
    alignment blocks of >= min_block on the same query contig.
    -> [(contig, qstart, qend, ref_pos)], ref_pos = the reference coordinate
    where the left block ends."""
    blocks = collections.defaultdict(list)
    for l in paf_lines:
        f = l.rstrip("\n").split("\t")
        if len(f) < 12 or "tp:A:P" not in l:
            continue
        qs, qe = int(f[2]), int(f[3])
        if qe - qs < min_block:
            continue
        ts, te = int(f[7]), int(f[8])
        blocks[f[0]].append((qs, qe, te if f[4] == "+" else ts))
    out = []
    for q, bs in blocks.items():
        bs.sort()
        end, rp = bs[0][1], bs[0][2]
        for qs, qe, r in bs[1:]:
            if qs - end >= min_len:
                out.append((q, end, qs, rp))
            if qe > end:
                end, rp = qe, r
    return out


def extract_inserts(acc, asm, ref, mm2, k8, paftools, min_len, min_block, out_fh,
                    paf_path=None):
    """Write every insert >= min_len as >acc|ref_pos|len|source. Returns the count."""
    res = subprocess.run([mm2, "-cx", "asm5", "--cs", "-t", "2", ref, asm],
                         capture_output=True, text=True)
    # a failed alignment (a missing --ref) read as "no inserts" and was then
    # marked done, so a rerun in the workdir reused the empty result
    if res.returncode != 0:
        sys.exit(f"FATAL: minimap2 failed on {acc} (exit {res.returncode}):\n"
                 f"{res.stderr[-2000:]}")
    paf = res.stdout
    if paf_path:
        with open(paf_path, "w") as fh:
            fh.write(paf)
    srt = subprocess.run(["sort", "-k6,6", "-k8,8n"], input=paf,
                         capture_output=True, text=True).stdout
    out = subprocess.run([k8, paftools, "call", "-L", str(min_block),
                          "-l", str(min_block), "-"],
                         input=srt, capture_output=True, text=True)
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
        out_fh.write(f">{acc}|{f[2]}|{len(a)-len(r)}|indel\n{a}\n")
    gaps = alignment_gaps(paf.splitlines(), min_len, min_block)
    if gaps:
        seqs = read_fasta(asm)
        for q, s, e, rp in gaps:
            n += 1
            out_fh.write(f">{acc}|{rp}|{e-s}|gap\n{seqs[q][s:e]}\n")
    return n


def coverage_by_genome(paf_lines, skip_self=True):
    """{query: {target sample: insert bases covered (merged)}}; the sample is
    the PanSN prefix of the target name, the insert's genome its first field."""
    iv = collections.defaultdict(lambda: collections.defaultdict(list))
    for l in paf_lines:
        f = l.split("\t")
        if len(f) < 10:
            continue
        q, sample = f[0], f[5].split("#")[0]
        if skip_self and sample == q.split("|")[0]:
            continue
        iv[q][sample].append((int(f[2]), int(f[3])))
    out = {}
    for q, by in iv.items():
        out[q] = {}
        for sample, xs in by.items():
            xs.sort(); tot, end = 0, -1
            for s, e in xs:
                if e > end:
                    tot += e - max(s, end); end = e
            out[q][sample] = tot
    return out


def score(length, per_genome, max_frac, min_genomes):
    """(best_bp, best_genome, n_genomes >= max_frac, verdict)"""
    if per_genome:
        g, best = max(per_genome.items(), key=lambda kv: (kv[1], kv[0]))
    else:
        g, best = "", 0
    n = sum(1 for v in per_genome.values() if v >= max_frac * length)
    if best < max_frac * length:
        v = "FOREIGN"
    elif n < min_genomes:
        v = "REVIEW_one_homologue"
    else:
        v = "native"
    return best, g, n, v


def recheck_verdict(length, verdict, per_genome, is_bp, min_genomes):
    if verdict == "native":
        return "native"
    if is_bp >= 0.5 * length:
        return "native_is6110"
    if sum(1 for v in per_genome.values() if v >= 0.5 * length) >= min_genomes:
        return "native_recheck"
    return verdict


def vecscreen_strength(score, qs, qe, qlen):
    """VecScreen's match categories: terminal = within 25 bases of an end."""
    terminal = min(qs, qe) <= 25 or max(qs, qe) >= qlen - 24
    for name, t, i in (("strong", 24, 30), ("moderate", 19, 25), ("weak", 16, 23)):
        if score >= (t if terminal else i):
            return name
    return ""


def univec_hits(blast_lines):
    """{insert: (strong bp covered, best strong hit title)} from BLAST_FMT
    output; inserts with no strong hit are absent."""
    iv, best = collections.defaultdict(list), {}
    for l in blast_lines:
        f = l.rstrip("\n").split("\t")
        if len(f) < 7:
            continue
        q, qs, qe, sc, ql = f[0], int(f[2]), int(f[3]), int(float(f[4])), int(f[5])
        if vecscreen_strength(sc, qs, qe, ql) != "strong":
            continue
        iv[q].append((min(qs, qe), max(qs, qe)))
        if q not in best or sc > best[q][0]:
            best[q] = (sc, f[6])
    out = {}
    for q, xs in iv.items():
        xs.sort(); tot, end = 0, 0
        for s, e in xs:
            if e > end:
                tot += e - max(s, end + 1) + 1; end = e
        out[q] = (tot, best[q][1])
    return out


def apply_univec(rows, hits):
    """A strong UniVec hit makes an insert VECTOR, native or not."""
    for q, (bp, hit) in hits.items():
        if q in rows:
            rows[q].update(univec_strong_bp=bp, univec_hit=hit, final_verdict="VECTOR")


def load_rank(path):
    """--rank-by, read as snp_nonredundant.py reads it: accession<TAB>score;
    lines whose score is not a number (a header) are skipped."""
    score = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            try:
                score[f[0]] = float(f[1])
            except ValueError:
                pass
    return score


def background_from_lineages(accs, lineages, score):
    """one accession per sublineage, the best by `score` (D27's rule: higher
    wins, a tie goes to the larger accession), plus every accession with no
    call"""
    lin = {r["accession"]: (r.get("strain") or "") for r in
           csv.DictReader(open(lineages), delimiter="\t")}
    groups, out = collections.defaultdict(list), []
    for acc in sorted(accs):
        s = lin.get(acc, "").split(":")[0]
        if not s:
            out.append(acc)
        else:
            groups[s].append(acc)
    unscored = sorted(a for m in groups.values() if len(m) > 1
                      for a in m if a not in score)
    if unscored:
        sys.exit(f"FATAL: --rank-by has no score for {len(unscored)} genome(s) "
                 f"competing for a background slot: {' '.join(unscored[:10])}")
    pick = [max(m, key=lambda a: (score.get(a, 0), a)) for m in groups.values()]
    return sorted(out + pick)


def check_background(path, max_genomes):
    samples = set()
    import gzip
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name.count("#") < 2:
                    sys.exit(f"{path}: '{name}' is not PanSN (sample#hap#contig); "
                             "self-hits could not be excluded")
                samples.add(name.split("#")[0])
    if len(samples) > max_genomes:
        sys.exit(f"{path} holds {len(samples)} genomes (> --max-background "
                 f"{max_genomes}). The background must be small and diverse: pass "
                 "--lineages to build one genome per sublineage")
    return samples


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accessions", required=True)
    ap.add_argument("--assembly-dir", required=True)
    ap.add_argument("--assembly-suffix", default=".dnaA_rotated.fasta")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--lineages",
                    help="lineages.all.tsv: build the one-per-sublineage "
                         "background from the candidates")
    ap.add_argument("--rank-by",
                    help="with --lineages, required: accession<TAB>score, higher "
                         "is a better background genome; the file "
                         "snp_nonredundant.py --rank-by takes (D27, D32)")
    ap.add_argument("--background", "--panel-fasta", dest="background",
                    help="a PanSN-named background FASTA instead of --lineages")
    ap.add_argument("--max-background", type=int, default=120)
    ap.add_argument("--is6110", required=True,
                    help="canonical IS6110 FASTA, for the re-check")
    ap.add_argument("--univec", required=True,
                    help="NCBI UniVec_Core FASTA (bin/fetch_univec.sh), for the "
                         "vector check every insert gets (D29)")
    ap.add_argument("--blastn", default=os.environ.get("MTB_BLASTN", "blastn"))
    ap.add_argument("--makeblastdb",
                    default=os.environ.get("MTB_MAKEBLASTDB", "makeblastdb"))
    ap.add_argument("--min-len", type=int, default=1000)
    ap.add_argument("--min-block", type=int, default=5000)
    ap.add_argument("--max-foreign-frac", type=float, default=0.10,
                    help="an insert with less than this fraction of its length "
                         "aligned to any single other genome is called foreign")
    ap.add_argument("--min-homologue-genomes", type=int, default=2,
                    help="other genomes that must carry the insert for 'native'")
    ap.add_argument("--minimap2", required=True)
    ap.add_argument("--k8", required=True)
    ap.add_argument("--paftools", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if bool(a.lineages) == bool(a.background):
        ap.error("give exactly one of --lineages or --background")
    if a.lineages and not a.rank_by:
        ap.error("--lineages needs --rank-by: the background keeps the best "
                 "genome of each sublineage (D32)")

    os.makedirs(a.workdir, exist_ok=True)
    accs = [l.strip() for l in open(a.accessions) if l.strip()]
    missing = [x for x in accs
               if not os.path.exists(os.path.join(a.assembly_dir, x + a.assembly_suffix))]
    if missing:
        sys.exit(f"no assembly for {len(missing)} accession(s): {' '.join(missing[:10])}")

    if a.lineages:
        bg_ids = background_from_lineages(accs, a.lineages, load_rank(a.rank_by))
        bg = os.path.join(a.workdir, "background.fasta")
        with open(bg + ".ids", "w") as fh:
            fh.write("\n".join(bg_ids) + "\n")
        with open(bg, "w") as fh:
            for acc in bg_ids:
                for c, s in read_fasta(os.path.join(
                        a.assembly_dir, acc + a.assembly_suffix)).items():
                    fh.write(f">{acc}#1#{c}\n{s}\n")
    else:
        bg = a.background
    samples = check_background(bg, a.max_background)
    print(f"  background: {len(samples)} genomes ({bg})")

    with open(a.univec, "rb") as fh:
        uv_sha = hashlib.sha256(fh.read()).hexdigest()
    key = hashlib.sha256("\n".join(
        [a.ref, a.assembly_dir, a.assembly_suffix, str(a.min_len), str(a.min_block),
         uv_sha]
        + accs + sorted(samples)).encode()).hexdigest()
    kpath = os.path.join(a.workdir, "inputs.key")
    ins_fa = os.path.join(a.workdir, "inserts.fa")
    if os.path.exists(kpath) and open(kpath).read().strip() != key:
        sys.exit(f"{a.workdir} was written for other inputs or settings; use a "
                 "new --workdir")
    with open(kpath, "w") as fh:
        fh.write(key + "\n")

    done = os.path.join(a.workdir, "inserts.done")
    if not os.path.exists(done):
        print(f"  extracting insertions >= {a.min_len} bp from {len(accs)} genomes ...")
        with open(ins_fa, "w") as fh:
            for i, acc in enumerate(accs, 1):
                asm = os.path.join(a.assembly_dir, acc + a.assembly_suffix)
                extract_inserts(acc, asm, a.ref, a.minimap2, a.k8, a.paftools,
                                a.min_len, a.min_block, fh)
                if i % 50 == 0:
                    print(f"    {i}/{len(accs)}", flush=True)
        open(done, "w").close()
    seqs = read_fasta(ins_fa)
    print(f"  {len(seqs)} insertions >= {a.min_len} bp")
    if not seqs:
        sys.exit("no insertions extracted")

    paf = os.path.join(a.workdir, "inserts_vs_background.paf")
    pdone = paf + ".done"
    if not os.path.exists(pdone):
        print("  aligning inserts against the background ...")
        with open(paf, "w") as fh:
            subprocess.run([a.minimap2] + ASM10 + ["-t", str(a.threads), bg, ins_fa],
                           stdout=fh, stderr=subprocess.DEVNULL, check=True)
        open(pdone, "w").close()
    cov = coverage_by_genome(open(paf))

    rows = {}
    for q, s in seqs.items():
        acc, pos, ln, src = q.split("|")
        ln = int(ln)
        best, g, n, v = score(ln, cov.get(q, {}), a.max_foreign_frac,
                              a.min_homologue_genomes)
        rows[q] = dict(accession=acc, ref_pos=int(pos), insert_len=ln, source=src,
                       best_homologue_bp=best, best_homologue_genome=g,
                       n_homologue_genomes=n, frac=round(best / max(ln, 1), 4),
                       verdict=v, recheck_best_bp="", recheck_n_genomes="",
                       is6110_bp="", univec_strong_bp=0, univec_hit="",
                       final_verdict=v)

    todo = [q for q, r in rows.items() if r["verdict"] != "native"]
    if todo:
        sub = os.path.join(a.workdir, "recheck.fa")
        with open(sub, "w") as fh:
            for q in todo:
                fh.write(f">{q}\n{seqs[q]}\n")
        rb = subprocess.run([a.minimap2] + ASM10 + ["-f", "1000000", "-t", str(a.threads),
                                                    bg, sub],
                            capture_output=True, text=True, check=True).stdout.splitlines()
        ri = subprocess.run([a.minimap2, "-c", "-N", "10", "-p", "0.1", a.is6110, sub],
                            capture_output=True, text=True, check=True).stdout.splitlines()
        cb = coverage_by_genome(rb)
        ci = coverage_by_genome(ri, skip_self=False)
        for q in todo:
            r = rows[q]
            pg = cb.get(q, {})
            isb = sum(ci.get(q, {}).values())
            r["recheck_best_bp"] = max(pg.values()) if pg else 0
            r["recheck_n_genomes"] = sum(1 for v in pg.values() if v >= 0.5 * r["insert_len"])
            r["is6110_bp"] = isb
            r["final_verdict"] = recheck_verdict(r["insert_len"], r["verdict"], pg, isb,
                                                 a.min_homologue_genomes)

    # every insert against UniVec, native ones included (D29)
    udb = os.path.join(a.workdir, "univec")
    subprocess.run([a.makeblastdb, "-in", a.univec, "-dbtype", "nucl", "-out", udb],
                   capture_output=True, text=True, check=True)
    ub = subprocess.run([a.blastn] + VECSCREEN + ["-db", udb, "-query", ins_fa,
                                                  "-num_threads", str(a.threads),
                                                  "-outfmt", BLAST_FMT],
                        capture_output=True, text=True, check=True).stdout.splitlines()
    apply_univec(rows, univec_hits(ub))

    out = sorted(rows.values(), key=lambda r: (r["final_verdict"].startswith("native"),
                                               -r["insert_len"]))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(out)

    c = collections.Counter(r["final_verdict"] for r in out)
    print(f"\n  first pass: {dict(collections.Counter(r['verdict'] for r in out))}")
    print(f"  after re-check: {dict(c)}")
    fo = [r for r in out if not r["final_verdict"].startswith("native")]
    bg_n = collections.Counter(r["accession"] for r in fo)
    print(f"  genomes with an insert to review: {len(bg_n)}\n")
    print(f"  {'accession':<16}{'n':>4}{'total bp':>10}{'largest':>9}  verdicts / positions")
    for acc, n in bg_n.most_common():
        sub = sorted((r for r in fo if r["accession"] == acc), key=lambda r: -r["insert_len"])
        tot = sum(r["insert_len"] for r in sub)
        pos = ", ".join(f"{r['ref_pos']}:{r['final_verdict']}" for r in sub[:4])
        print(f"  {acc:<16}{n:>4}{tot:>10,}{sub[0]['insert_len']:>9,}  {pos}")
    print(f"\n  wrote {a.out}")


if __name__ == "__main__":
    main()
