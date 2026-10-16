#!/usr/bin/env python3
"""One merged, multi-sample VCF per cohort: small variants, structural
variants and IS6110 insertions, in one coordinate frame, with genotypes that
distinguish no data from reference from a deleted region.

WHY THIS EXISTS. Every product of this pipeline up to now is either per-sample
or a TSV. Nothing a standard consumer reads carries the cohort, which was the
objection that started this work.

THE THREE THINGS IT HAS TO GET RIGHT.

1. Every class in one file. Small variants and their placements come from the
   P5 matrix, structural variants from the SV matrix, IS6110 insertions from
   the p1iv cohort key table. Accessory regions are declared in section 4 of
   refbias/ACCESSORY_IN_VCF_PLAN.md and are NOT emitted yet, because the
   presence rule is calibrated on five isolates and 21 positive pairs; a
   genotype written from that would be a number invented at the last step. The
   header says so rather than leaving a reader to assume completeness.

2. The node frame is carried, not dropped. 7,454 of scale100's 61,049 matrix
   keys and 388 of its 1,125 IS6110 sites have no H37Rv position. Dropping them
   would discard exactly the sequence the matched-reference arm exists to
   reach. They are emitted on their own contigs, `node_<id>`, declared in the
   header, with the graph node in INFO. This is the mechanism the accessory
   work arrived at independently and the one rGFA formalises with SN/SO/SR.

3. Missing, reference and deleted are three different genotypes.
       GT=0   the sample matches the reference here, with data to say so
       GT=1   the sample carries the alternate allele
       GT=2   the position is DELETED in this sample -- emitted against the
              `*` allele, which is what VCF 4.2 provides for a base covered by
              an upstream deletion, so a consumer that knows the spec reads it
              correctly without reading our documentation
       GT=.   no data: no coverage, or a call that did not survive filtering
   The verbatim pipeline state is also kept in FORMAT/ST, so nothing is lost to
   the mapping. A naive merge of per-sample files cannot make this distinction
   at all, because an absent record is silent about which of the three it was.
"""
import argparse, bisect, collections, csv, glob, gzip, json, os, subprocess, sys, datetime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harmonize_records as _harm

STATE_GT = {"ALT": "1", "REF": "0", "ABSENT": "2", "NOCALL": "."}
META = {"key", "frame", "region", "kind", "h37rv_pos", "node", "node_offset",
        "acc_locus", "canonical_ref", "canonical_alt", "panel_af",
        "h37rv_minor", "n_alt", "n_ref", "n_absent", "n_nocall",
        "svtype", "svlen", "size_band", "n_both_callers", "src", "max_qual",
        "qual_band", "stage11_ppv", "qual_from", "component"}

def read_fasta_seq(path):
    """The first sequence of a (optionally gzipped) FASTA, upper-cased."""
    op = gzip.open if path.endswith(".gz") else open
    name, buf = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].split()[0]
            else:
                buf.append(line.strip())
    return name, "".join(buf).upper()


def read_matrix(path):
    if not os.path.exists(path):
        return [], []
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        cols = rdr.fieldnames or []
        samples = [c for c in cols if c not in META]
        return list(rdr), samples


def iter_matrix(path):
    """Samples, then row dicts one at a time.

    read_matrix returns list(rdr), which at 997 samples is 173,660 dicts of
    1,005 entries -- about 15 GB before a single record is written. The merge
    reads each row once and never looks back, so it does not need them all at
    once. Kept beside read_matrix, which the small SV matrix still uses.
    """
    if not os.path.exists(path):
        return [], iter(())
    fh = open(path)
    rdr = csv.DictReader(fh, delimiter="\t")
    cols = rdr.fieldnames or []
    samples = [c for c in cols if c not in META]

    def gen():
        try:
            for r in rdr:
                yield r
        finally:
            fh.close()
    return samples, gen()


def read_sparse(state_dir, keys_path, refmap_path):
    """The same rows read_matrix returns, built from the sparse state files.

    The dense matrix.tsv is keys x samples of text: 783 MB at 997 isolates and
    26.5 GB projected at 10,000, of which about 98% is the words REF and NOCALL.
    Every field this merge reads is either a property of the key -- already in
    keys.tsv -- or a per-sample state, so the dense file never has to be built.
    """
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from p5_states_io import NAME, load_states, samples_in
    keys = list(csv.DictReader(open(keys_path, newline=""), delimiter="\t"))
    order = [r["sample"] for r in
             csv.DictReader(open(refmap_path, newline=""), delimiter="\t")]
    samples = samples_in(state_dir, order)
    M, EV, codes = load_states(state_dir, keys, samples, with_evidence=True)
    # keys.tsv carries two columns matrix.tsv silently dropped, and both are
    # withheld here so this path reproduces the dense one exactly rather than
    # changing published output as a side effect of a storage change:
    #
    #   node_offset   every node-frame record has therefore been emitted at
    #                 POS=1 of its contig regardless of its real offset. The
    #                 contigs are also DECLARED with length 1, so simply using
    #                 the offset would produce invalid VCF -- the fix is to
    #                 declare each node contig at its true length and then use
    #                 it, which changes every merged VCF and is its own change.
    #   canonical_alt the merge currently recovers the ALT by parsing it out of
    #                 the key string, with a comment saying no such column
    #                 exists. It exists in keys.tsv.
    DROP = ("node_offset", "canonical_alt")
    print(f"  streaming {len(keys):,} keys from {len(samples)} sparse state "
          f"files without a dense matrix")

    def gen():
        for i, k in enumerate(keys):
            col = M[i]
            if not (col == 0).any():      # no ALT carrier: not a site
                continue
            r = {kk: vv for kk, vv in k.items() if kk not in DROP}
            for j, sm in enumerate(samples):
                r[sm] = NAME[col[j]]
            _attach_evidence(r, EV[i], samples, codes.names)
            yield r
    return gen(), samples


def _attach_evidence(r, evrow, samples, names):
    """r["_ev"] = {sample: evidence code} for the cells that have one: the
    reason behind an ABSENT or NOCALL (HANDOFF ABS-1 item 4), written to
    FORMAT/EV. Rows of a cohort finished before the evidence array existed
    get none."""
    if evrow is None:
        return
    import numpy as np
    nz = np.flatnonzero(np.asarray(evrow))
    if len(nz):
        r["_ev"] = {samples[j]: names[evrow[j]] for j in nz}

def read_array(array_dir, keys_path, want=None):
    """Rows as read_sparse returns them, from the memory-mapped states array.

    `want(key_row)` is asked BEFORE the per-sample cells are built, so a shard
    pays only for its own keys: at 10,000 isolates building every row's cells
    is 540,000 x 10,000 strings, which is the cost sharding exists to divide.
    """
    import numpy as np
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from p5_states_io import NAME, open_array, open_evidence
    keys = list(csv.DictReader(open(keys_path, newline=""), delimiter="\t"))
    M, samples = open_array(array_dir, keys)
    EV, ev_names = open_evidence(array_dir, M.shape)
    DROP = ("node_offset", "canonical_alt")     # as read_sparse, see there

    def gen():
        for i, k in enumerate(keys):
            if want is not None and not want[0](k):
                continue
            col = np.asarray(M[i])
            if not (col == 0).any():              # no ALT carrier: not a site
                continue
            r = {kk: vv for kk, vv in k.items() if kk not in DROP}
            for j, sm in enumerate(samples):
                r[sm] = NAME[col[j]]
            if EV is not None:
                _attach_evidence(r, EV[i], samples, ev_names)
            yield r
    return gen(), samples, keys


def _header_body(w):
    """The fixed ALT/INFO/FORMAT header lines, verbatim from the single-file writer."""
    w('##ALT=<ID=*,Description="Position deleted in this sample, per VCF 4.2 '
      'spanning-deletion semantics. GT=2 means the region is ABSENT, which is '
      'a different statement from GT=. meaning no data">\n')
    w('##ALT=<ID=INS:ME:IS6110,Description="IS6110 mobile element insertion">\n')
    w('##ALT=<ID=UNKNOWN,Description="Alternate allele not represented in the '
      'cohort matrix; the genotype and FORMAT/ST still carry the call">\n')
    w('##INFO=<ID=CLASS,Number=1,Type=String,Description="Variant class: '
      'small, sv, is6110 or accessory_presence. accessory_presence records '
      'state whether each sample carries a catalogued accessory locus '
      '(level 1); variation inside those loci is CLASS=small on node '
      'contigs">\n')
    w('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of '
      'a symbolic structural variant: the last deleted base">\n')
    w('##INFO=<ID=FRAME,Number=1,Type=String,Description="h37rv or node. A '
      'node-frame record sits on sequence the H37Rv path does not carry and is '
      'emitted on its own node_<id> contig rather than dropped">\n')
    w('##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele at '
      'SNP sites, reconstructed by Fitch parsimony over the panel tree '
      '(bin/panel_ancestral.py): the state of the ancestral node, the most '
      'recent common ancestor of every non-outgroup panel genome, a tie there '
      'resolved by the downward pass over all outgroups together. A property '
      'of the panel rather than of this cohort. `.` where that node is still '
      'tied; see AA_FLAG. Indels, MNPs and SVs carry AA_STATE instead">\n')
    w('##INFO=<ID=AA_FLAG,Number=.,Type=String,Description="How far to '
      'trust the ancestral state (fix AA-1). TIED:<bases> where parsimony '
      'does not resolve it (never a coin toss); NODATA where the site had '
      'none; OUTGROUP where the ingroup alone was tied and the outgroups '
      'decided it; HOMOPLASTIC where the panel tree needs at least the build\'s '
      'threshold of changes (bin/panel_ancestral.py --homoplastic-min, '
      'default 2). HOMOPLASTIC combines with the others">\n')
    w('##INFO=<ID=AA_CHANGES,Number=1,Type=Integer,Description="Fitch '
      'change count of this key on the panel tree, over the graph genomes '
      '(outgroup kit genomes left out)">\n')
    w('##INFO=<ID=AA_STATE,Number=1,Type=String,Description="REF or ALT: '
      'the allele the panel reconstruction puts at the ancestral node, for a '
      'record with no AA tag (indel, MNP, SV) whose key the build\'s '
      'panel_polarity.tsv resolves. The event writer takes an ALT-ancestral '
      'state only where the panel carries the ALT at >= its --min-panel-af">\n')
    w('##INFO=<ID=AA_SOURCE,Number=1,Type=String,Description="On every '
      'record. reconstructed: the ancestral state comes from the panel '
      'reconstruction (AA or AA_STATE). assumed: no reconstructed state, so '
      'downstream polarity takes the ALT (presence) as derived by '
      'assumption">\n')
    w('##INFO=<ID=AA_INVERTED,Number=0,Type=Flag,Description="The REFERENCE '
      'carries the DERIVED allele at this site, so GT=0 is derived and GT=1 '
      'is ancestral, which is why AA is annotated rather than REF being '
      're-polarised -- VCF requires REF to match the reference base">\n')
    w('##INFO=<ID=IS6110PROX,Number=0,Type=Flag,Description="The interval '
      'lies within 500 bp of an IS6110 landmark, where depth counts '
      'mismapped copies of the element unless it is filtered">\n')
    w('##INFO=<ID=NODE,Number=1,Type=String,Description="Graph node, for '
      'node-frame records">\n')
    w('##INFO=<ID=ACCIS6110,Number=0,Type=Flag,Description="An accessory '
      'presence locus whose sequence is an IS6110 copy: the same event as the '
      'IS6110 arm\'s insertion record at the site, so count one of the two">\n')
    w('##INFO=<ID=UNCATALOGUED,Number=0,Type=Flag,Description="A caller '
      'deletion no interval of the catalogue covers, so it is written '
      'presence-only from the caller matrix instead of being dropped">\n')
    for t, d in (("KIND", "small-variant kind from P4"),
                 ("REGION", "core, pe_ppe, masked, off_path_near or off_path_accessory"),
                 ("SVTYPE", "structural variant type"),
                 ("SVLEN", "structural variant length"),
                 ("COMPONENT", "called or inherited; inherited rows carry no caller QUAL"),
                 ("QUALBAND", "caller quality band, or composed"),
                 ("MEINFO", "mobile element info: name,start,end,polarity"),
                 ("EVIDENCE", "IS6110 support: two_sided, one_sided, or mixed where carriers disagree"),
                 ("SVSOURCE", "graph or caller: where the interval came from"),
                 ("SVTIER", "interval support tier: A_multi_assembly, B_single_assembly or C_caller_only"),
                 ("PANELCARRIERS", "panel genomes carrying the deletion, from the graph"),
                 ("SITECLASS", "IS6110 relationship to the reference: ref_lacking, ref_shared, or mixed where carriers disagree"),
                 ("ACCLOCUS", "accessory locus id, for level-1 presence records"),
                 ("ACCLEN", "length of the accessory locus representative allele"),
                 ("ACCKLASS", "accessory locus class from the catalogue"),
                 ("ACCROUTE", "instrument the locus is routed to")):
        num = "1"
        typ = ("Integer" if t in ("SVLEN", "PANELCARRIERS", "ACCLEN")
               else "String")
        if t == "MEINFO":
            num = "4"
        w(f'##INFO=<ID={t},Number={num},Type={typ},Description="{d}">\n')
    w('##FORMAT=<ID=GT,Number=1,Type=String,Description="Haploid genotype. '
      '0 reference with data, 1 alternate, 2 the * allele meaning this region '
      'is deleted in this sample, . no data">\n')
    w('##FORMAT=<ID=ST,Number=1,Type=String,Description="Pipeline state '
      'verbatim: ALT, REF, ABSENT or NOCALL; DUP where this arm\'s call was '
      'absorbed into the record of the arm that owns the event '
      '(bin/harmonize_records.py). Kept so the GT mapping loses '
      'nothing">\n')
    w('##FORMAT=<ID=EV,Number=1,Type=String,Description="Evidence behind the '
      'cell, where the pipeline records one: why an ABSENT or NOCALL is what it '
      'is (e.g. masked, junction_joined, no_projection, frames_conflict, '
      'tandem), how an inherited SV was confirmed (two_frames, h37rv_only, '
      'matched_only), how an IS6110 site off the matched reference\'s path was '
      'read in the sample\'s H37Rv alignment (absent_checked_alt, '
      'absent_checked_ref, absent_masked, absent_low_depth), and '
      'node_not_on_r / inherited_unchecked for calls that '
      'rest only on the matched reference. Trailing field: absent where none '
      'is recorded">\n')
    for line in _harm.HEADER:
        w(line)


# ---- the VCF header, shared by the single-file and the assembled output ----
def write_header(w, a, nodes_used, order):
    today = datetime.date.today().strftime("%Y%m%d")
    w("##fileformat=VCFv4.2\n")
    w(f"##fileDate={today}\n")
    w("##source=bin/merge_cohort_vcf.py\n")
    w(f"##cohort={a.cohort_name}\n")
    w(f"##graphBuild={a.build_id}\n")
    w(f"##reference={a.h37rv_contig}\n")
    w(f"##contig=<ID={a.h37rv_contig},length={a.h37rv_length}>\n")
    for nid in sorted(nodes_used, key=lambda x: int(x) if x.isdigit() else 0):
        w(f"##contig=<ID=node_{nid},length={max(nodes_used[nid],1)}>\n")
    _header_body(w)
    w("#" + "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                       "FILTER", "INFO", "FORMAT"] + order) + "\n")


_EV_OK = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
             "0123456789_.|/+-")


def ev_value(text):
    """An evidence string as a FORMAT value: no tab, colon, comma, semicolon,
    space or `=`, which would break the cell or a reader that splits it."""
    return "".join(c if c in _EV_OK else "_" for c in str(text)).strip("_")


def vcf_cell(state, ev=""):
    """One sample's GT:ST[:EV] cell. EV is a trailing FORMAT field, written
    only where there is evidence (VCF allows trailing fields to be dropped),
    so a record with none keeps the plain GT:ST form."""
    c = f'{STATE_GT.get(state, ".")}:{state or "NA"}'
    ev = ev_value(ev) if ev else ""
    return c + ":" + ev if ev else c


def record_line(rec):
    """(VCF line, whether * was added) for one record tuple."""
    chrom, pos, ident, ref, alts, info, cells = rec
    # FORMAT is GT:ST, or GT:ST:EV when any cell carries evidence
    fmt = "GT:ST:EV" if any(c.count(":") >= 2 for c in cells.split("\t")) else "GT:ST"
    # GT=2 is the * allele (ABSENT). Every block maps ABSENT to 2, but only the
    # small-variant block used to add *, so the SV, interval, accessory and
    # IS6110 records wrote GT=2 against one ALT -- 976 records in scale200, and
    # `bcftools +fill-tags` aborts on the first. Decide it here, once, from the
    # cells actually written.
    added = False
    if "*" not in alts and ("\t" + cells).find("\t2:") != -1:
        alts = list(alts) + ["*"]
        added = True
    return ("\t".join([chrom, str(pos), ident, ref, ",".join(alts), ".",
                       "PASS", ";".join(info), fmt]) + "\t" + cells + "\n",
            added)


def assemble(a):
    """Header plus every shard's parts, streamed, into one bgzipped VCF.

    H37Rv parts first in shard order, then node parts in shard order: shards
    own contiguous H37Rv position ranges and contiguous (lexicographic) ranges
    of node contig names, so this concatenation is the sorted order a single
    process writes -- the assembled file is the unsharded one, byte for byte
    once decompressed. Memory is one read buffer, whatever the cohort size.
    """
    metas = []
    for i in range(a.n_shards):
        p = f"{a.part_prefix}.{i}.meta.json"
        if not os.path.exists(p):
            sys.exit(f"FATAL: shard {i} has no {p}; did every shard finish?")
        metas.append(json.load(open(p)))
    m0 = metas[0]
    for i, m in enumerate(metas):
        for k in ("order", "cohort_name", "build_id", "n_shards", "h37rv_contig"):
            if m[k] != m0[k]:
                sys.exit(f"FATAL: shard {i} disagrees with shard 0 on {k}")
        if m["shard"] != i:
            sys.exit(f"FATAL: {a.part_prefix}.{i}.meta.json is shard {m['shard']}")
    if m0["n_shards"] != a.n_shards:
        sys.exit(f"FATAL: shards were written for {m0['n_shards']}, not {a.n_shards}")
    nodes_used = {}
    for m in metas:
        for nid, ln in m["nodes_used"].items():
            nodes_used[nid] = max(nodes_used.get(nid, 0), ln)
    tmp = a.out + ".tmp"
    with open(tmp, "wb") as outfh:
        bz = subprocess.Popen([a.bgzip, "-c"], stdin=subprocess.PIPE, stdout=outfh)
        buf = []
        write_header(buf.append, a, nodes_used, m0["order"])
        bz.stdin.write("".join(buf).encode())
        for part in ("h37rv", "node"):
            for i in range(a.n_shards):
                src = f"{a.part_prefix}.{i}.{part}.gz"
                dz = subprocess.Popen([a.bgzip, "-dc", src], stdout=subprocess.PIPE)
                while True:
                    chunk = dz.stdout.read(1 << 22)
                    if not chunk:
                        break
                    bz.stdin.write(chunk)
                if dz.wait() != 0:
                    sys.exit(f"FATAL: could not decompress {src}")
        bz.stdin.close()
        if bz.wait() != 0:
            sys.exit("FATAL: bgzip failed while assembling")
    os.replace(tmp, a.out)
    subprocess.run([a.tabix, "-f", "-p", "vcf", a.out], check=True)
    tot = collections.Counter()
    for m in metas:
        tot.update(m["counts"])
    n_rec = tot["records"]
    print(f"  assembled {a.n_shards} shards: {len(m0['order'])} samples, "
          f"{n_rec:,} records")
    print(f"    small variants {tot['small']:,}   structural "
          f"{tot['sv'] + tot['iv']:,} ({tot['iv']:,} catalogued intervals, "
          f"{tot['sv']:,} caller rows)   IS6110 {tot['is6110']:,}   "
          f"accessory {tot['l1']:,}")
    print(f"    contigs: {a.h37rv_contig} plus {len(nodes_used):,} node_* contigs")
    print(f"    ancestral allele on {tot['aa']:,} records, of which "
          f"{tot['aa_inverted']:,} have the reference carrying the derived allele")
    if tot["star_added"]:
        print(f"    {tot['star_added']:,} symbolic records carry * for ABSENT cells")
    if tot["is6110_unplaced"]:
        print(f"    {tot['is6110_unplaced']:,} IS6110 sites EXCLUDED: no position "
              f"in either frame, so no VCF record is possible")
    print(f"  -> {a.out}")
    return 0


def load_panel_polarity(path):
    """{(pos, REF, ALT): (state, changes, flag)} from the build's
    panel_polarity.tsv (bin/panel_ancestral.py), H37Rv keys only; {} for no
    path. Rows whose `ancestral` is not REF or ALT are skipped."""
    if not path:
        return {}
    if not os.path.exists(path):
        sys.exit(f"FATAL: no panel polarity table at {path}")
    out = {}
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        need = {"pos", "ref", "alt", "ancestral"}
        if not need <= set(rd.fieldnames or []):
            sys.exit(f"FATAL: {path} lacks columns "
                     f"{sorted(need - set(rd.fieldnames or []))}; rebuild it "
                     f"with bin/panel_ancestral.py (P0 step ancestral)")
        for r in rd:
            if r["ancestral"] not in ("REF", "ALT"):
                continue
            out[(int(r["pos"]), r["ref"].upper(), r["alt"].upper())] = (
                r["ancestral"], r.get("changes", ""), r.get("flag", ""))
    return out


def annotate_aa_source(recs, contig, pol):
    """AA_SOURCE on every record (fix AA-1, point 3), so a polarity that is
    an assumption is visible in the VCF rather than only in the event
    writer's summary.

      reconstructed  the record has a resolved AA (SNPs: ancestral.tsv), or
                     its key is a resolved row of panel_polarity.tsv: an
                     indel, MNP or sequence-resolved SV by its own (pos, REF,
                     ALT); a symbolic SV by its SMALLALLELE key(s), all of
                     which must resolve and agree. Such a record also gets
                     AA_STATE (REF or ALT), AA_CHANGES and the row's AA_FLAG.
      assumed        anything else: tied or no-data SNPs, IS6110 and
                     accessory records, catalogued intervals and caller rows
                     without a small allele, node-frame records.

    Returns (records, Counter of AA_SOURCE values)."""
    n = collections.Counter()
    out = []
    for rec in recs:
        chrom, pos, ident, ref, alts, info = rec[:6]
        info = list(info)
        tags = {x.split("=", 1)[0]: x.split("=", 1)[1] for x in info if "=" in x}
        src = "assumed"
        if tags.get("AA", ".") not in (".", ""):
            src = "reconstructed"
        elif pol and chrom == contig and "AA" not in tags:
            real = [x for x in alts if x != "*"]
            keys = []
            if len(real) == 1 and not real[0].startswith("<"):
                keys = [(pos, ref.upper(), real[0].upper())]
            elif tags.get("SMALLALLELE"):
                for k in tags["SMALLALLELE"].split(","):
                    f = k.split(":")
                    if len(f) == 3 and f[0] == "h37rv" and ">" in f[2]:
                        r_, x_ = f[2].split(">", 1)
                        keys.append((int(f[1]), r_.upper(), x_.upper()))
                    else:
                        keys = []          # a key this cannot read: no claim
                        break
            hits = [pol.get(k) for k in keys]
            if hits and all(hits) and len({h[0] for h in hits}) == 1:
                src = "reconstructed"
                info.append(f"AA_STATE={hits[0][0]}")
                if len(hits) == 1:
                    if hits[0][1] != "":
                        info.append(f"AA_CHANGES={hits[0][1]}")
                    if hits[0][2]:
                        info.append(f"AA_FLAG={hits[0][2]}")
        info.append(f"AA_SOURCE={src}")
        n[src] += 1
        out.append(rec[:5] + (info,) + tuple(rec[6:]))
    return out, n


IS_COPY_LEN = (1300, 1420)      # harmonize_records.IS_LEN
IS_COPY_K, IS_COPY_FRAC = 20, 0.9


def is6110_copy_loci(catalogue, h37, gff, cat):
    """Ids of the catalogue loci in `cat` whose sequence is an IS6110 copy:
    IS_COPY_LEN long, with >= IS_COPY_FRAC of its 20-mers in H37Rv's first
    annotated IS6110 copy (either strand). Read from the catalogue's FASTA
    beside its TSV. Without the H37Rv FASTA, the IS6110 GFF or that FASTA,
    nothing is flagged and that is said."""
    fa = catalogue[:-4] + ".fasta" if catalogue.endswith(".tsv") else ""
    if not (h37 and gff and fa and os.path.exists(fa)):
        print("  NOTE level 1: no H37Rv FASTA, IS6110 GFF or catalogue FASTA; "
              "IS6110-copy loci are not flagged ACCIS6110")
        return set()
    copies = _harm.read_is6110_copies(gff)
    if not copies:
        sys.exit(f"FATAL: no IS6110 copy in {gff}")
    el = h37[copies[0][0] - 1:copies[0][1]]
    rc = el.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
    km = {s[i:i + IS_COPY_K] for s in (el, rc)
          for i in range(len(s) - IS_COPY_K + 1)}
    out, name, buf = set(), None, []

    def judge(n, s):
        if n in cat and IS_COPY_LEN[0] <= len(s) <= IS_COPY_LEN[1]:
            ks = [s[i:i + IS_COPY_K] for i in range(len(s) - IS_COPY_K + 1)]
            if ks and sum(k in km for k in ks) >= IS_COPY_FRAC * len(ks):
                out.add(n)
    for line in open(fa):
        if line.startswith(">"):
            if name:
                judge(name, "".join(buf).upper())
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        judge(name, "".join(buf).upper())
    return out


def _load_sv_intervals():
    """bin/sv_intervals.py, for its same-event rule; one definition, not two."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "sv_intervals", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     "sv_intervals.py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", default="",
                    help="the dense matrix. Optional now: --states-dir reads "
                         "the same information from the sparse per-sample "
                         "files, which is what should be used at scale.")
    ap.add_argument("--states-dir", default="",
                    help="read the small-variant block from the sparse "
                         "per-sample state files instead of a dense "
                         "matrix.tsv. Needs --keys.")
    ap.add_argument("--keys", default="",
                    help="the key table, with --states-dir or --states-array")
    ap.add_argument("--states-array", default="",
                    help="directory holding states.u8.npy from p5_matrix.py: "
                         "the small-variant block read from the memory-mapped "
                         "array. Needs --keys. Preferred over both other inputs")
    ap.add_argument("--n-shards", type=int, default=1,
                    help="split the output into this many shards; each shard "
                         "owns a contiguous H37Rv position range and a "
                         "contiguous range of node contigs")
    ap.add_argument("--shard", type=int, default=-1,
                    help="write only this shard (0-based) as parts under "
                         "--part-prefix; assemble them with --assemble")
    ap.add_argument("--part-prefix", default="",
                    help="path prefix for shard parts: <prefix>.<i>.h37rv.gz, "
                         "<prefix>.<i>.node.gz, <prefix>.<i>.meta.json")
    ap.add_argument("--assemble", action="store_true",
                    help="join --n-shards shards under --part-prefix into --out")
    ap.add_argument("--refmap", default="",
                    help="sample order, with --states-dir")
    ap.add_argument("--sv-matrix", default="")
    ap.add_argument("--accessory-presence", default="",
                    help="directory of <sample>.presence.tsv from "
                         "accessory/bin/locus_presence.py. Emits LEVEL 1: one "
                         "biallelic record per accessory locus, every sample "
                         "stated. Needs --accessory-catalogue.")
    # no default: accessory/assets/ held CX333's hand-placed catalogue;
    # p5_finish.sh passes the build's. Required with --accessory-presence.
    ap.add_argument("--accessory-catalogue", default="",
                    help="<build>/assets/accessory_catalogue.tsv")
    ap.add_argument("--sv-intervals", default="",
                    help="the interval catalogue from bin/sv_intervals.py. With "
                         "it the DELETION block is built from catalogued "
                         "intervals genotyped in every sample, instead of from "
                         "caller clusters. Needs --sv-interval-states.")
    ap.add_argument("--sv-interval-states", default="",
                    help="svgt_iv_states.tsv from p5_svgt.sh in interval mode, "
                         "keyed on `interval`")
    ap.add_argument("--is6110-keys", default="")
    ap.add_argument("--inherited-sv-states", default="",
                    help="inherited_sv_states.tsv from p5_svgt.sh --merge "
                         "(bin/p5_inherited_sv.py): each sample's matched-frame "
                         "genotype of the insertions it inherits from its "
                         "reference (HANDOFF ABS-2). With it, an inherited "
                         "insertion is written as measured -- ALT, REF or "
                         "NOCALL with its reason -- instead of present in every "
                         "sample matched to that reference. Needs --sv-members.")
    ap.add_argument("--sv-members", default="",
                    help="sv_matrix.members.tsv from bin/p5_sv_matrix.py: "
                         "which of a sample's P4b rows, called or inherited, "
                         "each caller-matrix cell came from")
    ap.add_argument("--sv-states", default="",
                    help="per-sample SV states from bin/p5_sv_genotype.py. The "
                         "SV matrix never writes REF, on the measured grounds "
                         "that a missing caller call is usually a missed call; "
                         "this table carries the deletions whose absence was "
                         "measured directly as read depth across the interval.")
    ap.add_argument("--is6110-states", default="",
                    help="stage-2 states from is6110/bin/is6110_p5_stage2.py. "
                         "The key table alone records only the cells the arm "
                         "reported, so every other sample lands on NOCALL -- "
                         "98.7% of the block on scale200, which makes the "
                         "insertion sites untestable by any method that needs "
                         "to contrast presence with absence. Stage 2 asks each "
                         "non-carrier's own element-free alignment whether it "
                         "looked and found nothing, and earns REF or ABSENT.")
    ap.add_argument("--ancestral", default="",
                    help="per-site ancestral allele from bin/panel_ancestral.py (<build>/assets/ancestral.tsv); "
                         "absent means no AA is emitted")
    ap.add_argument("--panel-polarity", default="",
                    help="<build>/assets/panel_polarity.tsv, from the same "
                         "reconstruction: AA_STATE, AA_CHANGES and AA_FLAG for "
                         "records with no AA (indels, MNPs, sequence-resolved "
                         "SVs, symbolic SVs through their SMALLALLELE key), so "
                         "AA_SOURCE says reconstructed there. Without it every "
                         "such record is AA_SOURCE=assumed")
    ap.add_argument("--is6110-gff", default="",
                    help="H37Rv IS6110 copies (<build>/assets/is6110_elements.gff), "
                         "for the copy-loss rule of the harmonization pass")
    ap.add_argument("--no-harmonize", action="store_true",
                    help="skip bin/harmonize_records.py (one event, one record)")
    ap.add_argument("--h37rv-contig", default="NC_000962.3")
    ap.add_argument("--h37rv-length", type=int, default=4411532)
    ap.add_argument("--h37rv-fasta", default="",
                    help="H37Rv FASTA, for the REF base of symbolic records "
                         "(SV, IS6110, accessory). Without it those records "
                         "carry REF=N, which does not match the reference")
    ap.add_argument("--cohort-name", required=True)
    # no default: it was 7713a8d71d8e (CX333), so a caller that omitted it
    # stamped another build's VCF with CX333's id
    ap.add_argument("--build-id", required=True)
    ap.add_argument("--bgzip", default=os.environ.get("MTB_BGZIP", "bgzip"))
    # MTB_TABIX is not defined in config/project_env.sh although MTB_BGZIP is,
    # so fall back to tabix beside bgzip rather than to bare `tabix`, which is
    # not on PATH here. The same class of omission as the MTB_ODGI one: a tool
    # the project uses, named nowhere, working only where PATH happens to help.
    _bg = os.environ.get("MTB_BGZIP", "")
    _tbx = os.path.join(os.path.dirname(_bg), "tabix") if _bg else "tabix"
    ap.add_argument("--tabix", default=os.environ.get("MTB_TABIX", _tbx))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.accessory_presence and not a.accessory_catalogue:
        sys.exit("FATAL: --accessory-presence needs --accessory-catalogue "
                 "(<build>/assets/accessory_catalogue.tsv)")
    if a.assemble:
        if not a.part_prefix:
            sys.exit("FATAL: --assemble needs --part-prefix")
        return assemble(a)
    sharded = a.shard >= 0
    if sharded and not (0 <= a.shard < a.n_shards and a.part_prefix):
        sys.exit("FATAL: --shard needs 0 <= shard < --n-shards and --part-prefix")
    if a.n_shards > 1 and not sharded:
        sys.exit("FATAL: --n-shards > 1 needs --shard (then --assemble)")
    if sharded and not a.keys:
        sys.exit("FATAL: sharding needs --keys, to place the node boundaries")

    # The ancestral allele is a property of the PANEL, reconstructed once over
    # the 333 genomes and their tree, so it is joined here as a lookup rather
    # than recomputed per cohort -- and adding it needs no cohort re-run.
    anc = {}
    if a.ancestral and not os.path.exists(a.ancestral):
        sys.exit(f"FATAL: no ancestral allele table at {a.ancestral}")
    if a.ancestral:
        for r in csv.DictReader(open(a.ancestral), delimiter="\t"):
            anc[(int(r["pos"]), r["ref"], r["alt"])] = (
                r["AA"], r["flag"], r.get("changes", ""))
    pol = load_panel_polarity(a.panel_polarity)

    if not (a.matrix or a.states_dir or a.states_array):
        sys.exit("FATAL: give --states-array, --states-dir or --matrix")
    want = [None]           # set below, once the shard is known; read lazily
    if a.states_array:
        if not a.keys:
            sys.exit("FATAL: --states-array needs --keys")
        rows, samples, _ = read_array(a.states_array, a.keys, want)
    elif a.states_dir:
        if not (a.keys and a.refmap):
            sys.exit("FATAL: --states-dir needs --keys and --refmap")
        rows, samples = read_sparse(a.states_dir, a.keys, a.refmap)
    else:
        samples, rows = iter_matrix(a.matrix)
    if not samples:
        sys.exit(f"no rows in {a.matrix or a.states_dir or a.states_array}")
    # The SV matrix is streamed, not loaded: at 10,000 isolates it is several
    # GB of dicts. Its header gives the samples now; rows are read per shard.
    sv_samples, sv_iter = iter_matrix(a.sv_matrix) if a.sv_matrix else ([], iter(()))
    order = sorted(set(samples) | set(sv_samples))
    # THE COHORT IS THE P5 SAMPLE SET, and nothing else may add a column. The
    # IS6110 tables used to append every sample they named, so a key table from
    # another cohort (the pilot's, via a fallback in p5_finish.sh) added that
    # cohort's isolates as extra columns. Rows for samples outside the cohort
    # are now counted and dropped.
    cohort = set(order)
    n_foreign = collections.Counter()

    H37 = ""
    if a.h37rv_fasta:
        hname, H37 = read_fasta_seq(a.h37rv_fasta)
        if hname != a.h37rv_contig:
            sys.exit(f"FATAL: --h37rv-fasta contig {hname!r} is not "
                     f"--h37rv-contig {a.h37rv_contig!r}")
        if len(H37) != a.h37rv_length:
            sys.exit(f"FATAL: --h37rv-fasta is {len(H37):,} bp, "
                     f"--h37rv-length says {a.h37rv_length:,}")

    def h37_base(pos):
        """REF base at a 1-based H37Rv position, or N without the FASTA."""
        return H37[pos - 1] if H37 and 0 < pos <= len(H37) else "N"

    # ---- which records this process writes ---------------------------------
    # One shard owns H37Rv positions (lo, hi] of an even split, and the node
    # contigs whose names fall in its share of the sorted node names in
    # keys.tsv. Both are contiguous in the output's sort order, which is what
    # lets --assemble concatenate the parts instead of re-sorting them.
    if sharded:
        L = a.h37rv_length
        h_lo = a.shard * L // a.n_shards
        h_hi = (a.shard + 1) * L // a.n_shards
        last = a.shard == a.n_shards - 1
        names = sorted({f"node_{k['node']}" for k in csv.DictReader(
            open(a.keys, newline=""), delimiter="\t")
            if k.get("frame") == "node" and k.get("node")})
        bounds = [names[(j * len(names)) // a.n_shards]
                  for j in range(1, a.n_shards)] if names else []

        def in_shard(chrom, pos):
            if chrom == a.h37rv_contig:
                return h_lo < pos <= h_hi or (last and pos > h_hi)
            return bisect.bisect_right(bounds, chrom) == a.shard
    else:
        def in_shard(chrom, pos):
            return True

    # Resolve the IS6110 key table from the cohort name when it is not given.
    # I have now passed the wrong path three times in one session -- the
    # pilot's table is `p1i_cohort_keys.tsv` with no cohort prefix while every
    # other cohort's is `<tag>_p1i_cohort_keys.tsv`, and a hand-written loop
    # silently produced a VCF with zero IS6110 records each time. The counts
    # looked plausible, so nothing failed. Resolving it here means the caller
    # cannot get it wrong.
    # THE UNPREFIXED TABLE BELONGS TO THE PILOT AND TO NOBODY ELSE. An earlier
    # version of this resolver fell back to is6110/results/p1i_cohort_keys.tsv
    # for any cohort whose own table was absent, which would have joined the
    # 23-isolate pilot's insertion sites onto a different cohort entirely --
    # worse than the zero-record bug it was written to prevent, because the
    # result would have looked populated. Only pilot and pilot_rerun map to it.
    KEYTAB = {"pilot": "is6110/results/p1i_cohort_keys.tsv",
              "pilot_rerun": "is6110/results/p1i_cohort_keys.tsv",
              "scale100": "is6110/results/scale_p1i_cohort_keys.tsv"}
    if not a.is6110_keys:
        cand = KEYTAB.get(a.cohort_name,
                          f"is6110/results/{a.cohort_name}_p1i_cohort_keys.tsv")
        if os.path.exists(cand):
            a.is6110_keys = cand
            print(f"  IS6110 keys resolved to {cand}")
        else:
            print(f"  no IS6110 key table for {a.cohort_name} at {cand}; "
                  f"emitting no insertion-site records")
    is6110 = collections.defaultdict(dict)
    is_meta = {}
    is_ev = collections.defaultdict(set)
    is_sc = collections.defaultdict(set)
    n_is_skip = collections.Counter()
    if a.is6110_keys and os.path.exists(a.is6110_keys):
        for r in csv.DictReader(open(a.is6110_keys), delimiter="\t"):
            k = r["key"]
            if r["sample"] not in cohort:
                n_foreign["IS6110 key table"] += 1
                continue
            # a site on a graph node its path visits more than once has no
            # key (is6110_write_vcf.py, audit P3IS-3); counted, not emitted
            if r.get("frame") == "repeat_node":
                n_is_skip["on a repeated graph node, not keyed"] += 1
                continue
            # Two rows of one sample on one key were settled by row order
            # (audit P3IS-5). The writer now refuses them; here the same
            # state is one cell and a disagreement stops the run.
            prev = is6110[k].get(r["sample"])
            if prev is not None:
                if prev != r["state"]:
                    sys.exit(f"FATAL: {r['sample']} has rows {prev} and "
                             f"{r['state']} on IS6110 key {k}; regenerate "
                             f"{a.is6110_keys}")
                n_is_skip["duplicate (sample, key) row, same state"] += 1
            # only the state is ever read back, so only the state is kept: a
            # copied row per cell is ~15 fields x keys x samples of dicts
            is6110[k][r["sample"]] = r["state"]
            is_meta.setdefault(k, r)
            is_ev[k].add(r.get("evidence", ""))
            is_sc[k].add(r.get("site_class", ""))
        # EVIDENCE and SITECLASS describe the key, not its first carrier
        # (audit P3IS-6): one value where every carrier agrees, else "mixed"
        for k, m in is_meta.items():
            is_meta[k] = dict(
                m, evidence=next(iter(is_ev[k])) if len(is_ev[k]) == 1 else "mixed",
                site_class=next(iter(is_sc[k])) if len(is_sc[k]) == 1 else "mixed")
        if n_is_skip:
            print("  IS6110 key table: " + ", ".join(
                f"{k} {v:,}" for k, v in sorted(n_is_skip.items())))

    # ---- stage-2 states, if they exist.
    # The arm's key table and P5's key space spell the same site differently:
    # `h37rv:1026906` there, `h37rv:1026906:C><INS>` here, because P5 keys carry
    # the allele. Joining is therefore on the key with its allele suffix
    # removed, which is the first two colon-fields for an H37Rv key and the
    # first three for a node key.
    STAGE2 = {"pilot": "is6110/results/p5_is6110_states.stage2.tsv",
              "pilot_rerun": "is6110/results/p5_is6110_states.stage2.tsv",
              "scale100": "is6110/results/scale_p5_is6110_states.stage2.tsv"}
    if not a.is6110_states:
        cand = STAGE2.get(
            a.cohort_name,
            f"is6110/results/{a.cohort_name}_p5_is6110_states.stage2.tsv")
        if os.path.exists(cand):
            a.is6110_states = cand

    def arm_key(k):
        f = k.split(":")
        n = 3 if f and f[0] == "node" else 2
        return ":".join(f[:n])

    def is6110_place(m):
        """(chrom, pos) of an IS6110 key-table row, or (None, None)."""
        if m.get("frame") == "node":
            raw = (m.get("node") or "").strip()
            if not raw:
                return None, None
            nid, _, off = raw.partition(":")
            return f"node_{nid}", int(off or 0) + 1
        p = m.get("h37rv_pos") or ""
        if not p or p in (".", "NA"):
            return None, None
        return a.h37rv_contig, int(p)

    # keys outside this shard are dropped before their stage-2 cells are read;
    # unplaceable ones are kept so exactly one shard (0) counts them
    for k in list(is_meta):
        c_, p_ = is6110_place(is_meta[k])
        if (c_ is None and a.shard > 0) or (c_ is not None and not in_shard(c_, p_)):
            del is_meta[k]; del is6110[k]

    n_promoted = collections.Counter()
    is_cell_ev = collections.defaultdict(dict)
    if a.is6110_states and os.path.exists(a.is6110_states):
        for r in csv.DictReader(open(a.is6110_states), delimiter="\t"):
            k = arm_key(r["key"])
            if r["sample"] not in cohort:
                n_foreign["IS6110 stage-2 states"] += 1
                continue
            if k not in is_meta:
                n_promoted["key not in the arm's table or not in this shard"] += 1
                continue
            prev = is6110[k].get(r["sample"])
            # Never overwrite a reported ALT with a computed state. Stage 2
            # can now state ALT itself (fix SV-EVAL-1: an off-path key whose
            # P1 H37Rv reads show IS6110 junctions, EV absent_checked_alt),
            # but only for a cell the arm left NOCALL, so a disagreement here
            # still means the projection landed somewhere the arm had called.
            if prev == "ALT" and r["state"] != "ALT":
                n_promoted["ALT kept over a computed state"] += 1
                continue
            is6110[k][r["sample"]] = r["state"]
            # FORMAT/EV: the compact code stage 2 writes beside its free-text
            # evidence (a table from before 2026-10-16 has none)
            if r.get("ev"):
                is_cell_ev[k][r["sample"]] = r["ev"]
            n_promoted[r["state"]] += 1
        print(f"  IS6110 stage-2 states from {a.is6110_states}: "
              + ", ".join(f"{k} {v:,}" for k, v in sorted(n_promoted.items())))
    else:
        print(f"  no IS6110 stage-2 states for {a.cohort_name}; every "
              f"non-reporting sample stays NOCALL, which leaves the insertion "
              f"sites untestable -- run is6110/bin/is6110_p5_stage2.py")

    # ---- records -----------------------------------------------------------
    # (chrom, pos, ident, ref, alts, info, per-sample {sample: (gt, st)})
    recs = []
    nodes_used = {}

    def place(r):
        """Return (chrom, pos) for a matrix or key row."""
        if r.get("frame") == "node":
            nid = r.get("node") or "0"
            off = int(r.get("node_offset") or 0)
            return f"node_{nid}", off + 1
        p = r.get("h37rv_pos") or ""
        if not p or p in (".", "NA"):
            return None, None
        return a.h37rv_contig, int(p)

    def add(rec):
        """Keep a record, declaring its node contig. Callers have already
        checked in_shard, before building the record's cells."""
        chrom, pos = rec[0], rec[1]
        if chrom != a.h37rv_contig:
            nid = chrom[len("node_"):]
            nodes_used[nid] = max(nodes_used.get(nid, 0), pos)
        recs.append(rec)

    # the small-variant reader skips out-of-shard keys before building cells
    want[0] = lambda k: (lambda c: c[0] is not None and in_shard(*c))(
        place({kk: vv for kk, vv in k.items() if kk != "node_offset"}))

    n_small = n_sv = n_is = 0
    for r in rows:
        chrom, pos = place(r)
        if chrom is None or not in_shard(chrom, pos):
            continue
        # The ALT lives in the KEY, not in a column. There is no
        # `canonical_alt` in the matrix -- `r.get("canonical_alt")` returned
        # None for every row, so every small variant was written with
        # ALT=<UNKNOWN>, and the ancestral-allele join could never match. Keys
        # are `h37rv:<pos>:<REF>><ALT>` and `node:<node>:<off>:<REF>><ALT>`,
        # so the allele pair is the last colon-separated field.
        ref = (r.get("canonical_ref") or "N").upper() or "N"
        alt = "."
        tail = r["key"].rsplit(":", 1)[-1]
        if ">" in tail:
            kref, kalt = tail.split(">", 1)
            if kref:
                ref = kref.upper()
            alt = kalt.upper() or "."
        states = {s: r.get(s, "") for s in samples}
        need_star = any(v == "ABSENT" for v in states.values())
        alts = [alt if alt and alt != "." else "<UNKNOWN>"]
        if need_star:
            alts.append("*")
        info = [f"CLASS=small", f"KIND={r.get('kind','')}",
                f"REGION={r.get('region','')}", f"FRAME={r.get('frame','')}"]
        hit = anc.get((pos, ref, alt)) if chrom == a.h37rv_contig else None
        if hit:
            aa, flag, nch = hit
            info.append(f"AA={aa}")
            if flag:
                info.append(f"AA_FLAG={flag}")
            if nch != "":
                info.append(f"AA_CHANGES={nch}")
            # not `elif` on the flag: since AA-1 a resolved site can carry
            # OUTGROUP or HOMOPLASTIC, and its inversion must still be marked
            if aa == alt:
                # The reference carries the DERIVED allele here, so GT=0 means
                # derived and GT=1 means ancestral -- inverted relative to the
                # naive reading. True of 8.9% of panel SNP sites.
                info.append("AA_INVERTED")
        if r.get("frame") == "node":
            info.append(f"NODE={r.get('node','')}")
        evs = r.get("_ev") or {}          # FORMAT/EV, from the evidence array
        cells = "\t".join(vcf_cell(states.get(s, ""), evs.get(s, "")) for s in order)
        add((chrom, pos, r["key"], ref, alts, info, cells))
        n_small += 1

    # ---- measured SV states, overlaid on the matrix's presence-only cells.
    # The matrix writes ALT, NOCALL and ABSENT and never REF; p5_sv_genotype.py
    # measures depth across each deletion's interval and earns the REF the
    # caller could not. A reported ALT is never overwritten.
    sv_rows = []
    for r in sv_iter:
        p = r.get("h37rv_pos") or ""
        if p and p not in (".", "NA") and in_shard(a.h37rv_contig, int(p)):
            sv_rows.append(r)
    sv_keys = {r["key"] for r in sv_rows}
    sv_state = collections.defaultdict(dict)

    # ---- inherited insertions, as measured (ABS-2) -------------------------
    # The caller matrix states an inherited insertion ALT in every sample whose
    # reference carries it: nothing ever read the sample. p5_inherited_sv.py
    # genotypes each one in the sample's own alignment to its reference, and
    # the membership table says which cells came from inheritance alone, so
    # only those are replaced; a cell any caller reached keeps its ALT.
    if bool(a.inherited_sv_states) != bool(a.sv_members):
        sys.exit("FATAL: --inherited-sv-states and --sv-members go together")
    members = collections.defaultdict(lambda: collections.defaultdict(list))
    inh_ins = {}
    if a.sv_members:
        for path in (a.sv_members, a.inherited_sv_states):
            if not os.path.exists(path):
                sys.exit(f"FATAL: no {path}")
        for q in csv.DictReader(open(a.sv_members, newline=""), delimiter="\t"):
            if q["key"] in sv_keys:
                members[q["key"]][q["sample"]].append((q["member"], q["component"]))
        for q in csv.DictReader(open(a.inherited_sv_states, newline=""),
                                delimiter="\t"):
            if q.get("kind") == "INS":
                inh_ins[(q["sample"], q["id"])] = (q["state"], q.get("evidence", ""))
        print(f"  inherited insertions genotyped in the matched frame: "
              f"{len(inh_ins):,} (sample, event) pairs from {a.inherited_sv_states}")
    PRIORITY = {"ALT": 0, "REF": 1, "NOCALL": 2}
    n_inh = collections.Counter()

    def inherited_cell(r, s, st):
        """(state, evidence) of a caller-matrix cell that inheritance alone
        made ALT; (st, "") for any other cell."""
        if st != "ALT":
            return st, ""
        mem = members.get(r["key"], {}).get(s) if members else None
        if mem is None:
            if not members and r.get("component") == "inherited":
                n_inh["inherited ALT kept, no membership table: unchecked"] += 1
                return st, "inherited_unchecked"
            return st, ""
        if any(c != "inherited" for _, c in mem):
            return st, ""                      # a caller reached it
        if r.get("svtype") != "INS":
            # an inherited deletion the catalogue does not hold: no check
            n_inh["inherited DEL outside the catalogue: unchecked"] += 1
            return st, "inherited_unchecked"
        got = [inh_ins[(s, m)] for m, _ in mem if (s, m) in inh_ins]
        if not got:
            n_inh["inherited INS with no matched-frame genotype: unchecked"] += 1
            return st, "inherited_unchecked"
        new, ev = min(got, key=lambda g: PRIORITY.get(g[0], 3))
        n_inh[f"inherited INS ALT -> {new}"] += 1
        return new, ev
    if a.sv_states and os.path.exists(a.sv_states):
        n_sv_state = collections.Counter()
        rdr = csv.DictReader(open(a.sv_states), delimiter="\t")
        if "key" not in (rdr.fieldnames or []):
            sys.exit(f"FATAL: {a.sv_states} has no `key` column; its columns are "
                     f"{rdr.fieldnames}. --sv-states overlays the CALLER matrix "
                     f"and is keyed on `key`. A table keyed on `interval` is the "
                     f"interval arm's and belongs to --sv-interval-states.")
        for r in rdr:
            if r["key"] not in sv_keys:
                continue                      # another shard's, or no record
            sv_state[r["key"]][r["sample"]] = r["state"]
            n_sv_state[r["state"]] += 1
        print(f"  SV states from {a.sv_states}: "
              + ", ".join(f"{k} {v:,}" for k, v in sorted(n_sv_state.items())))
    elif a.sv_matrix:
        print(f"  no measured SV states for {a.cohort_name}; the SV block will "
              f"carry no REF, which leaves it untestable -- run "
              f"bin/p5_svgt.sh")

    # ---- THE DELETION BLOCK, FROM THE INTERVAL CATALOGUE -------------------
    # WHY THIS EXISTS. The SV block used to be caller clusters with depth states
    # overlaid, and the rebuilt interval arm wrote to a key space -- `svi:` on
    # 3,658 catalogued intervals -- that this merge never read, because it reads
    # `sv:` on caller clusters. Zero keys in common. Every improvement to that
    # arm therefore landed in a file the VCF, the event matrix and the
    # association scan all ignored.
    #
    # With --sv-intervals the deletions come from the catalogue: one record per
    # interval, genotyped in EVERY sample by depth and, where the two-frame arm
    # ran, confirmed against the H37Rv frame. The caller matrix still supplies
    # the INSERTIONS, because depth cannot genotype an insertion and the
    # catalogue is deletions only -- dropping them would lose 2,178 records on
    # scale200 rather than improve anything. Those stay presence-only and say so
    # in COMPONENT.
    iv_states = collections.defaultdict(dict)
    iv_ev = collections.defaultdict(dict)
    ivs = []
    if a.sv_intervals and a.sv_interval_states:
        # catalogue first, so only this shard's intervals have cells loaded
        ivs = [r for r in csv.DictReader(open(a.sv_intervals, newline=""),
                                         delimiter="\t")
               if in_shard(a.h37rv_contig, max(int(r["start"]) - 1, 1))]
        iv_ids = {r["interval"] for r in ivs}
        for q in csv.DictReader(open(a.sv_interval_states, newline=""),
                                delimiter="\t"):
            if q["interval"] in iv_ids:
                iv_states[q["interval"]][q["sample"]] = q["state"]
                # FORMAT/EV: the reason behind an ABSENT or NOCALL, and how
                # an inherited deletion was decided (ABS-1, ABS-2)
                evq = q.get("evidence") or ""
                if evq and (q["state"] in ("ABSENT", "NOCALL")
                            or evq.startswith(("inherited", "absent"))):
                    iv_ev[q["interval"]][q["sample"]] = evq
        n_st = collections.Counter(v for d in iv_states.values()
                                   for v in d.values())
        t = sum(n_st.values()) or 1
        print(f"  deletion block from {len(ivs):,} catalogued intervals, "
              f"{t:,} genotyped cells: "
              + "  ".join(f"{k} {v:,} ({v/t:.1%})" for k, v in sorted(n_st.items())))
    elif a.sv_intervals or a.sv_interval_states:
        sys.exit("FATAL: --sv-intervals and --sv-interval-states go together")

    # WHICH CALLER DELETIONS THE CATALOGUE SUPERSEDES: only those an interval
    # actually covers, under the catalogue's own same-event rule. Every caller
    # deletion used to be dropped once a catalogue was given, so a deletion the
    # catalogue never took in -- outside its 50 bp-100 kb range, or from a
    # cohort whose caller set built no catalogue -- vanished from the VCF
    # without a count. The whole catalogue is read for this, not just this
    # shard's rows, so a deletion near a shard edge sees its interval.
    cat_starts, cat_rows = [], []
    if ivs:
        _ivm = _load_sv_intervals()
        cat_rows = sorted((int(r["start"]), abs(int(r.get("svlen") or 0)))
                          for r in csv.DictReader(open(a.sv_intervals, newline=""),
                                                  delimiter="\t"))
        cat_starts = [x[0] for x in cat_rows]

    def catalogued(r):
        try:
            p0, L = int(r["h37rv_pos"]), abs(int(r.get("svlen") or 0))
        except (ValueError, KeyError):
            return False
        lo = bisect.bisect_left(cat_starts, p0 - _ivm.TOL_CAP)
        hi = bisect.bisect_right(cat_starts, p0 + _ivm.TOL_CAP)
        # same_deletion: position, length AND half the bases shared (review
        # 2, R2-GENO-2); h37rv_pos is the anchor, the catalogue's start the
        # first deleted base
        return any(_ivm.same_deletion(p0 + 1, L, s0, l0)
                   for s0, l0 in cat_rows[lo:hi])
    n_uncat = 0

    n_iv = 0
    for r in ivs:
        st = iv_states.get(r["interval"], {})
        if not st:
            continue
        info = [f"CLASS=sv", f"SVTYPE={r.get('svtype','DEL')}",
                f"SVLEN=-{abs(int(r.get('svlen') or 0))}",
                f"COMPONENT=interval",
                f"SVSOURCE={r.get('source','')}",
                f"SVTIER={r.get('support_tier','')}",
                f"PANELCARRIERS={r.get('n_ref_carriers','0')}"]
        if str(r.get("is6110_prox", "0")) == "1":
            info.append("IS6110PROX")
        evs = iv_ev.get(r["interval"], {})
        cell = "\t".join(vcf_cell(st.get(s, ""), evs.get(s, "")) for s in order)
        # `start` is the first DELETED base (sv_intervals.py), and a symbolic
        # VCF deletion is anchored on the base BEFORE it, with END the last
        # deleted base. Writing POS=start put every deletion one base late and,
        # without END, a region query inside the deletion missed the record.
        st0, en0 = int(r["start"]), int(r.get("end") or r["start"])
        pos = max(st0 - 1, 1)
        info.append(f"END={en0}")
        add((a.h37rv_contig, pos, r["interval"], h37_base(pos),
             [f"<{r.get('svtype','DEL')}>"], info, cell))
        n_iv += 1
    if ivs:
        print(f"  {n_iv:,} deletion records written from the catalogue")

    # ---- LEVEL 1: does the sample carry the insert at all? ----------------
    # Michael's two-level framing. An accessory locus carries two different
    # characters and the pipeline had been emitting only the second: variation
    # INSIDE the insert, which for a sample that does not have the insert is
    # neither reference nor unknown but INAPPLICABLE. The first character --
    # whether the insert is there -- is stated here, one record per locus.
    #
    # Polarity follows the coordinate frame. Carrying the insert is the
    # DERIVED state and gets the ALT allele; a convergence test then looks for
    # repeated independent gains of the insert, which is the event of
    # interest.
    #
    # EVERY CATALOGUE LOCUS IS WRITTEN (fix SV-EVAL-2, 2026-10-16). Loci whose
    # sequence H37Rv carries (locus_presence.read_route_blind: 465 of the P200
    # build's 564) used to be skipped as unmeasurable, and with them every
    # record whose cells were all NOCALL, so new1k_p200's VCF held 99 loci.
    # Presence tables from locus_presence.py version 2 (a `route` column)
    # measure those loci at their H37Rv anchor; their states are written as
    # measured. A table from before has no measurement for them, so their
    # cells are NOCALL "read_route_blind" -- the same no-ABSENT guarantee as
    # audit P3IS-2 -- but the record is still written. Every NOCALL carries
    # its reason in FORMAT/EV: the table's evidence_code (version 2), else
    # "uncertain" / "read_route_blind", and "no_presence_row" where the
    # sample has no row for the locus.
    #
    # IS6110 COPIES. A locus whose sequence is IS6110 (is6110_copy_loci) is
    # the same event as the IS6110 arm's insertion record at that site. It is
    # written and flagged ACCIS6110 so a consumer can count it once; it is
    # not dropped here, because the IS6110 arm keys its sites independently
    # and a silent drop would lose any site only this block has (a decision
    # left to the user; see HANDOFF SV-EVAL-2).
    #
    # Emitted BEFORE the caller SV rows so those rows can tell which insertion
    # events this block already states (audit P4P5-6).
    n_l1 = n_l1_blind = n_l1_norow = 0
    acc_emitted = []                 # (pos, length) of every level-1 record
    if a.accessory_presence:
        import importlib.util
        _spec = importlib.util.spec_from_file_location(
            "locus_presence", os.path.join(
                os.path.dirname(os.path.abspath(__file__)), "..", "accessory",
                "bin", "locus_presence.py"))
        _lp = importlib.util.module_from_spec(_spec)
        _spec.loader.exec_module(_lp)
        cat_all = list(csv.DictReader(open(a.accessory_catalogue, newline=""),
                                      delimiter="\t"))
        if not cat_all or "locus_id" not in cat_all[0]:
            sys.exit(f"FATAL: {a.accessory_catalogue} has no locus_id rows")
        n_nopos = sum(1 for r in cat_all if not (r.get("pos") or "").isdigit())
        cat = {r["locus_id"]: r for r in cat_all
               if (r.get("pos") or "").isdigit()
               and in_shard(a.h37rv_contig, int(r["pos"]))}
        pres = collections.defaultdict(dict)   # locus -> sample -> (state, code, v2)
        pf = sorted(glob.glob(os.path.join(a.accessory_presence,
                                           "*.presence.tsv")))
        n_v2 = n_old = 0
        for f in pf:
            with open(f, newline="") as fh:
                rd = csv.DictReader(fh, delimiter="\t")
                v2 = "route" in (rd.fieldnames or [])
                n_v2 += v2
                n_old += not v2
                for r in rd:
                    if r["locus"] in cat:     # this shard's loci only
                        pres[r["locus"]][r["sample"]] = (
                            r["state"], r.get("evidence_code") or "", v2)
        got = {s2 for d in pres.values() for s2 in d}
        print(f"  level 1: {len(pf):,} presence tables ({n_v2:,} version 2, "
              f"{n_old:,} older), {len(pres):,} of {len(cat):,} loci, "
              f"{len(got):,} samples ({len(got & set(order)):,} in this cohort)")
        is_copy = is6110_copy_loci(a.accessory_catalogue, H37, a.is6110_gff, cat)
        L1 = {"PRESENT": "ALT", "ABSENT": "REF", "UNCERTAIN": "NOCALL"}
        n_st = collections.Counter()
        n_ev = collections.Counter()
        for lid, m in sorted(cat.items(), key=lambda x: (int(x[1]["pos"]), x[0])):
            per = pres.get(lid, {})
            blind = _lp.read_route_blind(m)
            cells, any_old_blind = [], False
            for s2 in order:
                hit = per.get(s2)
                if hit is None:
                    st, ev = "NOCALL", "no_presence_row"
                else:
                    raw, code, v2 = hit
                    if not v2 and blind:
                        # an older table cannot have measured it, whatever it says
                        st, ev = "NOCALL", "read_route_blind"
                        any_old_blind = True
                    else:
                        st = L1.get(raw, "NOCALL")
                        ev = (code or raw.lower() or "no_state") if st == "NOCALL" else ""
                n_st[st] += 1
                if ev:
                    n_ev[ev] += 1
                cells.append(vcf_cell(st, ev))
            n_l1_blind += any_old_blind
            n_l1_norow += not per
            acc_len = int(m.get("graph_len") or m.get("rep_len") or 0)
            info = ["CLASS=accessory_presence", "COMPONENT=level1",
                    f"ACCLOCUS={lid}",
                    f"ACCLEN={acc_len}",
                    f"ACCKLASS={m.get('klass','')}",
                    f"PANELCARRIERS={m.get('graph_carriers') or m.get('carriers_any') or 0}"]
            if str(m.get("route", "")):
                info.append(f"ACCROUTE={m['route']}")
            if lid in is_copy:
                info.append("ACCIS6110")
            add((a.h37rv_contig, int(m["pos"]), f"acc:{lid}",
                 h37_base(int(m["pos"])), ["<INS>"], info, "\t".join(cells)))
            acc_emitted.append((int(m["pos"]), acc_len))
            n_l1 += 1
        t = sum(n_st.values()) or 1
        print(f"  {n_l1:,} level-1 records ({len(is_copy):,} flagged ACCIS6110): "
              + "  ".join(f"{k} {v:,} ({v/t:.1%})"
                          for k, v in sorted(n_st.items())))
        if n_ev:
            print("    NOCALL reasons: " + ", ".join(
                f"{k} {v:,}" for k, v in sorted(n_ev.items())))
        if n_l1_blind:
            print(f"    {n_l1_blind:,} loci whose sequence H37Rv carries, read "
                  f"from tables older than presence version 2: written NOCALL "
                  f"read_route_blind; rerun the presence step to measure them")
        if n_l1_norow:
            print(f"    {n_l1_norow:,} loci no presence table has a row for")
        if n_nopos and a.shard <= 0:
            print(f"    {n_nopos:,} catalogue loci have no H37Rv position and "
                  f"cannot be written")
    acc_emitted.sort()
    acc_starts = [x[0] for x in acc_emitted]

    # ONE INSERTION, ONE RECORD (audit P4P5-6). A caller insertion that the
    # IS6110 arm or the level-1 block above already states is the same event
    # written a second time, and association and parsimony would count it
    # twice: scale200 had 143 IS6110-length caller INS records, 100 within
    # 20 bp of an IS6110 record, and 328 caller INS records at the very
    # position of an accessory locus. The dedicated arm's record is kept --
    # it states every sample, where the caller row is presence-only -- and
    # the caller row is dropped and counted. Positions are this shard's.
    IS_LEN, IS_WIN, ACC_WIN = (1340, 1370), 20, 10
    is_starts = sorted(int(m["h37rv_pos"]) for m in is_meta.values()
                       if m.get("frame") != "node"
                       and (m.get("h37rv_pos") or "").isdigit())

    def near(starts, p0, w):
        i = bisect.bisect_left(starts, p0 - w)
        return i < len(starts) and starts[i] <= p0 + w

    def dup_of(p0, L):
        if IS_LEN[0] <= L <= IS_LEN[1] and near(is_starts, p0, IS_WIN):
            return "an IS6110 record"
        i = bisect.bisect_left(acc_starts, p0 - ACC_WIN)
        while i < len(acc_emitted) and acc_emitted[i][0] <= p0 + ACC_WIN:
            al = acc_emitted[i][1]
            if al and L and max(al, L) / min(al, L) <= 1.5:
                return "an accessory-presence record"
            i += 1
        return ""
    n_dup = collections.Counter()

    for r in sv_rows:
        p = r.get("h37rv_pos") or ""
        if not p or p in (".", "NA"):
            continue
        svtype = r.get("svtype", "SV")
        if svtype == "INS":
            try:
                why = dup_of(int(p), abs(int(r.get("svlen") or 0)))
            except ValueError:
                why = ""
            if why:
                n_dup[why] += 1
                continue
        # With the catalogue in play the caller's deletions it covers are
        # superseded by it; insertions, and deletions it does not cover, are
        # still written from the caller matrix.
        uncat = False
        if ivs and svtype == "DEL":
            if catalogued(r):
                continue
            uncat = True
            n_uncat += 1
        info = [f"CLASS=sv", f"SVTYPE={svtype}",
                f"SVLEN={r.get('svlen','0')}",
                f"COMPONENT={r.get('component','')}",
                f"QUALBAND={r.get('qual_band','')}"]
        if uncat:
            info.append("UNCATALOGUED")
        measured = sv_state.get(r["key"], {})
        cell = []
        for s in order:
            st = r.get(s, "")
            if st != "ALT" and measured.get(s):
                st = measured[s]
            st, ev = inherited_cell(r, s, st)
            cell.append(vcf_cell(st, ev))
        if svtype == "DEL" and str(r.get("svlen") or "").lstrip("-").isdigit():
            info.append(f"END={int(p) + abs(int(r['svlen']))}")
        add((a.h37rv_contig, int(p), r["key"], h37_base(int(p)),
             [f"<{svtype}>"], info, "\t".join(cell)))
        n_sv += 1
    if ivs:
        print(f"  {n_uncat:,} caller deletions no catalogue interval covers, "
              f"kept presence-only (UNCATALOGUED)")
    for why, n in sorted(n_inh.items()):
        print(f"  {why}: {n:,} cells")
    for why, n in sorted(n_dup.items()):
        print(f"  {n:,} caller insertions dropped as the same event as {why}")



    n_unplaced = 0
    for k, per in is6110.items():
        m = is_meta[k]
        # The IS6110 key table's `node` column is "<node>:<offset>", not a bare
        # node id, and it has no `node_offset` column at all -- so passing the
        # row straight to place() produced contig names like `node_46502:0`,
        # which is not a legal VCF contig ID, for 175 of scale200's contigs.
        # And 319 of its 963 node-frame rows have an EMPTY node and a key of
        # literally `node:`: those sites are unplaceable in either frame, and
        # they were all collapsing onto one bogus `node_0` record at position 1.
        # A VCF record needs a position; a site with none is counted and
        # excluded rather than invented.
        if m.get("frame") == "node":
            raw = (m.get("node") or "").strip()
            if not raw:
                n_unplaced += 1
                continue
            nid, _, off = raw.partition(":")
            m = dict(m, node=nid, node_offset=off or "0")
        chrom, pos = place(m)
        if chrom is None:
            n_unplaced += 1
            continue
        info = ["CLASS=is6110", "SVTYPE=INS", "SVLEN=1355",
                f"MEINFO=IS6110,1,1355,.",     # strand not determined (P3IS-6)
                f"EVIDENCE={m.get('evidence','')}",
                f"SITECLASS={m.get('site_class','')}",
                f"FRAME={m.get('frame','')}"]
        if m.get("frame") == "node":
            info.append(f"NODE={m.get('node','')}")
        cell = []
        evk = is_cell_ev.get(k, {})
        for s in order:
            st = per.get(s) or "NOCALL"
            cell.append(vcf_cell(st, evk.get(s, "")))
        ref_b = h37_base(pos) if chrom == a.h37rv_contig else "N"
        add((chrom, pos, k, ref_b, ["<INS:ME:IS6110>"], info,
             "\t".join(cell)))
        n_is += 1

    # ONE EVENT, ONE RECORD across the arms (bin/harmonize_records.py,
    # approved 2026-10-11): HaplotypeCaller owns indels under 50 bp, the SV
    # block owns 50 bp and over, the IS6110 arm owns its insertions, and
    # deletions of an H37Rv IS6110 copy are labelled as copy losses.
    if not a.no_harmonize:
        copies = _harm.read_is6110_copies(a.is6110_gff) if a.is6110_gff else []
        if not a.is6110_gff:
            print("  no --is6110-gff: IS6110 copy losses not labelled")
        recs, n_harm = _harm.harmonize(recs, a.h37rv_contig, copies)
        for why, k in sorted(n_harm.items()):
            print(f"  harmonize: {why}: {k:,}")
    # AA-1: AA_SOURCE on every record (after harmonization, which relabels
    # >= 50 bp small alleles as CLASS=sv)
    recs, n_src = annotate_aa_source(recs, a.h37rv_contig, pol)
    print("  ancestral state source: " + ", ".join(
        f"{k} {v:,}" for k, v in sorted(n_src.items())))
    recs.sort(key=lambda x: (x[0] != a.h37rv_contig, x[0], x[1]))

    # ---- write -------------------------------------------------------------
    n_star_added = 0
    if sharded:
        # Record lines only; --assemble writes the one header. H37Rv and node
        # records go to separate parts so the assembled file keeps all H37Rv
        # records ahead of every node contig, as the single-file output does.
        parts = {}
        for part in ("h37rv", "node"):
            fh = open(f"{a.part_prefix}.{a.shard}.{part}.gz.tmp", "wb")
            parts[part] = (fh, subprocess.Popen([a.bgzip, "-c"],
                                                stdin=subprocess.PIPE, stdout=fh))
        for rec in recs:
            line, added = record_line(rec)
            n_star_added += added
            part = "h37rv" if rec[0] == a.h37rv_contig else "node"
            parts[part][1].stdin.write(line.encode())
        for part, (fh, bz) in parts.items():
            bz.stdin.close()
            if bz.wait() != 0:
                sys.exit(f"FATAL: bgzip failed writing shard {a.shard} {part}")
            fh.close()
            os.replace(f"{a.part_prefix}.{a.shard}.{part}.gz.tmp",
                       f"{a.part_prefix}.{a.shard}.{part}.gz")
    else:
        tmp = a.out[:-3] if a.out.endswith(".gz") else a.out
        with open(tmp, "w") as fh:
            write_header(fh.write, a, nodes_used, order)
            for rec in recs:
                line, added = record_line(rec)
                n_star_added += added
                fh.write(line)
        if a.out.endswith(".gz"):
            subprocess.run([a.bgzip, "-f", tmp], check=True)
            # check=True: an index that failed to build left a VCF that region
            # queries silently cannot read
            subprocess.run([a.tabix, "-f", "-p", "vcf", a.out], check=True)

    print(f"  cohort {a.cohort_name}: {len(order)} samples, {len(recs):,} records")
    print(f"    small variants {n_small:,}   structural {n_sv + n_iv:,} "
          f"({n_iv:,} catalogued intervals, {n_sv:,} caller rows)   "
          f"IS6110 {n_is:,}")
    n_aa = sum(1 for r in recs if any(i.startswith("AA=") for i in r[5]))
    n_inv = sum(1 for r in recs if "AA_INVERTED" in r[5])
    if sharded:
        meta = dict(shard=a.shard, n_shards=a.n_shards, order=order,
                    cohort_name=a.cohort_name, build_id=a.build_id,
                    h37rv_contig=a.h37rv_contig, nodes_used=nodes_used,
                    counts=dict(records=len(recs), small=n_small, sv=n_sv,
                                iv=n_iv, is6110=n_is, l1=n_l1, aa=n_aa,
                                aa_inverted=n_inv, star_added=n_star_added,
                                is6110_unplaced=n_unplaced))
        tmpm = f"{a.part_prefix}.{a.shard}.meta.json.tmp"
        with open(tmpm, "w") as fh:
            json.dump(meta, fh)
        os.replace(tmpm, f"{a.part_prefix}.{a.shard}.meta.json")
    print(f"    contigs: {a.h37rv_contig} plus {len(nodes_used):,} node_* contigs")
    print(f"    ancestral allele on {n_aa:,} records, of which {n_inv:,} have "
          f"the reference carrying the derived allele")
    if n_star_added:
        print(f"    {n_star_added:,} symbolic records carry * for ABSENT cells")
    if not H37:
        print("    WARNING: no --h37rv-fasta; symbolic records carry REF=N")
    for what, n in sorted(n_foreign.items()):
        print(f"    {n:,} rows DROPPED from the {what}: sample not in this "
              f"cohort's P5 sample set")
    if n_unplaced:
        print(f"    {n_unplaced:,} IS6110 sites EXCLUDED: no position in either "
              f"frame, so no VCF record is possible")
    print(f"  -> {f'{a.part_prefix}.{a.shard}.*' if sharded else a.out}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
