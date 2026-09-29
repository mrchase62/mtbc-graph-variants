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
import argparse, collections, csv, glob, gzip, os, subprocess, sys, datetime

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
    M = load_states(state_dir, keys, samples)
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
            yield r
    return gen(), samples

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
                    help="the key table, with --states-dir")
    ap.add_argument("--refmap", default="",
                    help="sample order, with --states-dir")
    ap.add_argument("--sv-matrix", default="")
    ap.add_argument("--accessory-presence", default="",
                    help="directory of <sample>.presence.tsv from "
                         "accessory/bin/locus_presence.py. Emits LEVEL 1: one "
                         "biallelic record per accessory locus, every sample "
                         "stated. Needs --accessory-catalogue.")
    ap.add_argument("--accessory-catalogue",
                    default="accessory/assets/accessory_catalogue.tsv")
    ap.add_argument("--sv-intervals", default="",
                    help="the interval catalogue from bin/sv_intervals.py. With "
                         "it the DELETION block is built from catalogued "
                         "intervals genotyped in every sample, instead of from "
                         "caller clusters. Needs --sv-interval-states.")
    ap.add_argument("--sv-interval-states", default="",
                    help="svgt_iv_states.tsv from p5_svgt.sh in interval mode, "
                         "keyed on `interval`")
    ap.add_argument("--is6110-keys", default="")
    ap.add_argument("--sv-states", default="",
                    help="per-sample SV states from bin/p5_sv_genotype.py. The "
                         "SV matrix never writes REF, on the measured grounds "
                         "that a missing delly call is usually a missed call; "
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
    ap.add_argument("--ancestral", default="data/trees/cx333.ancestral.tsv",
                    help="per-site ancestral allele from bin/ancestral_alleles.py; "
                         "absent means no AA is emitted")
    ap.add_argument("--h37rv-contig", default="NC_000962.3")
    ap.add_argument("--h37rv-length", type=int, default=4411532)
    ap.add_argument("--h37rv-fasta", default="",
                    help="H37Rv FASTA, for the REF base of symbolic records "
                         "(SV, IS6110, accessory). Without it those records "
                         "carry REF=N, which does not match the reference")
    ap.add_argument("--cohort-name", required=True)
    ap.add_argument("--build-id", default="7713a8d71d8e")
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

    # The ancestral allele is a property of the PANEL, reconstructed once over
    # the 333 genomes and their tree, so it is joined here as a lookup rather
    # than recomputed per cohort -- and adding it needs no cohort re-run.
    anc = {}
    if a.ancestral and os.path.exists(a.ancestral):
        for r in csv.DictReader(open(a.ancestral), delimiter="\t"):
            anc[(int(r["pos"]), r["ref"], r["alt"])] = (r["AA"], r["flag"])

    if not (a.matrix or a.states_dir):
        sys.exit("FATAL: give --matrix or --states-dir")
    if a.states_dir:
        if not (a.keys and a.refmap):
            sys.exit("FATAL: --states-dir needs --keys and --refmap")
        rows, samples = read_sparse(a.states_dir, a.keys, a.refmap)
    else:
        samples, rows = iter_matrix(a.matrix)
    if not samples:
        sys.exit(f"no rows in {a.matrix or a.states_dir}")
    sv_rows, sv_samples = read_matrix(a.sv_matrix) if a.sv_matrix else ([], [])
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
    if a.is6110_keys and os.path.exists(a.is6110_keys):
        for r in csv.DictReader(open(a.is6110_keys), delimiter="\t"):
            k = r["key"]
            if r["sample"] not in cohort:
                n_foreign["IS6110 key table"] += 1
                continue
            is6110[k][r["sample"]] = r
            is_meta.setdefault(k, r)

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

    n_promoted = collections.Counter()
    if a.is6110_states and os.path.exists(a.is6110_states):
        for r in csv.DictReader(open(a.is6110_states), delimiter="\t"):
            k = arm_key(r["key"])
            if r["sample"] not in cohort:
                n_foreign["IS6110 stage-2 states"] += 1
                continue
            if k not in is_meta:
                n_promoted["key not in the arm's table"] += 1
                continue
            prev = is6110[k].get(r["sample"])
            # Never overwrite a reported ALT with a computed absence; stage 2
            # promotes nothing to ALT, so a disagreement here would mean the
            # projection landed somewhere the arm had already called.
            if prev and prev.get("state") == "ALT" and r["state"] != "ALT":
                n_promoted["ALT kept over a computed state"] += 1
                continue
            is6110[k][r["sample"]] = dict(prev or is_meta[k],
                                          sample=r["sample"], state=r["state"])
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
        """Return (chrom, pos) for a matrix or key row, declaring node contigs."""
        if r.get("frame") == "node":
            nid = r.get("node") or "0"
            off = int(r.get("node_offset") or 0)
            nodes_used.setdefault(nid, 0)
            nodes_used[nid] = max(nodes_used[nid], off + 1)
            return f"node_{nid}", off + 1
        p = r.get("h37rv_pos") or ""
        if not p or p in (".", "NA"):
            return None, None
        return a.h37rv_contig, int(p)

    n_small = n_sv = n_is = 0
    for r in rows:
        chrom, pos = place(r)
        if chrom is None:
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
            aa, flag = hit
            info.append(f"AA={aa}")
            if flag:
                info.append(f"AA_FLAG={flag}")
            elif aa == alt:
                # The reference carries the DERIVED allele here, so GT=0 means
                # derived and GT=1 means ancestral -- inverted relative to the
                # naive reading. True of 8.9% of panel SNP sites.
                info.append("AA_INVERTED")
        if r.get("frame") == "node":
            info.append(f"NODE={r.get('node','')}")
        cells = "\t".join(f'{STATE_GT.get(states.get(s, ""), ".")}:'
                          f'{states.get(s, "") or "NA"}' for s in order)
        recs.append((chrom, pos, r["key"], ref, alts, info, cells))
        n_small += 1

    # ---- measured SV states, overlaid on the matrix's presence-only cells.
    # The matrix writes ALT, NOCALL and ABSENT and never REF; p5_sv_genotype.py
    # measures depth across each deletion's interval and earns the REF the
    # caller could not. A reported ALT is never overwritten.
    sv_state = collections.defaultdict(dict)
    if a.sv_states and os.path.exists(a.sv_states):
        n_sv_state = collections.Counter()
        rdr = csv.DictReader(open(a.sv_states), delimiter="\t")
        if "key" not in (rdr.fieldnames or []):
            sys.exit(f"FATAL: {a.sv_states} has no `key` column; its columns are "
                     f"{rdr.fieldnames}. --sv-states overlays the CALLER matrix "
                     f"and is keyed on `key`. A table keyed on `interval` is the "
                     f"interval arm's and belongs to --sv-interval-states.")
        for r in rdr:
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
    ivs = []
    if a.sv_intervals and a.sv_interval_states:
        for q in csv.DictReader(open(a.sv_interval_states, newline=""),
                                delimiter="\t"):
            iv_states[q["interval"]][q["sample"]] = q["state"]
        ivs = list(csv.DictReader(open(a.sv_intervals, newline=""),
                                  delimiter="\t"))
        n_st = collections.Counter(v for d in iv_states.values()
                                   for v in d.values())
        t = sum(n_st.values()) or 1
        print(f"  deletion block from {len(ivs):,} catalogued intervals, "
              f"{t:,} genotyped cells: "
              + "  ".join(f"{k} {v:,} ({v/t:.1%})" for k, v in sorted(n_st.items())))
    elif a.sv_intervals or a.sv_interval_states:
        sys.exit("FATAL: --sv-intervals and --sv-interval-states go together")

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
        cell = "\t".join(f'{STATE_GT.get(st.get(s, ""), ".")}:{st.get(s, "") or "NA"}'
                          for s in order)
        # `start` is the first DELETED base (sv_intervals.py), and a symbolic
        # VCF deletion is anchored on the base BEFORE it, with END the last
        # deleted base. Writing POS=start put every deletion one base late and,
        # without END, a region query inside the deletion missed the record.
        st0, en0 = int(r["start"]), int(r.get("end") or r["start"])
        pos = max(st0 - 1, 1)
        info.append(f"END={en0}")
        recs.append((a.h37rv_contig, pos, r["interval"], h37_base(pos),
                     [f"<{r.get('svtype','DEL')}>"], info, cell))
        n_iv += 1
    if ivs:
        print(f"  {n_iv:,} deletion records written from the catalogue")

    for r in sv_rows:
        p = r.get("h37rv_pos") or ""
        if not p or p in (".", "NA"):
            continue
        svtype = r.get("svtype", "SV")
        # With the catalogue in play the caller's deletions are superseded by
        # it; only the insertions it cannot reach are still needed.
        if ivs and svtype == "DEL":
            continue
        info = [f"CLASS=sv", f"SVTYPE={svtype}",
                f"SVLEN={r.get('svlen','0')}",
                f"COMPONENT={r.get('component','')}",
                f"QUALBAND={r.get('qual_band','')}"]
        measured = sv_state.get(r["key"], {})
        cell = []
        for s in order:
            st = r.get(s, "")
            if st != "ALT" and measured.get(s):
                st = measured[s]
            cell.append(f'{STATE_GT.get(st, ".")}:{st or "NA"}')
        if svtype == "DEL" and str(r.get("svlen") or "").lstrip("-").isdigit():
            info.append(f"END={int(p) + abs(int(r['svlen']))}")
        recs.append((a.h37rv_contig, int(p), r["key"], h37_base(int(p)),
                     [f"<{svtype}>"], info, "\t".join(cell)))
        n_sv += 1

    # ---- LEVEL 1: does the sample carry the insert at all? ----------------
    # Michael's two-level framing. An accessory locus carries two different
    # characters and the pipeline had been emitting only the second: variation
    # INSIDE the insert, which for a sample that does not have the insert is
    # neither reference nor unknown but INAPPLICABLE. The first character --
    # whether the insert is there -- is stated for every sample and was not
    # being emitted at all.
    #
    # Polarity follows the coordinate frame. These loci are sequence present in
    # panel genomes and absent from H37Rv, so carrying the insert is the DERIVED
    # state and gets the ALT allele. A convergence test then looks for repeated
    # independent gains of the insert, which is the event of interest.
    n_l1 = 0
    if a.accessory_presence:
        cat = {r["locus_id"]: r for r in csv.DictReader(
            open(a.accessory_catalogue, newline=""), delimiter="\t")}
        pres = collections.defaultdict(dict)
        pf = sorted(glob.glob(os.path.join(a.accessory_presence,
                                           "*.presence.tsv")))
        for f in pf:
            with open(f, newline="") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    pres[r["locus"]][r["sample"]] = r["state"]
        got = {s2 for d in pres.values() for s2 in d}
        print(f"  level 1: {len(pf):,} presence tables, {len(pres):,} loci, "
              f"{len(got):,} samples ({len(got & set(order)):,} in this cohort)")
        L1 = {"PRESENT": "ALT", "ABSENT": "REF", "UNCERTAIN": "NOCALL"}
        n_st = collections.Counter()
        for lid, per in sorted(pres.items()):
            m = cat.get(lid)
            if not m or not (m.get("pos") or "").isdigit():
                continue
            st = {s2: L1.get(per.get(s2, ""), "NOCALL") for s2 in order}
            n_st.update(st.values())
            info = ["CLASS=accessory_presence", "COMPONENT=level1",
                    f"ACCLOCUS={lid}",
                    f"ACCLEN={m.get('graph_len') or m.get('rep_len') or 0}",
                    f"ACCKLASS={m.get('klass','')}",
                    f"PANELCARRIERS={m.get('graph_carriers') or m.get('carriers_any') or 0}"]
            if str(m.get("route", "")):
                info.append(f"ACCROUTE={m['route']}")
            cell = "\t".join(f'{STATE_GT[st[s2]]}:{st[s2]}' for s2 in order)
            recs.append((a.h37rv_contig, int(m["pos"]), f"acc:{lid}",
                         h37_base(int(m["pos"])), ["<INS>"], info, cell))
            n_l1 += 1
        t = sum(n_st.values()) or 1
        print(f"  {n_l1:,} level-1 records: "
              + "  ".join(f"{k} {v:,} ({v/t:.1%})"
                          for k, v in sorted(n_st.items())))

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
                f"MEINFO=IS6110,1,1355,+",
                f"EVIDENCE={m.get('evidence','')}",
                f"SITECLASS={m.get('site_class','')}",
                f"FRAME={m.get('frame','')}"]
        if m.get("frame") == "node":
            info.append(f"NODE={m.get('node','')}")
        cell = []
        for s in order:
            r = per.get(s)
            st = r["state"] if r else "NOCALL"
            cell.append(f'{STATE_GT.get(st, ".")}:{st}')
        ref_b = h37_base(pos) if chrom == a.h37rv_contig else "N"
        recs.append((chrom, pos, k, ref_b, ["<INS:ME:IS6110>"], info,
                     "\t".join(cell)))
        n_is += 1

    recs.sort(key=lambda x: (x[0] != a.h37rv_contig, x[0], x[1]))

    # ---- write -------------------------------------------------------------
    today = datetime.date.today().strftime("%Y%m%d")
    tmp = a.out[:-3] if a.out.endswith(".gz") else a.out
    with open(tmp, "w") as fh:
        w = fh.write
        w("##fileformat=VCFv4.2\n")
        w(f"##fileDate={today}\n")
        w("##source=bin/merge_cohort_vcf.py\n")
        w(f"##cohort={a.cohort_name}\n")
        w(f"##graphBuild={a.build_id}\n")
        w(f"##reference={a.h37rv_contig}\n")
        w(f"##contig=<ID={a.h37rv_contig},length={a.h37rv_length}>\n")
        for nid in sorted(nodes_used, key=lambda x: int(x) if x.isdigit() else 0):
            w(f"##contig=<ID=node_{nid},length={max(nodes_used[nid],1)}>\n")
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
        w('##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele, '
          'reconstructed by Fitch parsimony over the 333-genome CX333 tree '
          'rooted on the canettii outgroup GCF_035581225. A property of the '
          'panel rather than of this cohort. `.` where the root state is tied; '
          'see AA_FLAG">\n')
        w('##INFO=<ID=AA_FLAG,Number=1,Type=String,Description="TIED:<bases> '
          'where parsimony does not resolve the root state, NODATA where the '
          'site had none. 5,130 of 72,986 panel sites are tied and are NOT '
          'resolved by a coin toss">\n')
        w('##INFO=<ID=AA_INVERTED,Number=0,Type=Flag,Description="The REFERENCE '
          'carries the DERIVED allele at this site, so GT=0 is derived and GT=1 '
          'is ancestral. True of 6,525 of 72,986 panel SNP sites, 8.9%, which is '
          'why AA is annotated rather than REF being re-polarised -- VCF requires '
          'REF to match the reference base">\n')
        w('##INFO=<ID=IS6110PROX,Number=0,Type=Flag,Description="The interval '
          'lies within 500 bp of an IS6110 landmark, where depth counts '
          'mismapped copies of the element unless it is filtered">\n')
        w('##INFO=<ID=NODE,Number=1,Type=String,Description="Graph node, for '
          'node-frame records">\n')
        for t, d in (("KIND", "small-variant kind from P4"),
                     ("REGION", "core, pe_ppe, masked, off_path_near or off_path_accessory"),
                     ("SVTYPE", "structural variant type"),
                     ("SVLEN", "structural variant length"),
                     ("COMPONENT", "called or inherited; inherited rows carry no caller QUAL"),
                     ("QUALBAND", "caller quality band, or composed"),
                     ("MEINFO", "mobile element info: name,start,end,polarity"),
                     ("EVIDENCE", "IS6110 support: two_sided or one_sided"),
                     ("SVSOURCE", "graph or caller: where the interval came from"),
                     ("SVTIER", "interval support tier: A_multi_assembly, B_single_assembly or C_caller_only"),
                     ("PANELCARRIERS", "panel genomes carrying the deletion, from the graph"),
                     ("SITECLASS", "IS6110 relationship to the reference: ref_lacking or ref_shared"),
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
          'verbatim: ALT, REF, ABSENT or NOCALL. Kept so the GT mapping loses '
          'nothing">\n')
        w("#" + "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                           "FILTER", "INFO", "FORMAT"] + order) + "\n")
        n_star_added = 0
        for chrom, pos, ident, ref, alts, info, cells in recs:
            # GT=2 is the * allele (ABSENT). Every block maps ABSENT to 2, but
            # only the small-variant block used to add *, so the SV, interval,
            # accessory and IS6110 records wrote GT=2 against one ALT -- 976
            # records in scale200, and `bcftools +fill-tags` aborts on the
            # first. Decide it here, once, from the cells actually written.
            if "*" not in alts and ("\t" + cells).find("\t2:") != -1:
                alts = list(alts) + ["*"]
                n_star_added += 1
            w("\t".join([chrom, str(pos), ident, ref, ",".join(alts), ".",
                         "PASS", ";".join(info), "GT:ST"]) + "\t"
              + cells + "\n")

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
    print(f"  -> {a.out}")

if __name__ == "__main__":
    main()
