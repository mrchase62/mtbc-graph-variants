#!/usr/bin/env python3
"""Read P5's per-sample state files, sparse or dense, into one uint8 array.

One reader, used by every consumer, so the sparse format is defined in a single
place. The array is (keys x samples) of the four codes below; at 10,000
isolates and 540,000 keys that is 5.4 GB, which fits, while the dense TSV of
the same information is 26.5 GB of text, which is the thing not to build.

THE KEY SET IS THE INDEX in the sparse form, so it is checked rather than
assumed: the file header carries the key count and a SHA-1 of the key list, and
a mismatch is fatal. A sparse file read against the wrong key set is not
slightly wrong, it is wrong at every position, and this project has twice paid
for a stale intermediate that looked complete.
"""
import csv, hashlib, os, sys
import numpy as np

CODE = {"ALT": 0, "REF": 1, "ABSENT": 2, "NOCALL": 3}
NAME = ["ALT", "REF", "ABSENT", "NOCALL"]


def keys_sha1(keys):
    return hashlib.sha1("\n".join(k["key"] for k in keys).encode()).hexdigest()


class EvidenceCodes:
    """The evidence vocabulary of one states array: code 0 is "no evidence
    recorded", every other code a reason string from p5_states.py (masked,
    junction_joined, no_projection, ...). Built as the files are read, so a
    new reason needs no table here; written beside the array so a reader
    never depends on the order codes were met in."""

    def __init__(self, names=None):
        self.names = list(names) if names else [""]
        self.index = {n: i for i, n in enumerate(self.names)}

    def code(self, name):
        if not name:
            return 0
        i = self.index.get(name)
        if i is None:
            if len(self.names) >= 255:
                sys.exit(f"FATAL: more than 254 distinct evidence codes; "
                         f"{name!r} would not fit a uint8")
            i = self.index[name] = len(self.names)
            self.names.append(name)
        return i


def load_states(state_dir, keys, samples, quiet=False, with_evidence=False):
    """(len(keys) x len(samples)) uint8 array of CODE values.

    With `with_evidence`, (M, EV, codes): EV is the same shape, holding each
    cell's evidence code (the sparse files' fourth column; 0 where none, and
    everywhere for a file written before the column existed), and `codes`
    the EvidenceCodes that name them."""
    ksha = keys_sha1(keys)
    kindex = {k["key"]: i for i, k in enumerate(keys)}
    M = np.full((len(keys), len(samples)), CODE["NOCALL"], dtype=np.uint8)
    EV = (np.zeros((len(keys), len(samples)), dtype=np.uint8)
          if with_evidence else None)
    codes = EvidenceCodes()
    n_sparse = 0
    for j, s in enumerate(samples):
        path = os.path.join(state_dir, f"{s}.states.tsv")
        with open(path, newline="") as fh:
            first = fh.readline()
            if first.startswith("#format\tsparse"):
                n_sparse += 1
                hdr, line = {}, fh.readline()
                while line.startswith("#"):
                    kk, _, vv = line[1:].rstrip("\n").partition("\t")
                    hdr[kk] = vv
                    line = fh.readline()
                if hdr.get("keys_sha1") and hdr["keys_sha1"] != ksha:
                    sys.exit(f"FATAL: {path} was written against a different "
                             f"key set ({hdr['keys_sha1'][:12]} vs "
                             f"{ksha[:12]}); re-run p5 --states")
                if hdr.get("n_keys") and int(hdr["n_keys"]) != len(keys):
                    sys.exit(f"FATAL: {path} has {hdr['n_keys']} keys, this "
                             f"key set has {len(keys)}")
                M[:, j] = CODE.get(hdr.get("default", "REF"), CODE["REF"])
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) >= 2 and f[0].isdigit():
                        M[int(f[0]), j] = CODE.get(f[1], CODE["NOCALL"])
                        if EV is not None and len(f) >= 4 and f[3]:
                            EV[int(f[0]), j] = codes.code(f[3])
            else:
                fh.seek(0)
                for x in csv.DictReader(fh, delimiter="\t"):
                    i = kindex.get(x["key"])
                    if i is not None:
                        M[i, j] = CODE.get(x["state"], CODE["NOCALL"])
                        if EV is not None and x.get("evidence"):
                            EV[i, j] = codes.code(x["evidence"])
    if n_sparse and not quiet:
        print(f"  {n_sparse} of {len(samples)} state files in the sparse form")
    if with_evidence:
        return M, EV, codes
    return M


# ---- the cohort states array ------------------------------------------------
#
# Written ONCE by p5_matrix.py, read by every cohort-level consumer. Each of
# them used to rebuild the same information from the dense text matrix, or from
# the 10,000 sparse files, which at 10,000 isolates is hundreds of millions of
# lines parsed per consumer. Saved as an uncompressed .npy it is memory-mapped:
# a consumer that needs a slice of keys -- one shard of the merged VCF -- touches
# only those rows, and nothing is parsed at all.
#
# Row i is keys.tsv row i; column j is the j-th sample in the meta file. The
# meta file carries the key count and the key-list checksum, and open_array
# refuses a mismatch for the same reason the sparse files do.
ARRAY_NAME = "states.u8.npy"
ARRAY_META = "states.meta.tsv"


def save_array(out_dir, M, samples, keys):
    """Write M (keys x samples, uint8 CODE values) and its meta file."""
    tmp = os.path.join(out_dir, "." + ARRAY_NAME + ".tmp")
    with open(tmp, "wb") as fh:
        np.save(fh, np.ascontiguousarray(M, dtype=np.uint8))
    os.replace(tmp, os.path.join(out_dir, ARRAY_NAME))
    tmp = os.path.join(out_dir, "." + ARRAY_META + ".tmp")
    with open(tmp, "w") as fh:
        fh.write(f"#n_keys\t{len(keys)}\n#keys_sha1\t{keys_sha1(keys)}\n"
                 f"#n_samples\t{len(samples)}\n")
        for s in samples:
            fh.write(s + "\n")
    os.replace(tmp, os.path.join(out_dir, ARRAY_META))


def open_array(array_dir, keys):
    """(memory-mapped keys x samples uint8 array, samples), checked against keys."""
    meta, samples = {}, []
    with open(os.path.join(array_dir, ARRAY_META)) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                k, _, v = line[1:].partition("\t")
                meta[k] = v
            elif line:
                samples.append(line)
    if meta.get("keys_sha1") != keys_sha1(keys) or int(meta.get("n_keys", -1)) != len(keys):
        sys.exit(f"FATAL: {array_dir}/{ARRAY_NAME} was written against a "
                 f"different key set; re-run p5 --matrix")
    M = np.load(os.path.join(array_dir, ARRAY_NAME), mmap_mode="r")
    if M.shape != (len(keys), len(samples)):
        sys.exit(f"FATAL: {ARRAY_NAME} is {M.shape}, expected "
                 f"({len(keys)}, {len(samples)})")
    return M, samples


# ---- the evidence array, beside the states array -----------------------------
# Same shape and row order as states.u8.npy, one code per cell (0: none), and
# the code table. Written by p5_matrix.py with the states array; read by
# merge_cohort_vcf.py for FORMAT/EV (HANDOFF ABS-1 item 4: each ABSENT/NOCALL
# cell's evidence in the merged VCF). Its meta is the states array's: the two
# are written together and checked against the same key set.
EV_NAME = "states.ev.u8.npy"
EV_CODES = "states.ev.codes.tsv"


def save_evidence(out_dir, EV, codes):
    tmp = os.path.join(out_dir, "." + EV_NAME + ".tmp")
    with open(tmp, "wb") as fh:
        np.save(fh, np.ascontiguousarray(EV, dtype=np.uint8))
    os.replace(tmp, os.path.join(out_dir, EV_NAME))
    tmp = os.path.join(out_dir, "." + EV_CODES + ".tmp")
    with open(tmp, "w") as fh:
        fh.write("code\tevidence\n")
        for i, n in enumerate(codes.names):
            fh.write(f"{i}\t{n}\n")
    os.replace(tmp, os.path.join(out_dir, EV_CODES))


def open_evidence(array_dir, shape):
    """(memory-mapped evidence array, code names) or (None, None) when the
    cohort was finished before the array existed. A present array of another
    shape is FATAL: it belongs to another generation of the states."""
    p = os.path.join(array_dir, EV_NAME)
    if not os.path.exists(p):
        return None, None
    st = os.path.join(array_dir, ARRAY_NAME)
    if os.path.exists(st) and os.path.getmtime(p) < os.path.getmtime(st):
        sys.exit(f"FATAL: {p} is older than {st}; it is another generation's "
                 f"evidence. Re-run p5 --matrix")
    names = [""]
    with open(os.path.join(array_dir, EV_CODES)) as fh:
        next(fh, None)
        for line in fh:
            i, _, n = line.rstrip("\n").partition("\t")
            if int(i) == len(names):
                names.append(n)
            elif int(i) != 0:
                sys.exit(f"FATAL: {array_dir}/{EV_CODES} is not numbered 0..n")
    EV = np.load(p, mmap_mode="r")
    if EV.shape != tuple(shape):
        sys.exit(f"FATAL: {EV_NAME} is {EV.shape}, the states array "
                 f"{tuple(shape)}; re-run p5 --matrix")
    return EV, names


def samples_in(state_dir, order):
    return [s for s in order
            if os.path.exists(os.path.join(state_dir, f"{s}.states.tsv"))]
