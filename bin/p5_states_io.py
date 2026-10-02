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


def load_states(state_dir, keys, samples, quiet=False):
    """(len(keys) x len(samples)) uint8 array of CODE values."""
    ksha = keys_sha1(keys)
    kindex = {k["key"]: i for i, k in enumerate(keys)}
    M = np.full((len(keys), len(samples)), CODE["NOCALL"], dtype=np.uint8)
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
            else:
                fh.seek(0)
                for x in csv.DictReader(fh, delimiter="\t"):
                    i = kindex.get(x["key"])
                    if i is not None:
                        M[i, j] = CODE.get(x["state"], CODE["NOCALL"])
    if n_sparse and not quiet:
        print(f"  {n_sparse} of {len(samples)} state files in the sparse form")
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


def samples_in(state_dir, order):
    return [s for s in order
            if os.path.exists(os.path.join(state_dir, f"{s}.states.tsv"))]
