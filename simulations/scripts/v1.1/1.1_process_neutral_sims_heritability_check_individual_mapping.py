#!/usr/bin/env python3
"""
Locate the individual-indexing bug.

Variant alignment is already confirmed perfect, so the scrambling is in WHICH
person each VCF column corresponds to. This script fingerprints individuals on
a few hundred variants and finds, empirically, which tree-sequence individual
each .bed column actually is.

Tests the two competing conventions:
  PAIR : g_all index i  ==  mts.samples()[2i], mts.samples()[2i+1]   (step 8)
  IND  : individual table ID i                                        (write_vcf)

Run:  python find_individual_mapping.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import tskit

SIM_PATH = Path(
    "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.1"
)
SIM_VERSION = "1.1"

PREFIX = SIM_PATH / f"{SIM_VERSION}_neutral_out"
MUT_TREE_FILE = SIM_PATH / f"{SIM_VERSION}_neutral_out.recap.mut.trees"
PHENO_CSV = SIM_PATH / f"{SIM_VERSION}_phenotypes.csv"

N_FINGERPRINT = 600      # variants used to fingerprint each individual
N_TEST = 12              # .bed individuals to trace back
FREQ_LO, FREQ_HI = 0.05, 0.45

ok = lambda b: "PASS" if b else "**FAIL**"


# ------------------------------------------------------ .bed single-variant read
_CODE = {0: 2, 1: -9, 2: 1, 3: 0}
_LUT = np.empty((256, 4), dtype=np.int8)
for _b in range(256):
    for _j in range(4):
        _LUT[_b, _j] = _CODE[(_b >> (2 * _j)) & 3]


def read_variants(bed_path, var_indices, n_ind):
    """Read specific variant rows by seeking. Returns (n_var, n_ind) int8."""
    bpv = (n_ind + 3) // 4
    out = np.empty((len(var_indices), n_ind), dtype=np.int8)
    with open(bed_path, "rb") as f:
        for k, vi in enumerate(var_indices):
            f.seek(3 + vi * bpv)
            raw = np.frombuffer(f.read(bpv), dtype=np.uint8)
            out[k] = _LUT[raw].reshape(-1)[:n_ind]
    return out


# ---------------------------------------------------------------------- load
bim = pd.read_csv(f"{PREFIX}.bim", sep=r"\s+", header=None,
                  names=["chr", "snpid", "cm", "pos", "a1", "a2"])
fam = pd.read_csv(f"{PREFIX}.fam", sep=r"\s+", header=None,
                  names=["fid", "iid", "pat", "mat", "sex", "phe"])
n_ind_bed, n_var_bed = len(fam), len(bim)

pheno = pd.read_csv(PHENO_CSV)
sample_idx = pheno["sample_idx"].to_numpy()

print("Loading tree sequence...", flush=True)
mts = tskit.load(MUT_TREE_FILE)
n_dip_all = mts.num_samples // 2


# ================================================ 1. the mapping, stated plainly
print("\n" + "=" * 64)
print("1. NODE -> INDIVIDUAL STRUCTURE")
print("=" * 64)

samples = mts.samples()
ind_of_node = np.array([mts.node(n).individual for n in samples])

print(f"num_samples        {mts.num_samples}")
print(f"num_samples // 2   {n_dip_all}")
print(f"num_individuals    {mts.num_individuals}")
print(f"individuals with sample nodes: {len(np.unique(ind_of_node))}")
extra = mts.num_individuals - n_dip_all
print(f"extra individuals in table: {extra}   {ok(extra == 0)}")
if extra:
    print("   -> SLiM retained non-sample individuals. Individual table IDs are")
    print("      NOT the same thing as the step-8 pair index.")

print(f"any sample node with individual == -1: "
      f"{bool((ind_of_node == -1).any())}   {ok(not (ind_of_node == -1).any())}")

pairs_consistent = np.array_equal(ind_of_node[::2], ind_of_node[1::2])
print(f"adjacent node pairs share an individual: {pairs_consistent}   "
      f"{ok(pairs_consistent)}")

sample_ind_ids = np.unique(ind_of_node)
identity = np.array_equal(ind_of_node[::2], np.arange(n_dip_all))
print(f"pair index == individual ID (identity): {identity}   {ok(identity)}")
print(f"sample individual ID range: {sample_ind_ids.min()} .. {sample_ind_ids.max()}")
contiguous = np.array_equal(sample_ind_ids,
                            np.arange(sample_ind_ids.min(), sample_ind_ids.max() + 1))
print(f"sample individual IDs contiguous: {contiguous}")
if not identity and contiguous:
    print(f"   -> constant offset of {sample_ind_ids.min()}")


# ==================================== 2. pick fingerprint variants & pull dosages
print("\n" + "=" * 64)
print("2. BUILDING FINGERPRINTS")
print("=" * 64)

pos_to_bimrow = {p: i for i, p in enumerate(bim["pos"].to_numpy())}
stride = max(1, mts.num_sites // (N_FINGERPRINT * 6))

fp_bimrows, fp_pair, fp_ind = [], [], []

# argsort groups nodes by individual; stable keeps within-individual node order
order = np.argsort(ind_of_node, kind="stable")
if len(order) != 2 * len(sample_ind_ids):
    raise SystemExit("Individuals do not all carry exactly 2 sample nodes.")
ind_node_pairs = order.reshape(len(sample_ind_ids), 2)

for si, var in enumerate(mts.variants()):
    if len(fp_bimrows) >= N_FINGERPRINT:
        break
    if si % stride:
        continue
    p = int(var.site.position)
    row = pos_to_bimrow.get(p)
    if row is None:
        continue
    gt = var.genotypes
    f = gt.mean()
    if not (FREQ_LO < f < FREQ_HI):
        continue

    fp_bimrows.append(row)
    fp_pair.append(gt.reshape(-1, 2).sum(axis=1))          # PAIR convention
    fp_ind.append(gt[ind_node_pairs].sum(axis=1))          # IND convention

K = len(fp_bimrows)
print(f"using {K} fingerprint variants")
if K < 50:
    raise SystemExit("Too few fingerprint variants; widen the frequency window.")

G_pair = np.asarray(fp_pair, dtype=np.int8)   # (K, n_dip_all), step-8 ordering
G_ind = np.asarray(fp_ind, dtype=np.int8)     # (K, n_sample_individuals), by ID
G_bed = read_variants(Path(f"{PREFIX}.bed"), fp_bimrows, n_ind_bed)

print(f"G_bed {G_bed.shape}   G_pair {G_pair.shape}   G_ind {G_ind.shape}")


# ======================================================== 3. targeted check
print("\n" + "=" * 64)
print("3. DOES .bed COLUMN j MATCH sample_idx[j] UNDER EITHER CONVENTION?")
print("=" * 64)

test_j = np.linspace(0, n_ind_bed - 1, N_TEST).astype(int)
hit_pair = hit_ind = 0
for j in test_j:
    v = G_bed[:, j]
    mm_pair = int((G_pair[:, sample_idx[j]] != v).sum())
    id_pos = np.searchsorted(sample_ind_ids, sample_idx[j])
    mm_ind = (int((G_ind[:, id_pos] != v).sum())
              if id_pos < len(sample_ind_ids)
              and sample_ind_ids[id_pos] == sample_idx[j] else -1)
    hit_pair += mm_pair == 0
    hit_ind += mm_ind == 0
    print(f"  bed col {j:>6} (sample_idx {sample_idx[j]:>6})  "
          f"mismatches PAIR {mm_pair:>4}   IND {mm_ind:>4}   of {K}")

print(f"\nPAIR convention exact hits: {hit_pair}/{N_TEST}")
print(f"IND  convention exact hits: {hit_ind}/{N_TEST}")


# ============================== 4. if neither, search for who they really are
if hit_pair < N_TEST and hit_ind < N_TEST:
    print("\n" + "=" * 64)
    print("4. UNTARGETED SEARCH - who IS each .bed column?")
    print("=" * 64)
    print("(scanning all tree-sequence individuals; best and runner-up shown)\n")

    for j in test_j[:6]:
        v = G_bed[:, j][:, None]
        mism = (G_pair != v).sum(axis=0)
        best = int(np.argmin(mism))
        srt = np.partition(mism, 1)[:2]
        print(f"  bed col {j:>6}: best PAIR index {best:>6} "
              f"({int(mism[best])} mismatches, runner-up {int(srt.max())})   "
              f"expected {sample_idx[j]}   delta {best - sample_idx[j]:+d}")

    print("\nA constant non-zero delta means an offset; a scattered delta means")
    print("a genuine permutation. Either way, step 8 and step 10 disagree.")

print("\nDone.")