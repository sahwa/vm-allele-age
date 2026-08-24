#!/usr/bin/env python3
"""
Confirm the stale-VCF diagnosis.

Recovers, for ALL 50,000 .bed columns, which tree-sequence individual they are,
using a hash of each individual's genotype fingerprint (O(n), not pairwise).

Distinguishes:
  STALE DRAW  - bed individuals are a sorted subset, but a DIFFERENT one than
                sample_idx. Files came from different runs.
  PERMUTATION - bed individuals are the right set in the wrong order.
  Neither     - something else is going on.

Run:  python confirm_stale_draw.py
"""

import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import tskit

SIM_PATH = Path(
    "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.1"
)
SIM_VERSION = "1.1"

PREFIX = SIM_PATH / f"{SIM_VERSION}_neutral_out"
VCF_FILE = SIM_PATH / f"{SIM_VERSION}_neutral_out.vcf"
MUT_TREE_FILE = SIM_PATH / f"{SIM_VERSION}_neutral_out.recap.mut.trees"
PHENO_CSV = SIM_PATH / f"{SIM_VERSION}_phenotypes.csv"

N_FINGERPRINT = 600
FREQ_LO, FREQ_HI = 0.05, 0.45


# ------------------------------------------------------------ mtimes first
print("=" * 66)
print("0. FILE TIMESTAMPS")
print("=" * 66)
for p in [MUT_TREE_FILE, VCF_FILE, Path(f"{PREFIX}.bed"), Path(f"{PREFIX}.bim"),
          Path(f"{PREFIX}.fam"), PHENO_CSV,
          SIM_PATH / f"{SIM_VERSION}_variant_info.csv"]:
    if p.exists():
        t = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(p.stat().st_mtime))
        print(f"  {t}   {p.name}")
    else:
        print(f"  {'MISSING':<19}   {p.name}")
print("\nA VCF older than phenotypes.csv by more than one runtime = stale.\n")


# ------------------------------------------------------------ bed decoding
_CODE = {0: 2, 1: -9, 2: 1, 3: 0}
_LUT = np.empty((256, 4), dtype=np.int8)
for _b in range(256):
    for _j in range(4):
        _LUT[_b, _j] = _CODE[(_b >> (2 * _j)) & 3]


def read_variants(bed_path, var_indices, n_ind):
    bpv = (n_ind + 3) // 4
    out = np.empty((len(var_indices), n_ind), dtype=np.int8)
    with open(bed_path, "rb") as f:
        for k, vi in enumerate(var_indices):
            f.seek(3 + vi * bpv)
            raw = np.frombuffer(f.read(bpv), dtype=np.uint8)
            out[k] = _LUT[raw].reshape(-1)[:n_ind]
    return out


# ------------------------------------------------------------------- load
bim = pd.read_csv(f"{PREFIX}.bim", sep=r"\s+", header=None,
                  names=["chr", "snpid", "cm", "pos", "a1", "a2"])
fam = pd.read_csv(f"{PREFIX}.fam", sep=r"\s+", header=None,
                  names=["fid", "iid", "pat", "mat", "sex", "phe"])
n_ind_bed = len(fam)

sample_idx = pd.read_csv(PHENO_CSV)["sample_idx"].to_numpy()

print("Loading tree sequence...", flush=True)
mts = tskit.load(MUT_TREE_FILE)
n_dip_all = mts.num_samples // 2

pos_to_bimrow = {p: i for i, p in enumerate(bim["pos"].to_numpy())}
stride = max(1, mts.num_sites // (N_FINGERPRINT * 6))

fp_rows, fp_cols = [], []
for si, var in enumerate(mts.variants()):
    if len(fp_rows) >= N_FINGERPRINT:
        break
    if si % stride:
        continue
    row = pos_to_bimrow.get(int(var.site.position))
    if row is None:
        continue
    gt = var.genotypes
    if not (FREQ_LO < gt.mean() < FREQ_HI):
        continue
    fp_rows.append(row)
    fp_cols.append(gt.reshape(-1, 2).sum(axis=1))

K = len(fp_rows)
print(f"fingerprinting on {K} variants\n")

G_pop = np.ascontiguousarray(np.asarray(fp_cols, dtype=np.int8).T)   # (n_dip_all, K)
G_bed = np.ascontiguousarray(read_variants(Path(f"{PREFIX}.bed"), fp_rows,
                                           n_ind_bed).T)             # (n_bed, K)


# ------------------------------------------------- hash-based full mapping
print("=" * 66)
print("1. RECOVERING FULL MAPPING")
print("=" * 66)

pop_hash = {}
collisions = 0
for i in range(n_dip_all):
    h = G_pop[i].tobytes()
    if h in pop_hash:
        collisions += 1
    pop_hash[h] = i
print(f"fingerprint collisions in population: {collisions} "
      f"({'ok' if collisions < n_dip_all * 0.001 else 'TOO MANY - raise K'})")

mapping = np.full(n_ind_bed, -1, dtype=np.int64)
for j in range(n_ind_bed):
    mapping[j] = pop_hash.get(G_bed[j].tobytes(), -1)

n_unmatched = int((mapping == -1).sum())
print(f"bed columns matched to a population individual: "
      f"{n_ind_bed - n_unmatched}/{n_ind_bed}")
if n_unmatched:
    print("   -> unmatched columns mean the .bed individuals are not in this")
    print("      tree sequence at all. Check that the .trees file matches too.")


# ------------------------------------------------------------- interpret
print("\n" + "=" * 66)
print("2. INTERPRETATION")
print("=" * 66)

matched = mapping[mapping >= 0]
is_sorted = bool(np.all(np.diff(matched) > 0))
print(f"bed individuals are strictly increasing: {is_sorted}")

set_bed = set(matched.tolist())
set_pheno = set(sample_idx.tolist())
overlap = len(set_bed & set_pheno)
print(f"set overlap with phenotypes.csv sample_idx: {overlap}/{len(set_pheno)} "
      f"({100 * overlap / len(set_pheno):.1f}%)")
print(f"expected overlap for two independent draws of {len(set_pheno)} "
      f"from {n_dip_all}: ~{100 * len(set_pheno) / n_dip_all:.1f}%")

exact_order = np.array_equal(matched, np.sort(sample_idx))
print(f"identical to sample_idx: {exact_order}")

print()
if is_sorted and overlap < len(set_pheno) * 0.9 and not exact_order:
    print("VERDICT: STALE DRAW.")
    print("  The .bed holds a valid, correctly-ordered set of individuals - just")
    print("  not the ones the phenotypes were computed for. Your indexing logic")
    print("  is fine; the files are from different executions.")
    print("  Fix: re-run the pipeline end-to-end into a clean directory.")
elif not is_sorted and overlap > len(set_pheno) * 0.99:
    print("VERDICT: PERMUTATION.")
    print("  Right people, wrong order - a real indexing bug in the VCF write.")
else:
    print("VERDICT: inconclusive - read the numbers above directly.")

# Save the recovered mapping so the existing files can be rescued if desired.
out = SIM_PATH / f"{SIM_VERSION}_recovered_individual_mapping.csv"
pd.DataFrame({"bed_column": np.arange(n_ind_bed),
              "iid": fam["iid"].to_numpy(),
              "true_pair_index": mapping}).to_csv(out, index=False)
print(f"\nRecovered mapping written to {out.name}")
print("If the .trees file is intact you can rebuild phenotypes against this")
print("mapping instead of regenerating the genotypes.")