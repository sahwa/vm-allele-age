#!/usr/bin/env python
"""
Singleton phenotype generation for the diagGREML sweep.

For one (n, L) cell, sweeps over selection strength S and writes, per S:
  - per-individual phenotype
  - per-individual singleton counts by age bin
  - true per-bin genetic variance

Effects are drawn per singleton with variance SIGMA2_B_REF * exp(-S * age).
SIGMA2_B_REF and SIGMA2_E are fixed from a reference cell so realised h2
varies across the sweep rather than being pinned.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import tskit

# =================================================================
# 0. Configuration
# =================================================================
p = argparse.ArgumentParser()
p.add_argument("--tree", required=True, help="msprime .trees for this cell")
p.add_argument("--ref-json", required=True, help="reference_constants.json")
p.add_argument("--out-dir", required=True)
p.add_argument("--rep", type=int, default=0)
p.add_argument("--seed", type=int, default=42)
p.add_argument("--mu", type=float, default=1.25e-8)
args = p.parse_args()

S_VALUES = [0.0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2]
BINS = np.concatenate([[0.0], np.geomspace(1, 205, 7), [np.inf]])

OUT = Path(args.out_dir)
OUT.mkdir(parents=True, exist_ok=True)

with open(args.ref_json) as fh:
    REF = json.load(fh)
SIGMA2_B_REF = REF["SIGMA2_B_REF"]
SIGMA2_E     = REF["SIGMA2_E"]
print(f"Reference: SIGMA2_B={SIGMA2_B_REF:.6g}  SIGMA2_E={SIGMA2_E:.4g}",
      flush=True)

# =================================================================
# 1. Load and extract singletons from the tables
# =================================================================
ts = tskit.load(args.tree)
N_IND = ts.num_samples // 2
print(f"{ts.num_sites:,} sites, {N_IND:,} diploids, "
      f"L = {ts.sequence_length:.0e}", flush=True)

tabs = ts.tables
mut_node = tabs.mutations.node

is_sample = np.zeros(ts.num_nodes, dtype=bool)
is_sample[ts.samples()] = True
sing = is_sample[mut_node]

carrier   = tabs.nodes.individual[mut_node[sing]]
_mt       = tabs.mutations.time[sing]
sing_ages = np.where(np.isnan(_mt), tabs.nodes.time[mut_node[sing]], _mt)
M_SING    = int(sing.sum())

assert (carrier >= 0).all(), "sample nodes not mapped to individuals"
print(f"{M_SING:,} singletons  (mean {M_SING / N_IND:.1f} per person)", flush=True)

P_SING  = 1.0 / ts.num_samples          # allele frequency of a singleton
TWO_PQ  = 2 * P_SING * (1 - P_SING)
N_BINS  = len(BINS) - 1

BIN_LABELS = [f"{lo:.0f}+" if np.isinf(hi) else f"{lo:.0f}-{hi:.0f}"
              for lo, hi in zip(BINS[:-1], BINS[1:])]
bin_idx = np.searchsorted(BINS, sing_ages, side="right") - 1

# per-individual singleton counts by bin (independent of S)
sing_counts = np.bincount(
    carrier * N_BINS + bin_idx, minlength=N_IND * N_BINS
).reshape(N_IND, N_BINS)
assert sing_counts.sum() == M_SING

pd.DataFrame(sing_counts, columns=BIN_LABELS).to_csv(
    OUT / "singleton_counts.csv", index=False)

# V_M: mutations per generation x per-mutation effect variance
V_M_TRUE = 2 * args.mu * ts.sequence_length * SIGMA2_B_REF

iids = [f"ind{i}" for i in range(N_IND)]

# =================================================================
# 2. Sweep over selection strength
# =================================================================
summary = []

for S in S_VALUES:
    tag = f"S{S:.0e}"
    rng = np.random.default_rng(
        np.random.SeedSequence([args.seed, args.rep, int(S * 1e6)]))
    rng_beta, rng_noise = (np.random.default_rng(s) for s in rng.spawn(2))

    sd = np.sqrt(SIGMA2_B_REF * np.exp(-S * sing_ages))
    beta_sing = rng_beta.normal(0.0, sd)

    g = np.bincount(carrier, weights=beta_sing, minlength=N_IND)
    g -= g.mean()
    y = g + rng_noise.normal(0, np.sqrt(SIGMA2_E), size=N_IND)

    var_g = g.var()
    h2    = var_g / y.var()
    V_A_emp = float(np.sum(TWO_PQ * beta_sing**2))

    print(f"\n=== S = {S:.0e} ===")
    print(f"var(g) = {var_g:.5g}   V_E = {SIGMA2_E:.4g}   realised h2 = {h2:.4f}")
    print(f"V_A (sum 2pq*beta^2) = {V_A_emp:.5g}")
    print(f"V_M (analytic)       = {V_M_TRUE:.5g}")
    print(f"implied T = V_A/V_M  = {V_A_emp / V_M_TRUE:.4g} generations", flush=True)

    # ---- true per-bin genetic variance ----
    rows = []
    print(f"\n{'Bin (gens)':<16}{'n_sing':>10}{'mean_beta2':>13}"
          f"{'V_observed':>13}{'share':>9}")
    for b, label in enumerate(BIN_LABELS):
        m = bin_idx == b
        if not m.any():
            continue
        contribs = TWO_PQ * beta_sing[m] ** 2
        V_bin = contribs.sum()
        ss = (contribs ** 2).sum()
        n_eff = V_bin ** 2 / ss if ss > 0 else 0.0
        mean_b2 = (beta_sing[m] ** 2).mean()

        print(f"{label:<16}{int(m.sum()):>10,}{mean_b2:>13.4g}"
              f"{V_bin:>13.4g}{V_bin / V_A_emp:>9.3f}")
        rows.append({
            "S": S, "bin": label, "bin_lo": BINS[b], "bin_hi": BINS[b + 1],
            "n_variants": int(m.sum()), "n_eff": n_eff,
            "mean_beta2": mean_b2, "V_observed": V_bin,
            "true_share": V_bin / V_A_emp,
        })

    pd.DataFrame(rows).to_csv(OUT / f"bin_truth_{tag}.csv", index=False)

    # ---- phenotypes ----
    pd.DataFrame({"iid": iids, "y": y, "g": g}).to_csv(
        OUT / f"phenotypes_{tag}.csv", index=False)
    pd.DataFrame({"FID": 0, "IID": iids, "PHENO": y}).to_csv(
        OUT / f"phenotypes_{tag}.GENIE.txt", sep="\t", index=False)

    summary.append({
        "rep": args.rep, "S": S, "n_ind": N_IND, "L": ts.sequence_length,
        "M_sing": M_SING, "mean_c": M_SING / N_IND,
        "var_g": var_g, "V_E": SIGMA2_E, "h2_realised": h2,
        "V_A_empirical": V_A_emp, "V_M_true": V_M_TRUE,
        "T_implied": V_A_emp / V_M_TRUE,
    })

pd.DataFrame(summary).to_csv(OUT / "sweep_summary.csv", index=False)
print(f"\nWrote outputs to {OUT}", flush=True)
