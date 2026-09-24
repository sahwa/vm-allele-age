#!/usr/bin/env python
"""Reference cell: derive SIGMA2_B and SIGMA2_E for the diagGREML sweep.

Anchors per-mutation effect variance so that singleton h2 = H2_REF at
UKB scale (n = N_REF diploids, L = L_REF bp, no selection). Both constants
are then held fixed across every (n, L, S) cell of the sweep.
"""

import json
from pathlib import Path

import numpy as np
import tskit

# --- anchor definition ---
H2_REF = 0.10          # singleton h2 at UKB scale
L_REF  = 3e9           # UKB WGS sequence length

TREE_FILE = Path("/exafs1/well/visscher-wray/users/uwu199/projects/"
                 "vm-allele-age/simulations/data/v1.2/n500000_L1e+08/sim.trees")
OUT_JSON  = TREE_FILE.parent / "reference_constants.json"

RNG_SEED, REP, S = 42, 0, 0.0          # reference must have S = 0

_ss = np.random.SeedSequence([RNG_SEED, REP])
seed_beta, seed_noise = _ss.spawn(2)
rng_beta  = np.random.default_rng(seed_beta)
rng_noise = np.random.default_rng(seed_noise)

# --- load; no multiallelic stripping needed ---
ts = tskit.load(str(TREE_FILE))
N_IND = ts.num_samples // 2
print(f"{ts.num_sites:,} sites, {N_IND:,} diploids, "
      f"L = {ts.sequence_length:.0e}", flush=True)

tabs = ts.tables
mut_node  = tabs.mutations.node
is_sample = np.zeros(ts.num_nodes, dtype=bool)
is_sample[ts.samples()] = True
sing = is_sample[mut_node]

carrier   = tabs.nodes.individual[mut_node[sing]]
mt        = tabs.mutations.time[sing]
sing_ages = np.where(np.isnan(mt), tabs.nodes.time[mut_node[sing]], mt)
M_SING    = int(sing.sum())

assert (carrier >= 0).all(), "sample nodes not mapped to individuals"
print(f"{M_SING:,} singletons in this run", flush=True)

# --- per-mutation variance, anchored at UKB scale ---
M_REF = M_SING * (L_REF / ts.sequence_length)
SIGMA2_B_REF = H2_REF * N_IND / M_REF
print(f"M_REF (extrapolated to L={L_REF:.0e}): {M_REF:,.0f}")
print(f"SIGMA2_B_REF = {SIGMA2_B_REF:.6g}", flush=True)

# --- residual variance, fixed from the same anchor ---
SIGMA2_E = H2_REF * (1 - H2_REF) / H2_REF      # = V_A_ref * (1-h2)/h2, V_A_ref = H2_REF
print(f"SIGMA2_E = {SIGMA2_E:.6g}", flush=True)

# --- sanity check in this cell ---
sd = np.sqrt(SIGMA2_B_REF * np.exp(-S * sing_ages))
beta_sing = rng_beta.normal(0.0, sd)

g = np.bincount(carrier, weights=beta_sing, minlength=N_IND)
g -= g.mean()
y = g + rng_noise.normal(0, np.sqrt(SIGMA2_E), size=N_IND)

h2_here = g.var() / (g.var() + SIGMA2_E)
print(f"\nthis cell:  V_A = {g.var():.6f}   realised h2 = {h2_here:.5f}")
print(f"expected  :  V_A ~ {H2_REF * ts.sequence_length / L_REF:.6f}", flush=True)

json.dump({"SIGMA2_B_REF": SIGMA2_B_REF, "SIGMA2_E": SIGMA2_E,
           "H2_REF": H2_REF, "L_REF": L_REF,
           "N_IND_ref": N_IND, "M_SING_this_run": M_SING},
          open(OUT_JSON, "w"), indent=2)
print(f"\nWrote {OUT_JSON}", flush=True)
