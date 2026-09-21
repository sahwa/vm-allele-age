#!/usr/bin/env python
"""Production singleton-count runs at UKB-matched ascertainment.

Usage: python run_production.py <n_dip> [L]

N_final per sample size is calibrated so the singleton fraction matches UKB
(~0.46). The 500k value is extrapolated from N_final ~ n^0.699 fitted to the
four calibrated points.
"""

import copy
import sys
from pathlib import Path

import daiquiri
import msprime
import numpy as np
import pandas as pd
import stdpopsim

# =================================================================
# Configuration
# =================================================================
SIM_VERSION = "1.2"
SIM_PATH = Path(
    "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data"
) / f"v{SIM_VERSION}"

POP = "EUR"
DEMOG_MODEL = "OutOfAfrica_2T12"
EUR_GROWTH_RATE = 0.0195
MU = 1.25e-8
REC = 1e-8
SEED = 1
YOUNG_MAX = 205.0                  # 2T12 growth onset
BIN_SCHEMES = [1, 2, 3, 4, 6]

N_FINAL = {
    20_000: 563_593,
    50_000: 941_197,
    100_000: 1_583_484,
    200_000: 2_822_145,
    500_000: 5_350_000,            # extrapolated
}

if len(sys.argv) < 2:
    sys.exit(f"usage: {sys.argv[0]} <n_dip> [L]\n"
             f"n_dip must be one of {sorted(N_FINAL)}")

N_DIP = int(sys.argv[1])
if N_DIP not in N_FINAL:
    sys.exit(f"no calibrated N_final for n={N_DIP}; "
             f"have {sorted(N_FINAL)}")

L = float(sys.argv[2]) if len(sys.argv) > 2 else 1e8

OUT = SIM_PATH / f"n{N_DIP}_L{L:.0e}"
OUT.mkdir(parents=True, exist_ok=True)

daiquiri.setup(level="INFO")       # DEBUG floods the log at large n


# =================================================================
# Demography
# =================================================================
def scale_eur_final_size(base_demog, target_final):
    d = copy.deepcopy(base_demog)
    d.add_population_parameters_change(
        time=0, population=POP, initial_size=target_final,
        growth_rate=EUR_GROWTH_RATE,
    )
    d.sort_events()
    return d

species = stdpopsim.get_species("HomSap")
base = species.get_demographic_model(DEMOG_MODEL).model
demog = scale_eur_final_size(base, N_FINAL[N_DIP])

print(f"n={N_DIP:,}  N_final={N_FINAL[N_DIP]:,}  L={L:.0e}  "
      f"sample_frac={N_DIP / N_FINAL[N_DIP]:.4f}", flush=True)


# =================================================================
# Simulate (cached)
# =================================================================
tree_file = OUT / "sim.trees"
if tree_file.is_file():
    import tskit
    print(f"Loading cached {tree_file}", flush=True)
    ts = tskit.load(str(tree_file))
else:
    ts = msprime.sim_ancestry(
        samples={POP: N_DIP},
        demography=demog,
        sequence_length=L,
        recombination_rate=REC,
        model=[msprime.DiscreteTimeWrightFisher(duration=100),
               msprime.StandardCoalescent()],
        random_seed=SEED,
    )
    ts = msprime.sim_mutations(ts, rate=MU, random_seed=SEED)
    ts.dump(str(tree_file))

print(f"{ts.num_sites:,} sites, {ts.num_mutations:,} mutations", flush=True)


# =================================================================
# Singletons straight from the tables
# =================================================================
tabs = ts.tables
mut_node = tabs.mutations.node
mt = tabs.mutations.time
ages_all = np.where(np.isnan(mt), tabs.nodes.time[mut_node], mt)

is_sample = np.zeros(ts.num_nodes, dtype=bool)
is_sample[ts.samples()] = True
sing = is_sample[mut_node]

carrier = tabs.nodes.individual[mut_node[sing]]
sing_ages = ages_all[sing]

if (carrier < 0).any():
    sys.exit("sample nodes have no individual — cannot map carriers")

n_sing = int(sing.sum())
print(f"{n_sing:,} singletons  frac={n_sing / ts.num_sites:.3f}")
print("age_q (50/90/95/99/99.9):",
      np.percentile(sing_ages, [50, 90, 95, 99, 99.9]).round(0), flush=True)

np.save(OUT / "sing_ages.npy", sing_ages)
np.save(OUT / "carrier.npy", carrier)


# =================================================================
# Old-variant supply: the diagnostic that predicts the lever arm
# =================================================================
rows = []
print(f"\n{'threshold':>10}{'n_old':>12}{'mean/person':>13}")
for th in [50, 100, 200, 500, 1000, 2000]:
    n_old = int((sing_ages > th).sum())
    print(f"{th:>10}{n_old:>12,}{n_old / N_DIP:>13.3f}")
    rows.append({"n_dip": N_DIP, "L": L, "threshold": th,
                 "n_old": n_old, "mean_per_person": n_old / N_DIP})
pd.DataFrame(rows).to_csv(OUT / "old_variant_supply.csv", index=False)


# =================================================================
# Equal-count bins over the young range, plus a catch-all
# =================================================================
a_young = sing_ages[sing_ages < YOUNG_MAX]
if a_young.size == 0:
    sys.exit(f"no singletons younger than {YOUNG_MAX} generations")

summaries = []
for nb in BIN_SCHEMES:
    q = np.quantile(a_young, np.linspace(0, 1, nb + 1))
    q[0], q[-1] = 0.0, YOUNG_MAX
    q = np.unique(q)                       # guard against duplicate edges
    bins = np.concatenate([q, [np.inf]])
    n_bins = len(bins) - 1

    idx = np.searchsorted(bins, sing_ages, side="right") - 1
    counts = np.bincount(
        carrier * n_bins + idx, minlength=N_DIP * n_bins
    ).reshape(N_DIP, n_bins)

    labels = [f"{lo:.0f}+" if np.isinf(hi) else f"{lo:.0f}-{hi:.0f}"
              for lo, hi in zip(bins[:-1], bins[1:])]

    assert counts.sum() == n_sing, "binning lost singletons"

    pd.DataFrame(counts, columns=labels).to_csv(
        OUT / f"singleton_counts_nb{nb}.csv", index=False)

    print(f"\nnb={nb}  edges={np.round(bins[:-1], 1)}")
    print(f"{'bin':<16}{'n_sing':>12}{'mean':>10}{'var/mean':>10}"
          f"{'n_eff':>10}{'frac0':>8}")
    for b, lab in enumerate(labels):
        c = counts[:, b].astype(float)
        if c.sum() == 0:
            continue
        m = c.mean()
        print(f"{lab:<16}{int(c.sum()):>12,}{m:>10.2f}"
              f"{c.var() / m:>10.2f}"
              f"{c.sum() ** 2 / (c ** 2).sum():>10.0f}"
              f"{(c == 0).mean():>8.3f}")
        summaries.append({
            "n_dip": N_DIP, "L": L, "nb": nb, "bin": lab,
            "bin_lo": bins[b], "bin_hi": bins[b + 1],
            "n_sing": int(c.sum()), "mean_c": m,
            "var_mean": c.var() / m,
            "n_eff": c.sum() ** 2 / (c ** 2).sum(),
            "frac_zero": (c == 0).mean(),
        })

pd.DataFrame(summaries).to_csv(OUT / "bin_summary.csv", index=False)
print(f"\nWrote outputs to {OUT}", flush=True)
