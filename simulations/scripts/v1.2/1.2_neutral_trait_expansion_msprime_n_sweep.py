#!/usr/bin/env python
"""Production singleton-count runs at UKB-matched ascertainment.
Usage: python run_production.py <n_dip> [L]
"""
import sys, copy
import numpy as np, pandas as pd
import msprime, stdpopsim, tskit
import daiquiri

N_DIP = int(sys.argv[1])
L = float(sys.argv[2]) if len(sys.argv) > 2 else 1e8

# calibrated N_final; 500k from the n^0.699 fit
N_FINAL = {20_000: 563_593, 50_000: 941_197, 100_000: 1_583_484,
           200_000: 2_822_145, 500_000: 5_350_000}[N_DIP]

POP, MU, REC = "EUR", 1.25e-8, 1e-8
EUR_GROWTH_RATE, SEED = 0.0195, 1

species = stdpopsim.get_species("HomSap")
base = species.get_demographic_model("OutOfAfrica_2T12").model
daiquiri.setup(level="DEBUG")

def scale_eur_final_size(d0, target):
    d = copy.deepcopy(d0)
    d.add_population_parameters_change(
        time=0, population=POP, initial_size=target,
        growth_rate=EUR_GROWTH_RATE)
    d.sort_events()
    return d

print(f"n={N_DIP:,}  N_final={N_FINAL:,}  L={L:.0e}  "
      f"sample_frac={N_DIP/N_FINAL:.4f}", flush=True)

ts = msprime.sim_ancestry(
    samples={POP: N_DIP},
    demography=scale_eur_final_size(base, N_FINAL),
    sequence_length=L,
    recombination_rate=REC,
    model=[msprime.DiscreteTimeWrightFisher(duration=100),
           msprime.StandardCoalescent()],
    random_seed=SEED,
)
ts = msprime.sim_mutations(ts, rate=MU, random_seed=SEED)
ts.dump(str(SIM_PATH / f"n{N_DIP}_L{L:.0e}.trees"))
print(f"{ts.num_sites:,} sites, {ts.num_mutations:,} mutations", flush=True)

# --- singletons straight from the tables ---
tabs = ts.tables
mut_node = tabs.mutations.node
mt = tabs.mutations.time
ages_all = np.where(np.isnan(mt), tabs.nodes.time[mut_node], mt)

is_sample = np.zeros(ts.num_nodes, dtype=bool)
is_sample[ts.samples()] = True
sing = is_sample[mut_node]
carrier = tabs.nodes.individual[mut_node[sing]]
sing_ages = ages_all[sing]

print(f"{sing.sum():,} singletons  frac={sing.sum()/ts.num_sites:.3f}")
print("age_q:", np.percentile(sing_ages, [50, 90, 95, 99, 99.9]).round(0))

# --- the diagnostic that predicts the lever arm ---
print(f"\n{'threshold':>10}{'n_old':>10}{'mean/person':>13}")
for th in [50, 100, 200, 500, 1000, 2000]:
    n_old = int((sing_ages > th).sum())
    print(f"{th:>10}{n_old:>10,}{n_old/N_DIP:>13.3f}")

# --- equal-count bins over the young range, plus a catch-all ---
YOUNG_MAX = 205.0
for nb in [1, 2, 3, 4, 6]:
    a = sing_ages[sing_ages < YOUNG_MAX]
    q = np.quantile(a, np.linspace(0, 1, nb + 1))
    q[0], q[-1] = 0.0, YOUNG_MAX
    bins = np.concatenate([q, [np.inf]])

    idx = np.searchsorted(bins, sing_ages, side="right") - 1
    counts = np.bincount(carrier * (len(bins) - 1) + idx,
                         minlength=N_DIP * (len(bins) - 1)
                         ).reshape(N_DIP, len(bins) - 1)
    labels = [f"{lo:.0f}+" if np.isinf(hi) else f"{lo:.0f}-{hi:.0f}"
              for lo, hi in zip(bins[:-1], bins[1:])]

    pd.DataFrame(counts, columns=labels).to_csv(
        SIM_PATH / f"n{N_DIP}_singleton_counts_nb{nb}.csv", index=False)

    print(f"\nnb={nb}")
    print(f"{'bin':<16}{'n_sing':>10}{'mean':>9}{'var/mean':>10}"
          f"{'n_eff':>9}{'frac0':>8}")
    for b, lab in enumerate(labels):
        c = counts[:, b].astype(float)
        if c.sum() == 0:
            continue
        print(f"{lab:<16}{int(c.sum()):>10,}{c.mean():>9.2f}"
              f"{c.var()/c.mean():>10.2f}"
              f"{c.sum()**2/(c**2).sum():>9.0f}{(c==0).mean():>8.3f}")

np.save(SIM_PATH / f"n{N_DIP}_sing_ages.npy", sing_ages)
np.save(SIM_PATH / f"n{N_DIP}_carrier.npy", carrier)
