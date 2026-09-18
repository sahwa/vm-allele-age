import copy
from pathlib import Path

import msprime
import numpy as np
import pandas as pd
import stdpopsim
import tskit
import daiquiri
# =================================================================
# Configuration
# =================================================================
SIM_VERSION = "1.2.1"
SIM_TYPE = "neutral_trait_expansion"
REP = 0
RNG_SEED = 42

# Demography / simulation
DEMOG_MODEL = "OutOfAfrica_2T12"
POP = "EUR"
N_FINAL_TARGET = 500_000      # present-day EUR size (calibrated to UKB singleton frac)
EUR_GROWTH_RATE = 0.0195      # from the base model; kept when rescaling the endpoint
N_SAMPLE = 20_000             # diploids drawn as samples
L = 1e9                       # sequence length (bp)
MU = 1.25e-8                  # per-bp per-generation mutation rate
REC = 1e-8                    # per-bp per-generation recombination rate

# Trait
H2 = 0.5                      # target narrow-sense heritability
SIGMA_BETA = 1.0              # SD of causal effect sizes
PI_CAUSAL = 0.01              # fraction of segregating sites that are causal

# Age bins (log-spaced; edges chosen to span the observed age range)
BINS = np.concatenate([[0.0], np.geomspace(10, 1e6, 9), [np.inf]])

# Streaming
CHUNK_BP = 8e6                # genomic window per chunk; tune to available RAM

# Paths
SIM_PATH = Path(
    "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data"
) / f"v{SIM_VERSION}"
SIM_PATH_REP = SIM_PATH / "replicates" / f"rep{REP}"
SIM_PATH_REP.mkdir(parents=True, exist_ok=True)

TREE_FILE = SIM_PATH / f"{SIM_VERSION}_{SIM_TYPE}.recapitated.trees"
GENOME_PREFIX = SIM_PATH_REP / f"{SIM_VERSION}_{SIM_TYPE}"

# RNG streams — independent and reproducible per replicate
_ss = np.random.SeedSequence([RNG_SEED, REP])
seed_beta, seed_causal, seed_noise = _ss.spawn(3)
rng_beta = np.random.default_rng(seed_beta)
rng_causal = np.random.default_rng(seed_causal)
rng_noise = np.random.default_rng(seed_noise)

daiquiri.setup(level="DEBUG")


species = stdpopsim.get_species("HomSap")
base = species.get_demographic_model("OutOfAfrica_2T12").model

SAMPLE_SIZES = [20_000, 50_000, 100_000, 200_000, 500_000]
TARGET_FRAC = 0.46          # UKB singleton fraction
L_CAL = 1e6                 # small: fraction is a ratio, converges fast
N_BISECT = 6

CAL_FILE = SIM_PATH / f"{SIM_VERSION}_Nfinal_calibration.csv"

def scale_eur_final_size(base_demog, target_final):
    demog = copy.deepcopy(base_demog)
    demog.add_population_parameters_change(
        time=0, population="EUR",
        initial_size=target_final, growth_rate=0.0195,
    )
    demog.sort_events()
    return demog


def sim_and_summarise(n_dip, target_final, L=L_CAL, seed=1):
    """One calibration run. Returns a dict of summary stats."""
    demog = scale_eur_final_size(base, target_final)
    ts = msprime.sim_ancestry(
        samples={POP: n_dip},
        demography=demog,
        sequence_length=L,
        recombination_rate=REC,
        model=[msprime.DiscreteTimeWrightFisher(duration=100),
               msprime.StandardCoalescent()],
        random_seed=seed,
    )
    ts = msprime.sim_mutations(ts, rate=MU, random_seed=seed)

    afs = ts.allele_frequency_spectrum(polarised=True, span_normalise=False)
    n_sing, n_sites = int(afs[1]), ts.num_sites

    # singleton ages straight from the tables (no variant traversal)
    tabs = ts.tables
    mut_node = tabs.mutations.node
    mt = tabs.mutations.time
    ages_all = np.where(np.isnan(mt), tabs.nodes.time[mut_node], mt)
    is_sample = np.zeros(ts.num_nodes, dtype=bool)
    is_sample[ts.samples()] = True
    sing_ages = ages_all[is_sample[mut_node]]

    q = np.percentile(sing_ages, [50, 90, 95, 99, 99.9]) if sing_ages.size else [np.nan]*5

    return {
        "n_dip": n_dip,
        "N_final": target_final,
        "sample_frac": n_dip / target_final,
        "n_sites": n_sites,
        "n_sing": n_sing,
        "sing_frac": n_sing / n_sites if n_sites else np.nan,
        "age_p50": q[0], "age_p90": q[1], "age_p95": q[2],
        "age_p99": q[3], "age_p999": q[4],
        "n_mutations": ts.num_mutations,
        "L": L,
    }


rows = []
chosen = {}

for n_dip in SAMPLE_SIZES:
    lo, hi = n_dip * 5, n_dip * 200
    print(f"\n=== n = {n_dip:,} ===", flush=True)

    for it in range(N_BISECT):
        mid = np.sqrt(lo * hi)                      # geometric bisection
        r = sim_and_summarise(n_dip, mid)
        r["iter"] = it
        rows.append(r)

        print(f"  N={r['N_final']:>12,.0f}  frac_samp={r['sample_frac']:.4f}  "
              f"sites={r['n_sites']:>8,}  sing={r['n_sing']:>8,}  "
              f"sing_frac={r['sing_frac']:.3f}  "
              f"age_q={[int(r[k]) for k in ('age_p50','age_p90','age_p99')]}",
              flush=True)

        if r["sing_frac"] < TARGET_FRAC:
            lo = mid
        else:
            hi = mid

        # keep the run closest to target for this n
        best = chosen.get(n_dip)
        if best is None or abs(r["sing_frac"] - TARGET_FRAC) < abs(best["sing_frac"] - TARGET_FRAC):
            chosen[n_dip] = r

    c = chosen[n_dip]
    print(f"  -> chosen N_final = {c['N_final']:,.0f} "
          f"(sing_frac {c['sing_frac']:.3f}, multiplier {c['N_final']/n_dip:.1f}x)",
          flush=True)


cal = pd.DataFrame(rows)
cal.to_csv(CAL_FILE, index=False)

best = pd.DataFrame(chosen.values()).assign(
    multiplier=lambda d: d.N_final / d.n_dip)
best.to_csv(SIM_PATH / f"{SIM_VERSION}_Nfinal_chosen.csv", index=False)

print("\n=== chosen N_final per sample size ===")
print(best[["n_dip", "N_final", "multiplier", "sample_frac",
            "sing_frac", "age_p50", "age_p99"]].to_string(
    index=False, float_format=lambda v: f"{v:,.3f}"))
print(f"\nFull calibration trace -> {CAL_FILE}")
