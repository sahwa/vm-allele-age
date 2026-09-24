#!/usr/bin/env python
"""
v1.2 — neutral trait under Tennessen-style European expansion.

Simulates a neutral polygenic trait on an msprime tree sequence, bins variants
by true allele age, and writes a PLINK fileset + GENIE annotation matrix.

Memory note: the genotype matrix is never held in full. Variants are streamed in
genomic chunks; each chunk is packed straight into the .bed file and its
contribution to each individual's genetic value is accumulated on the fly.
Only per-variant metadata (O(M) floats) and per-individual vectors (O(N)) persist.
"""

import copy
from pathlib import Path

import msprime
import numpy as np
import pandas as pd
import stdpopsim
import tskit

# =================================================================
# Configuration
# =================================================================
SIM_VERSION = "1.2"
SIM_TYPE = "neutral_trait_expansion"
REP = 0
RNG_SEED = 42

# Demography / simulation
DEMOG_MODEL = "OutOfAfrica_2T12"
POP = "EUR"
N_FINAL_TARGET = 500_000      # present-day EUR size (calibrated to UKB singleton frac)
EUR_GROWTH_RATE = 0.0195      # from the base model; kept when rescaling the endpoint
N_SAMPLE = 20_000             # diploids drawn as samples
L = 1e8                       # sequence length (bp)
MU = 1.25e-8                  # per-bp per-generation mutation rate
REC = 1e-8                    # per-bp per-generation recombination rate

# Trait
H2_SING = 0.1                      # target narrow-sense heritability
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


# =================================================================
# 1. Load caches and strip out any multi-allelics
# =================================================================
if TREE_FILE.is_file():
    print(f"Loading cached tree sequence: {TREE_FILE}", flush=True)
    ts = tskit.load(str(TREE_FILE))
else:
    print(f"Tree doesn't exist. Exiting.....")
    exit()

multiallelic_site_ids = np.array(
    [site.id for site in ts.sites() if len(site.mutations) > 1]
)

print(f"Removing {len(multiallelic_site_ids)} multiallelic sites "
      f"out of {ts.num_sites}", flush=True)

if len(multiallelic_site_ids) > 0:
    tables = ts.dump_tables()
    tables.delete_sites(multiallelic_site_ids)
    tables.sort()
    ts = tables.tree_sequence()
print(f"mts now has {ts.num_sites} biallelic sites", flush=True)

N_IND = ts.num_samples // 2
M_RAW = ts.num_sites
assert ts.num_mutations == M_RAW, (
    "Multi-mutation sites present — table-column alignment below assumes "
    "exactly one mutation per site."
)
print(f"{M_RAW} sites, {N_IND} diploids", flush=True)

tabs = ts.tables
mut_node = tabs.mutations.node                  # node each mutation sits above

is_sample = np.zeros(ts.num_nodes, dtype=bool)
is_sample[ts.samples()] = True
sing = is_sample[mut_node]                      # singleton iff on a sample node
carrier = tabs.nodes.individual[mut_node[sing]] # node -> individual ID

# =================================================================
# 2. Effect sizes (neutral trait: beta independent of age and frequency)
# =================================================================
afs = ts.allele_frequency_spectrum(polarised = True, span_normalise = False)
n_sing, n_sites = int(afs[1]), ts.num_sites

mut_time = ts.tables.mutations.time
site_of_mut = ts.tables.mutations.site
counts = np.bincount(site_of_mut, minlength=ts.num_sites)

ok = counts == 1
ac = np.array([v.genotypes.sum() for v in ts.variants()])
sing_sites = ok & (ac == 1)
M_SING = sing_sites.sum()
# ages = mut_time[np.isin(site_of_mut, np.flatnonzero(sing_sites))]

causal = rng_causal.random(M_SING) < PI_CAUSAL
beta = np.zeros(M_RAW)
SIGMA2_B = H2_SING * N_IND / M_RAW    # from n=200k, L=1e9, say
beta_sing = rng_beta.normal(0, np.sqrt(SIGMA2_B), size=len(carrier))
g_sing = np.bincount(carrier, weights=beta_sing, minlength=N_IND)
g_sing -= g_sing.mean()

V_A_ref = g_sing.var()
SIGMA2_E = V_A_ref - (1-H2_SING / H2_SING)

y = g_sing + rng_noise.normal(0, np.sqrt(SIGMA2_E), size=N_IND)
np.bincount(carrier, weights=BETA_SING)

# =================================================================
# 3. Streaming pass: pack .bed, accumulate genetic values, collect metadata
# =================================================================
g_sing = np.bincount(carrier, weights=BETA_SING, minlength=N_IND)

iids = [f"ind{i}" for i in range(N_IND)]

# =================================================================
# 4. Phenotype: y = g + e, scaled to the target heritability
# =================================================================
g -= g.mean()
var_g = g.var()
sigma_e = np.sqrt(var_g * (1 - H2) / H2)
y = g + rng_noise.normal(0, sigma_e, size=N_IND)

print(f"\nvar(g)={var_g:.4g}  sigma_e={sigma_e:.4g}  "
      f"realised h2={var_g / y.var():.3f}", flush=True)


# =================================================================
# 5. Ground truth: V_M and V_A
# =================================================================
u_causal = MU * PI_CAUSAL * L                 # causal mutations per generation
V_M_true = 2 * u_causal * SIGMA_BETA**2
V_A_empirical = np.sum(2 * freqs * (1 - freqs) * beta**2)

print(f"\nV_M (true, analytic):            {V_M_true:.4g}")
print(f"V_A (empirical, sum 2pq*beta^2): {V_A_empirical:.4g}")
print(f"Implied T = V_A / V_M:           {V_A_empirical / V_M_true:.4g} generations",
      flush=True)


N_BINS = len(BINS) - 1
sing_counts = np.zeros((N_IND, N_BINS), dtype=np.int32)   # c_i^(t)
all_counts  = np.zeros((N_IND, N_BINS), dtype=np.int32)   # optional: all variants

# =================================================================
# 6. Singleton Counts
# =================================================================

ac = np.round(freqs * ts.num_samples).astype(int)
sing = (ac == 1) & keep

bin_idx = np.searchsorted(BINS, ages, side="right") - 1

print(f"\n{'Bin (gens)':<20}{'n_variants':>12}{'n_singleton':>13}{'sing_frac':>11}")
for b, (lo, hi) in enumerate(zip(BINS[:-1], BINS[1:])):
    m = (bin_idx == b) & keep
    n_var = m.sum()
    n_sing = (m & sing).sum()
    if n_var == 0:
        continue
    label = f"{lo:.0f}+" if np.isinf(hi) else f"{lo:.0f}-{hi:.0f}"
    print(f"{label:<20}{n_var:>12}{n_sing:>13}{n_sing/n_var:>11.3f}")

print(f"\nTotal: {keep.sum()} variants, {sing.sum()} singletons "
      f"({sing.sum()/keep.sum():.3f})")

# =================================================================
# 7. True per-bin genetic variance
# =================================================================
rows = []
print(f"\n{'Bin (gens)':<20}{'n':>9}{'n_eff':>10}{'V_observed':>13}{'share':>9}")
for lo, hi in zip(BINS[:-1], BINS[1:]):
    m = (ages >= lo) & (ages < hi)
    if not m.any():
        continue
    contribs = 2 * freqs[m] * (1 - freqs[m]) * beta[m]**2
    V_bin = contribs.sum()
    ss = (contribs**2).sum()
    n_eff = V_bin**2 / ss if ss > 0 else 0.0   # Kish effective n

    label = f"{lo:.0f}+" if np.isinf(hi) else f"{lo:.0f}-{hi:.0f}"
    print(f"{label:<20}{m.sum():>9}{n_eff:>10.1f}{V_bin:>13.4g}"
          f"{V_bin / V_A_empirical:>9.3f}")

    rows.append({
        "bin_lo": lo, "bin_hi": hi, "n_variants": int(m.sum()),
        "n_eff": n_eff, "V_observed": V_bin,
        "true_share": V_bin / V_A_empirical,
    })

pd.DataFrame(rows).to_csv(
    SIM_PATH_REP / f"{SIM_VERSION}_bin_truth.csv", index=False)


# =================================================================
# 8. GENIE annotation matrix (row i aligns with row i of the .bim)
# =================================================================
bin_labels = [
    f"{lo:.0f}+" if np.isinf(hi) else f"{lo:.0f}-{hi:.0f}"
    for lo, hi in zip(BINS[:-1], BINS[1:])
]
bins_all = pd.cut(ages, bins=BINS, labels=bin_labels, right=False)

annotations = pd.get_dummies(
    bins_all[keep], prefix="bin", prefix_sep="_"
).reindex(columns=[f"bin_{b}" for b in bin_labels], fill_value=0).astype(int)

assert len(annotations) == keep.sum(), "annotation rows != retained variants"
assert (annotations.sum(axis=1) == 1).all(), "each variant must fall in exactly one bin"

annotations.to_csv(
    SIM_PATH_REP / f"{SIM_VERSION}_annotations_age_bins.txt",
    sep=" ", index=False, header=False,
)
pd.DataFrame({"column_name": annotations.columns, "age_bin": bin_labels}).to_csv(
    SIM_PATH_REP / f"{SIM_VERSION}_annotations_legend.txt", sep=" ", index=False,
)

# =================================================================
# 9. Variant info and phenotypes
# =================================================================
pd.DataFrame({
    "site_id": np.arange(M_RAW),
    "position": positions,
    "age": ages,
    "bin": bins_all,
    "freq": freqs,
    "beta": beta,
    "causal": causal,
    "kept": keep,
}).to_csv(SIM_PATH_REP / f"{SIM_VERSION}_variant_info.csv", index=False)

pd.DataFrame({"FID": 0, "IID": iids, "PHENO": y}).to_csv(
    SIM_PATH_REP / f"{SIM_VERSION}_phenotypes.GENIE.txt",
    sep="\t", index=False,
)
pd.DataFrame({"iid": iids, "y": y, "g": g}).to_csv(
    SIM_PATH_REP / f"{SIM_VERSION}_phenotypes.csv", index=False)

print(f"\nWrote outputs to {SIM_PATH_REP}", flush=True)
