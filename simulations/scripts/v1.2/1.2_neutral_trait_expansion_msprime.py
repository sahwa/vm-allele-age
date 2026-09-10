
import sys
from pathlib import Path
import copy

import msprime
import stdpopsim
import numpy as np
import pandas as pd
import tskit


SIM_VERSION = "1.2"
SIM_TYPE = "neutral_trait_expansion"

#### stuff #####
REP = 0

RNG_SEED = 42

ss = np.random.SeedSequence([RNG_SEED, REP])
seed_mutations, seed_beta, seed_sample, seed_noise = ss.spawn(4)

rng_beta = np.random.default_rng(seed_beta)
rng_sample = np.random.default_rng(seed_sample)
rng_noise = np.random.default_rng(seed_noise)

RECOMBINATION_RATE = 1e-8   # must match initializeRecombinationRate() in the .slim script

MU = 1.44e-8                # per base per generation mutation rate
PI_TARGET = 0.01            # fraction of the genome that is mutational target

SIM_PATH = Path(
    f"/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v{SIM_VERSION}"
)

SIM_PATH_REP = SIM_PATH / "replicates" / f"rep{REP}"
SIM_PATH_REP.mkdir(parents=True, exist_ok=True)

# N_SAMPLE_TARGET = 50000     # diploids to sample for the GREML/GENIE analysis
H2 = 0.5
SIGMA_BETA = 1.0

BINS = [0, 1e2, 1e3, 1e4, 5e4, 1e5, 2e5, 5e5, np.inf]

# The mutation-overlaid tree sequence is large and only needed if you want to
# re-analyse this replicate later (e.g. tsdate). Off by default for big runs.
SAVE_MUT_TREES = True

TREE_FILE = SIM_PATH / f"{SIM_VERSION}_{SIM_TYPE}.trees"   # shared
RTS_FILE = SIM_PATH / f"{SIM_VERSION}_{SIM_TYPE}.recapitated.trees"
MUT_TREE_FILE = SIM_PATH_REP / f"{SIM_VERSION}_{SIM_TYPE}.recap.mut.trees"
GENOME_PREFIX = SIM_PATH_REP / f"{SIM_VERSION}_{SIM_TYPE}"

#### demographic parameters 
MU = 1.25e-8
L = 1e8
REC = 1e-8
N_FINAL_TARGET = 500_000   # population size, goes into scale_eur_final_size
N_SAMPLE = 20_000          # diploids actually simulated as samples
base = species.get_demographic_model("OutOfAfrica_2T12").model


def scale_eur_final_size(base_demog, target_final):
    demog = copy.deepcopy(base_demog)
    demog.add_population_parameters_change(
        time=0, population="EUR",
        initial_size=target_final, growth_rate=0.0195,
    )
    return demog

demog = scale_eur_final_size(base, N_FINAL_TARGET)

species = stdpopsim.get_species("HomSap")
name = "OutOfAfrica_2T12"
pop = "EUR"

demog = species.get_demographic_model(name).model

ts = msprime.sim_ancestry(
    samples = {pop: N_SAMPLE},
    demography = demog,
    sequence_length = L,
    recombination_rate = REC,
    model = [msprime.DiscreteTimeWrightFisher(duration=100),
                msprime.StandardCoalescent()],
    random_seed = 1
)


# save here since this step takes a long time


#### sim the mutations on the treee ###

ts = msprime.sim_mutations(ts, rate=MU, random_seed=1)

## Do some checks to infer 

afs = ts.allele_frequency_spectrum(polarised = True, span_normalise = False)
n_sing, n_sites = int(afs[1]), ts.num_sites

mut_time = ts.tables.mutations.time
site_of_mut = ts.tables.mutations.site
counts = np.bincount(site_of_mut, minlength=ts.num_sites)

ok = counts == 1
ac = np.array([v.genotypes.sum() for v in ts.variants()])
sing_sites = ok & (ac == 1)
ages = mut_time[np.isin(site_of_mut, np.flatnonzero(sing_sites))]

print(f"{name:18s} n={n_dip:6d} sites={n_sites:7d} "
        f"sing={n_sing:7d} frac={n_sing/n_sites:.3f} "
        f"age_q={np.percentile(ages, [50, 90, 95, 99, 99.9, 99.99]).round(0)}")


# =================================================================
# 6. Effect sizes (neutral trait: independent of age and frequency)
# =================================================================
beta = rng_beta.normal(0, SIGMA_BETA, size=M_raw)

# =================================================================
# 7. Ground-truth V_M and V_A checks
# =================================================================
u_target = MU * PI_TARGET * ts.sequence_length
V_M_true = 2 * u_target * SIGMA_BETA**2
V_A_analytic = 2 * SIM_NE * V_M_true
V_A_empirical = np.sum(2 * freqs * (1 - freqs) * beta**2)

print(f"\nV_M (true, analytic):            {V_M_true:.4g}")
print(f"V_A (analytic 2*Ne*V_M):         {V_A_analytic:.4g}")
print(f"V_A (empirical, sum 2pq*beta^2): {V_A_empirical:.4g}")
print(f"Ratio empirical/analytic: {V_A_empirical / V_A_analytic:.3f}  "
      f"(should be close to 1 if at equilibrium)", flush=True)


# =================================================================
# 5. Single pass over variants.
#    Collects, in one traversal: age, position, ref/alt alleles, population
#    allele frequency, and the dosage matrix for the sampled individuals.
#    This is the expensive step; nothing else re-traverses the tree sequence.
# =================================================================
n_dip_all = ts.num_samples // 2
n_sample = min(N_SAMPLE_TARGET, n_dip_all)
sample_idx = np.sort(rng_sample.choice(n_dip_all, n_sample, replace=False))

# haplotype indices for the sampled diploids, interleaved (2i, 2i+1)
hap_idx = np.stack([2 * sample_idx, 2 * sample_idx + 1], axis=1).ravel()

M_raw = mts.num_sites
dosages = np.empty((M_raw, n_sample), dtype=np.uint8)
ages = np.empty(M_raw)
freqs = np.empty(M_raw)
positions = np.empty(M_raw)
ref = np.empty(M_raw, dtype="<U1")
alt = np.empty(M_raw, dtype="<U1")

print(f"\nSingle pass over {M_raw} variants "
      f"({n_sample} of {n_dip_all} diploids sampled)...", flush=True)

for i, var in enumerate(mts.variants()):
    gt = var.genotypes
    mut = mts.mutation(var.site.mutations[0].id)

    ages[i] = mut.time if not np.isnan(mut.time) else mts.node(mut.node).time
    freqs[i] = gt.mean()                        # population allele frequency
    positions[i] = var.site.position
    ref[i] = var.alleles[0]
    alt[i] = var.alleles[1]
    dosages[i] = gt[hap_idx].reshape(-1, 2).sum(axis=1)

print(f"Age range: {ages.min():.1f} - {ages.max():.1f} generations")
print(f"Freq range: {freqs.min():.5f} - {freqs.max():.5f}", flush=True)

# =================================================================
# 9. True per-bin variance, against the analytic neutral kernel.
#    Computed over the full population (pre-MAF-filter): this is the
#    quantity GENIE's per-bin estimates are being asked to recover.
# =================================================================
bin_labels = []
for i in range(len(BINS) - 1):
    lo = int(BINS[i])
    hi = "inf" if np.isinf(BINS[i + 1]) else int(BINS[i + 1])
    bin_labels.append(f"{lo}-{hi}")

bins_assigned = pd.cut(ages, bins=BINS, labels=bin_labels, right=False)

print(f"\n{'Bin (gens)':<18}{'n':>7}{'observed':>12}{'predicted':>12}{'ratio':>9}")
bin_rows = []
for lo, hi in zip(BINS[:-1], BINS[1:]):
    m = (ages >= lo) & (ages < hi)
    if m.sum() == 0:
        continue
    contribs = 2 * freqs[m] * (1 - freqs[m]) * beta[m]**2
    V_bin = contribs.sum()

    # Integral of K(t) = exp(-t/2Ne) across the bin
    hi_term = 0.0 if np.isinf(hi) else np.exp(-hi / (2 * SIM_NE))
    K_int = 2 * SIM_NE * (np.exp(-lo / (2 * SIM_NE)) - hi_term)
    V_pred = V_M_true * K_int

    # Kish effective n: bin variance is dominated by a few large contributors
    n_eff = contribs.sum()**2 / (contribs**2).sum()

    label = f"{lo:.0f}+" if np.isinf(hi) else f"{lo:.0f}-{hi:.0f}"
    print(f"{label:<18}{m.sum():>7}{V_bin:>12.4g}{V_pred:>12.4g}"
          f"{V_bin / V_pred:>9.3f}")

    bin_rows.append({
        "bin_lo": lo, "bin_hi": hi, "n_variants": int(m.sum()),
        "n_eff": n_eff, "V_observed": V_bin, "V_predicted": V_pred,
    })

pd.DataFrame(bin_rows).to_csv(
    SIM_PATH_REP / f"{SIM_VERSION}_bin_truth.csv", index=False)


# =================================================================
# 10. MAF filter in the sample, applied to every per-variant array at once.
#     GENIE requires MAF > 0, and a variant monomorphic in the sampled 50,000
#     carries no information. Filtering all arrays together is what keeps the
#     annotation aligned to the .bim without any position matching.
# =================================================================
sample_freq = dosages.mean(axis=1) / 2
keep = (sample_freq > 0) & (sample_freq < 1)

print(f"\nAfter MAF filtering: {keep.sum()} of {M_raw} variants retained "
      f"({M_raw - keep.sum()} dropped as monomorphic in the sample)", flush=True)

dosages_f = dosages[keep]
positions_f = positions[keep]
ages_f = ages[keep]
ref_f, alt_f = ref[keep], alt[keep]

indv_names = [f"ind{i}" for i in range(n_sample)]

print("Writing plink fileset...", flush=True)
write_plink_bed(dosages_f, GENOME_PREFIX, positions_f, ref_f, alt_f, indv_names)


# =================================================================
# 11. GENIE annotation matrix.
#     Row i corresponds to row i of the .bim by construction, since both come
#     from the same `keep` mask applied to the same arrays.
# =================================================================
bins_assigned_filtered = pd.cut(ages_f, bins=BINS, labels=bin_labels, right=False)

annotations = pd.get_dummies(
    bins_assigned_filtered, prefix="bin", prefix_sep="_"
).astype(int)

assert len(annotations) == keep.sum(), "annotation rows != retained variants"
assert (annotations.sum(axis=1) == 1).all(), "each variant must fall in exactly one bin"

annotations.to_csv(
    SIM_PATH_REP / f"{SIM_VERSION}_annotations_age_bins.txt",
    sep=" ", index=False, header=False,
)

pd.DataFrame({
    "column_name": annotations.columns,
    "age_bin": bin_labels,
}).to_csv(
    SIM_PATH_REP / f"{SIM_VERSION}_annotations_legend.txt",
    sep=" ", index=False,
)


# =================================================================
# 12. Variant info and phenotype files
# =================================================================
pd.DataFrame({
    "site_id": np.arange(M_raw),
    "position": positions,
    "age": ages,
    "bin": bins_assigned,
    "freq": freqs,
    "beta": beta,
    "kept": keep,
}).to_csv(SIM_PATH_REP / f"{SIM_VERSION}_variant_info.csv", index=False)

pd.DataFrame({"FID": 0, "IID": indv_names, "y": y}).to_csv(
    SIM_PATH_REP / f"{SIM_VERSION}_phenotypes.GENIE.txt",
    sep="\t", index=False, header=["FID", "IID", "PHENO"]
)

pd.DataFrame({
    "sample_idx": sample_idx,
    "iid": indv_names,
    "y": y,
    "g": g,
}).to_csv(SIM_PATH_REP / f"{SIM_VERSION}_phenotypes.csv", index=False)

print(f"\nWrote outputs to {SIM_PATH_REP}", flush=True)

