"""
Full pipeline: SLiM tree sequence -> recapitate -> overlay mutations ->
extract true ages, frequencies, effect sizes -> phenotypes -> per-bin plink
files -> GENIE annotation matrix.

For a NEUTRAL trait simulation (single population, mutation-drift equilibrium).
"""
# %

import os
import subprocess
from pathlib import Path
import sys
import itertools

import msprime
import numpy as np
import pandas as pd
import pyslim
import tskit
import scipy

# =================================================================
# 0. Parameters
# =================================================================

# %

# REP = int(sys.argv[1])
REP = 1

SIM_NE = 20000
RNG_SEED = 42

ss = np.random.SeedSequence([RNG_SEED, REP])
seed_mutations, seed_beta, seed_sample, seed_noise = ss.spawn(4)

rng_beta = np.random.default_rng(seed_beta)
rng_sample = np.random.default_rng(seed_sample)
rng_noise = np.random.default_rng(seed_noise)

RECOMBINATION_RATE = 1e-8   # must match initializeRecombinationRate() in the .slim script

MU = 1.44e-8                # per base per generation mutation rate
PI_TARGET = 0.01             # fraction of the genome that is mutational target

V_S = 5

SIM_PATH = Path(
    "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v2.0"
)
SIM_VERSION = "2.0"

SIM_PATH_REP = SIM_PATH / f"rep{REP}"
SIM_PATH_REP.mkdir(exist_ok=True)

N_SAMPLE_TARGET = 50000      # diploids to sample for the GREML/GENIE analysis
H2 = 0.5
SIGMA_BETA = 1.0

BINS = [0, 1e2, 1e3, 1e4, 5e4, 1e5, 2e5, 5e5, np.inf]

TREE_FILE = SIM_PATH / f"{SIM_VERSION}_stabilising_selection_VS_{V_S}.trees"  # shared, not per-rep
# MUT_TREE_FILE = SIM_PATH_REP / f"{SIM_VERSION}_neutral_out.recap.mut.trees"
# VCF_FILE = SIM_PATH_REP / f"{SIM_VERSION}_neutral_out.vcf"
# GENOME_PREFIX = SIM_PATH_REP / f"{SIM_VERSION}_neutral_out"

# =================================================================
# 1. Load the SLiM tree sequence
# =================================================================
# Mutations are overlaid *after* the forward simulation rather than during it,
# because the trait is neutral: nothing in the forward sim depends on where
# the mutations landed, so msprime can paint them on afterwards far more
# cheaply than tracking them in SLiM itself.
ts = tskit.load(TREE_FILE)

n_multiroot = sum(1 for t in ts.trees() if t.num_roots > 1)
print(f"Loaded: {ts.num_samples} samples, {ts.num_trees} trees, "
      f"{n_multiroot} not yet coalesced", flush=True)

total_span = 0
weighted_tmrca = 0
for tree in ts.trees():
    root_time = tree.time(tree.root) if tree.num_roots == 1 else None
    if root_time is not None:
        weighted_tmrca += root_time * tree.span
        total_span += tree.span

mean_tmrca = weighted_tmrca / total_span
print(f"Mean TMRCA: {mean_tmrca:.0f} (expect ~4*Ne = {4*SIM_NE} under neutrality)")


# now is a good place to calculate the empirical Ne
# we can do this based on the distribution of coalescence times
# as mean(TMRCA) = 2*ne*mu
rng_ne = np.random.default_rng(42)

sample_nodes = ts.samples()   # the actual valid sample node IDs, not a raw range

rng_ne = np.random.default_rng(42)
idx_pairs = rng_ne.choice(len(sample_nodes), size=(10000, 2), replace=True)
sample_pairs = [(int(sample_nodes[a]), int(sample_nodes[b]))
                for a, b in idx_pairs if a != b]

print(f"Using {len(sample_pairs)} pairs from {len(sample_nodes)} true sample nodes")
print(f"Sample node ID range: {sample_nodes.min()} - {sample_nodes.max()}")

tmrcas = []
n_failed = 0

for tree in itertools.islice(ts.trees(), 0, None, 1000):
    for a, b in sample_pairs:
        try:
            tmrcas.append(tree.tmrca(a, b))
        except ValueError:
            n_failed += 1

tmrcas = np.array(tmrcas)
print(f"Collected {len(tmrcas)}, failed {n_failed}")

if len(tmrcas) > 0:
    print(f"Mean pairwise TMRCA: {tmrcas.mean():.0f}")
    print(f"Ne estimate: {tmrcas.mean()/2:.0f}")

# =================================================================
# 4. Strip multiallelic sites directly from mts.
#    Doing this once, here, means every downstream .variants() loop,
#    the VCF write, and the plink/GENIE outputs are automatically in
#    lockstep - no separate masking or post-hoc row-matching needed.
# # =================================================================
multiallelic_site_ids = np.array(
    [site.id for site in ts.sites() if len(site.mutations) > 1]
)
print(f"Removing {len(multiallelic_site_ids)} multiallelic sites "
      f"out of {ts.num_sites}", flush=True)

tables = ts.dump_tables()
tables.delete_sites(multiallelic_site_ids)
tables.sort()
mts = tables.tree_sequence()
print(f"mts now has {mts.num_sites} biallelic sites", flush=True)

# mts.dump(MUT_TREE_FILE)

# =================================================================
# 5. Extract true ages, frequencies, positions
#    mts is guaranteed biallelic at this point, so no per-variant
#    filtering is needed here.
# =================================================================

ages, freqs, positions, beta = [], [], [], []
for var in ts.variants():
    mut = ts.mutation(var.site.mutations[0].id)
    md = mut.metadata["mutation_list"][0]
    beta.append(md["selection_coeff"])   # SLiM stores the drawn effect here
    age = mut.time if not np.isnan(mut.time) else mts.node(mut.node).time
    ages.append(age)
    freqs.append(var.genotypes.mean())
    positions.append(var.site.position)


ages = np.array(ages)
freqs = np.array(freqs)
positions = np.array(positions)
beta = np.array(beta)

M = len(ages)
print(f"Extracted {M} variants", flush=True)
print(f"Age range: {ages.min():.1f} - {ages.max():.1f} generations")
print(f"Freq range: {freqs.min():.5f} - {freqs.max():.5f}")
print(f"Beta range: {beta.min():.5f} - {beta.max():.5f}")

# =================================================================
# 7. Ground-truth V_M and V_A checks
# =================================================================
# u_target: expected number of target mutations per gamete per generation.

NE_empricial = tmrcas.mean()/2
SIGMA_BETA = 0.1 # we can empirically calculate the SIGMA - but we should use the actual SLiM value
                 # as but the variance will be lower in the sims as the 
                 # deleterious alleles will have been removed by selection (that's the point)


V_A_empirical = np.sum(2 * freqs * (1-freqs) * beta**2)

print(f"SD of extracted betas: {beta.std():.4f} (SLiM param: {0.1})")
u_target = MU * PI_TARGET * ts.sequence_length # mutational target stays the same as before 

V_M_analytic_true = 2 * u_target * SIGMA_BETA**2

PERSISTENCE_TIME = V_A_empirical / V_M_analytic_true


S_BAR = SIGMA_BETA**2 / (2*V_S)

REGIME_PARAMETER = S_BAR * NE_empricial

V_A_analytic_neutral = 2 * NE_empricial * V_M_analytic_true
V_A_analytic_selection = V_M_analytic_true / S_BAR

print(f"N_e × s̄ = {REGIME_PARAMETER:.2f}")
if REGIME_PARAMETER > 1:
    print("  → Selection dominates; expect V_A ≈ V_M / s̄")
elif REGIME_PARAMETER < 1:
    print("  → Drift dominates; expect V_A ≈ 2·N_e·V_M")
else:
    print("  → Intermediate regime; expect neither prediction to be exact")


ratio_to_msb = V_A_empirical / V_A_analytic_selection
print(f"V_A empirical / V_A_MSB: {ratio_to_msb:.3f}  (should be ≈ 1.0)")

ratio_to_neutral = V_A_empirical / V_A_analytic_neutral
print(f"V_A empirical / V_A neutral: {ratio_to_neutral:.3f}  (should be << 1.0 if selection dominates)")

print(f"Persistence time observed: {PERSISTENCE_TIME:.0f} generations")
print(f"Persistence time 1/s̄ (MSB):  {1/S_BAR:.0f} generations")
print(f"Persistence time 2Ne (neutral): {2*NE_empricial:.0f} generations")


# # =================================================================
# # 8. Build genetic values, then sample individuals
# #
# #    NOTE: build g over ALL individuals first, then subset. Passing
# #    samples=hap_idx to .variants() drops variants that are monomorphic
# #    in the subset, which desynchronises the variant enumeration from
# #    beta[i] and silently corrupts g.
# # =================================================================
n_dip_all = ts.num_samples // 2
print(f"\nBuilding genetic values for all {n_dip_all} diploids...", flush=True)
print(ts.num_individuals, n_dip_all)  # do these match?

g_all = np.zeros(n_dip_all)
for i, var in enumerate(ts.variants()):
    geno = var.genotypes.reshape(-1, 2).sum(axis=1)   # diploid dosage 0/1/2
    g_all += beta[i] * geno

print(f"V_A across all individuals: {np.var(g_all):.4g} "
      f"(expect ~{V_A_empirical:.4g}; small gap is LD)", flush=True)

n_sample = min(N_SAMPLE_TARGET, n_dip_all)
sample_idx = np.sort(rng_sample.choice(n_dip_all, n_sample, replace=False))
g = g_all[sample_idx]

V_E = np.var(g) * (1 - H2) / H2
y = g + rng_noise.normal(0, np.sqrt(V_E), n_sample)

print(f"Sampled {n_sample} of {n_dip_all} diploids")
print(f"V_A in sample: {np.var(g):.4g},  V_E: {V_E:.4g}")
print(f"Realised h2: {np.var(g) / np.var(y):.3f} (target {H2})", flush=True)

# =================================================================
# 9. True per-bin variance, compared against the analytic neutral kernel
#    This is the ground truth the REML/GENIE estimates must recover.
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

    # Integral of K(t) = exp(-t/2Ne) across the bin.
    hi_term = 0.0 if np.isinf(hi) else np.exp(-hi / (2 * NE_empricial))
    K_int = 2 * NE_empricial * (
        np.exp(-lo * (1/(2*NE_empricial) + S_BAR)) - 
        (0.0 if np.isinf(hi) else np.exp(-hi * (1/(2*NE_empricial) + S_BAR)))
    )
    V_pred = V_M_analytic_true * K_int

    # Kish effective n: bin variance is dominated by a few large contributors,
    # so the naive sqrt(n) error bar is far too tight.
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
# 10. Write VCF for the sampled individuals
#     individual_names avoids plink's "Sample ID ends with _0" error.
#     mts is already biallelic-only, so no site_mask is needed here.
# =================================================================
indv_names = [f"ind{i}" for i in range(n_sample)]

print(f"\nWriting VCF for {n_sample} sampled individuals...", flush=True)
with open(VCF_FILE, "wt") as f:
    mts.write_vcf(f, individuals=sample_idx, individual_names=indv_names)

# =================================================================
# 11. Genome-wide plink fileset, MAF-filtered.
#     Filtering removes any variant monomorphic in the *sampled*
#     50,000 individuals (GENIE requires MAF > 0). This changes
#     which variants survive relative to `ages`/`bins_assigned`,
#     which were built over the full population - so the annotation
#     must be rebuilt from this .bim, matched by position, not
#     reused from the pre-filter arrays.
# =================================================================
cmd = (f"plink2 --vcf {VCF_FILE} --maf 0.0000001 "
       f"--make-bed --out {GENOME_PREFIX}")
subprocess.run(cmd, shell=True, check=True)

bim = pd.read_csv(f"{GENOME_PREFIX}.bim", sep="\t", header=None,
                   names=["chr", "snpid", "cm", "pos", "a1", "a2"])
print(f"After MAF filtering: {len(bim)} of {M} variants retained "
      f"({M - len(bim)} dropped as monomorphic in the sample)", flush=True)

# Reindex ages/bins to match the filtered .bim, by position.
pos_to_age = pd.Series(ages, index=positions)
ages_filtered = pos_to_age.loc[bim["pos"].values].to_numpy()
bins_assigned_filtered = pd.cut(ages_filtered, bins=BINS,
                                 labels=bin_labels, right=False)


# =================================================================
# 12. Per-bin plink files (for GCTA / GRM-based diagnostics)
#     tskit writes site.id into the VCF ID column, and site.id is the
#     0-based index into the variant enumeration, so no offset is needed
#     when extracting.
# =================================================================
for lo, hi in zip(BINS[:-1], BINS[1:]):
    m = (ages >= lo) & (ages < hi)
    if m.sum() == 0:
        continue

    variant_ids = np.where(m)[0]

    label = f"{lo:.0f}_{hi:.0f}" if not np.isinf(hi) else f"{lo:.0f}_inf"
    var_file = SIM_PATH_REP / f"{SIM_VERSION}_bin_{label}_vars.txt"
    prefix = SIM_PATH_REP / f"{SIM_VERSION}_bin_{label}"

    np.savetxt(var_file, variant_ids, fmt="%d")

    cmd = (f"plink2 --vcf {VCF_FILE} --extract {var_file} "
           f"--make-bed --out {prefix}")
    subprocess.run(cmd, shell=True, check=True)

    print(f"Wrote {prefix}.bed/bim/fam with {m.sum()} variants", flush=True)
    os.remove(var_file)

# =================================================================
# 13. GENIE annotation matrix — built from the MAF-filtered .bim,
#     so row count and order match the genome-wide genotype file
#     GENIE will actually load.
# =================================================================
annotations = pd.get_dummies(bins_assigned_filtered, prefix="bin", prefix_sep="_").astype(int)

annotations.to_csv(
    SIM_PATH_REP / f"{SIM_VERSION}_annotations_age_bins.txt",
    sep=" ", index=False, header=False,
)

legend = pd.DataFrame({
    "column_name": annotations.columns,
    "age_bin": bin_labels,
})
legend.to_csv(
    SIM_PATH_REP / f"{SIM_VERSION}_annotations_legend.txt",
    sep=" ", index=False,
)


# =================================================================
# 14. Variant info and phenotype files for downstream analysis
# =================================================================
pd.DataFrame({
    "site_id": np.arange(M),
    "position": positions,
    "age": ages,
    "bin": bins_assigned,
    "freq": freqs,
    "beta": beta,
}).to_csv(SIM_PATH_REP / f"{SIM_VERSION}_variant_info.csv", index=False)

# FID/IID matching the VCF sample names, so GCTA can join on them.
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