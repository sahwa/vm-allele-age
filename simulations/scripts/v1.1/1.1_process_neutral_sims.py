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

import msprime
import numpy as np
import pandas as pd
import pyslim
import tskit

# =================================================================
# 0. Parameters
# =================================================================

# %

REP = int(sys.argv[1])

SIM_NE = 100000
RNG_SEED = 42

ss = np.random.SeedSequence([RNG_SEED, REP])
seed_mutations, seed_beta, seed_sample, seed_noise = ss.spawn(4)

rng_beta = np.random.default_rng(seed_beta)
rng_sample = np.random.default_rng(seed_sample)
rng_noise = np.random.default_rng(seed_noise)

RECOMBINATION_RATE = 1e-8   # must match initializeRecombinationRate() in the .slim script

MU = 1.44e-8                # per base per generation mutation rate
PI_TARGET = 0.01             # fraction of the genome that is mutational target

SIM_PATH = Path(
    "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.1"
)
SIM_VERSION = "1.1"

SIM_PATH_REP = SIM_PATH / f"rep{REP}"
SIM_PATH_REP.mkdir(exist_ok=True)

N_SAMPLE_TARGET = 50000      # diploids to sample for the GREML/GENIE analysis
H2 = 0.5
SIGMA_BETA = 1.0

BINS = [0, 1e2, 1e3, 1e4, 5e4, 1e5, 2e5, 5e5, np.inf]

TREE_FILE = SIM_PATH / f"{SIM_VERSION}_neutral_trait.n_100000.trees"  # shared, not per-rep
MUT_TREE_FILE = SIM_PATH_REP / f"{SIM_VERSION}_neutral_out.recap.mut.trees"
VCF_FILE = SIM_PATH_REP / f"{SIM_VERSION}_neutral_out.vcf"
GENOME_PREFIX = SIM_PATH_REP / f"{SIM_VERSION}_neutral_out"

# =================================================================
# 1. Load the SLiM tree sequence
# =================================================================
# Mutations are overlaid *after* the forward simulation rather than during it,
# because the trait is neutral: nothing in the forward sim depends on where
# the mutations landed, so msprime can paint them on afterwards far more
# cheaply than tracking them in SLiM itself.
# ts = tskit.load(TREE_FILE)

# n_multiroot = sum(1 for t in ts.trees() if t.num_roots > 1)
# print(f"Loaded: {ts.num_samples} samples, {ts.num_trees} trees, "
#       f"{n_multiroot} not yet coalesced", flush=True)

# total_span = 0
# weighted_tmrca = 0
# for tree in ts.trees():
#     root_time = tree.time(tree.root) if tree.num_roots == 1 else None
#     if root_time is not None:
#         weighted_tmrca += root_time * tree.span
#         total_span += tree.span

# mean_tmrca = weighted_tmrca / total_span
# print(f"Mean TMRCA: {mean_tmrca:.0f} (expect ~4*Ne = {4*SIM_NE} under neutrality)")

# =================================================================
# 2. Recapitate (complete the ancient history with the msprime coalescent)
#    ancestral_Ne must match the SLiM population size.
# =================================================================

RTS_FILE = SIM_PATH / f"{SIM_VERSION}_neutral_trait.n_100000.recapitated.trees"

if not RTS_FILE.exists():
    ts = tskit.load(TREE_FILE)
    rts = pyslim.recapitate(ts, recombination_rate=RECOMBINATION_RATE,
                            ancestral_Ne=SIM_NE, random_seed=RNG_SEED)
    rts.dump(RTS_FILE)
    n_multiroot_after = sum(1 for t in rts.trees() if t.num_roots > 1)
    print(f"After recapitation: {n_multiroot_after} not yet coalesced (should be 0)", flush=True)
    assert n_multiroot_after == 0, "Recapitation incomplete - check ancestral_Ne"
else:
    rts = tskit.load(RTS_FILE)

# =================================================================
# 3. Overlay neutral QTL mutations
#    The pi-scaling trick: rather than designating 1% of sites as targets and
#    mutating at the full rate, mutate every site at 1% of the rate. Same
#    expected count of trait mutations, no need to track which sites are
#    targets. JC69 gives real A/C/G/T alleles, which plink and later dating
#    tools need.
# =================================================================

print("Overlaying mutations...", flush=True)
mts = msprime.sim_mutations(
    rts, rate=MU * PI_TARGET,
    random_seed=seed_mutations.generate_state(1)[0],  # SeedSequence -> int
    model=msprime.JC69(), keep=True,
)

print(f"After mutation overlay: {mts.num_sites} sites, "
      f"{mts.num_mutations} mutations", flush=True)

# =================================================================
# 4. Strip multiallelic sites directly from mts.
#    Doing this once, here, means every downstream .variants() loop,
#    the VCF write, and the plink/GENIE outputs are automatically in
#    lockstep - no separate masking or post-hoc row-matching needed.
# =================================================================
multiallelic_site_ids = np.array(
    [site.id for site in mts.sites() if len(site.mutations) > 1]
)
print(f"Removing {len(multiallelic_site_ids)} multiallelic sites "
      f"out of {mts.num_sites}", flush=True)

tables = mts.dump_tables()
tables.delete_sites(multiallelic_site_ids)
tables.sort()
mts = tables.tree_sequence()
print(f"mts now has {mts.num_sites} biallelic sites", flush=True)

mts.dump(MUT_TREE_FILE)

# =================================================================
# 5. Extract true ages, frequencies, positions
#    mts is guaranteed biallelic at this point, so no per-variant
#    filtering is needed here.
# =================================================================

mts = tskit.load(MUT_TREE_FILE)

ages, freqs, positions = [], [], []
for var in mts.variants():
    mut = mts.mutation(var.site.mutations[0].id)
    age = mut.time if not np.isnan(mut.time) else mts.node(mut.node).time
    ages.append(age)
    freqs.append(var.genotypes.mean())
    positions.append(var.site.position)

ages = np.array(ages)
freqs = np.array(freqs)
positions = np.array(positions)
M = len(ages)
print(f"Extracted {M} variants", flush=True)
print(f"Age range: {ages.min():.1f} - {ages.max():.1f} generations")
print(f"Freq range: {freqs.min():.5f} - {freqs.max():.5f}")

# =================================================================
# 6. Assign effect sizes (neutral trait: independent of age and frequency)
# =================================================================
beta = rng_beta.normal(0, SIGMA_BETA, size=M)

# =================================================================
# 7. Ground-truth V_M and V_A checks
# =================================================================
# u_target: expected number of target mutations per gamete per generation.
u_target = MU * PI_TARGET * mts.sequence_length

# V_M = 2 * U * E[beta^2]. The 2 is diploidy; E[beta^2] = SIGMA_BETA^2 for a
# mean-zero normal. Computed purely from parameters, no data involved.
V_M_true = 2 * u_target * SIGMA_BETA**2

# Under mutation-drift balance with no selection, standing V_A is the input
# rate times the persistence time, and under pure drift persistence is 2*Ne.
V_A_analytic = 2 * SIM_NE * V_M_true

# The same quantity measured from the realised data.
V_A_empirical = np.sum(2 * freqs * (1 - freqs) * beta**2)

print(f"\nV_M (true, analytic):            {V_M_true:.4g}")
print(f"V_A (analytic 2*Ne*V_M):         {V_A_analytic:.4g}")
print(f"V_A (empirical, sum 2pq*beta^2): {V_A_empirical:.4g}")
print(f"Ratio empirical/analytic: {V_A_empirical / V_A_analytic:.3f}  "
      f"(should be close to 1 if at equilibrium)", flush=True)

# =================================================================
# 8. Build genetic values, then sample individuals
#
#    NOTE: build g over ALL individuals first, then subset. Passing
#    samples=hap_idx to .variants() drops variants that are monomorphic
#    in the subset, which desynchronises the variant enumeration from
#    beta[i] and silently corrupts g.
# =================================================================
n_dip_all = mts.num_samples // 2
print(f"\nBuilding genetic values for all {n_dip_all} diploids...", flush=True)
print(mts.num_individuals, n_dip_all)  # do these match?

g_all = np.zeros(n_dip_all)
for i, var in enumerate(mts.variants()):
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
    hi_term = 0.0 if np.isinf(hi) else np.exp(-hi / (2 * SIM_NE))
    K_int = 2 * SIM_NE * (np.exp(-lo / (2 * SIM_NE)) - hi_term)
    V_pred = V_M_true * K_int

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
cmd = (f"plink2 --vcf {VCF_FILE} --maf 0.0000001 --set-all-var-ids '@:#:$r:$a' "
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
# 12. Per-bin plink files, built from the MAF-filtered genome-wide
#     .bed - so variant sets and IDs stay consistent with GENOME_PREFIX
#     and with snp_info.csv (all keyed by chr:pos:ref:alt).
# =================================================================
bins_assigned_filtered = pd.Series(bins_assigned_filtered, index=bim["snpid"].values)

for label_str in bin_labels:
    snp_ids_in_bin = bins_assigned_filtered[bins_assigned_filtered == label_str].index

    if len(snp_ids_in_bin) == 0:
        continue

    var_file = SIM_PATH_REP / f"{SIM_VERSION}_bin_{label_str}_vars.txt"
    prefix = SIM_PATH_REP / f"{SIM_VERSION}_bin_{label_str}"

    pd.Series(snp_ids_in_bin).to_csv(var_file, index=False, header=False)

    cmd = (f"plink2 --bfile {GENOME_PREFIX} --extract {var_file} "
           f"--make-bed --out {prefix}")
    subprocess.run(cmd, shell=True, check=True)

    print(f"Wrote {prefix}.bed/bim/fam with {len(snp_ids_in_bin)} variants",
          flush=True)
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


bim = pd.read_csv(f"{GENOME_PREFIX}.bim", sep="\t", header=None,
                   names=["chr", "snpid", "cm", "pos", "a1", "a2"])

# reindex ages/freqs to match, by position (as you're already doing)
pos_to_age = pd.Series(ages, index=positions)
ages_filtered = pos_to_age.loc[bim["pos"].values].to_numpy()

snp_info = pd.DataFrame({"snp_id": bim["snpid"].values})  # <- use plink's own IDs

pos_to_freq = pd.Series(freqs, index=positions)
freqs_filtered = pos_to_freq.loc[bim["pos"].values].to_numpy()
p = freqs_filtered

bins_assigned_clean = pd.cut(ages_filtered, bins=BINS, labels=bin_labels, right=False)

annotations = pd.get_dummies(bins_assigned_clean, prefix="bin", prefix_sep="_").astype(int)
annotations.index = snp_info.index  # ensure positional alignment, not label alignment

for col in annotations.columns:
    snp_info[col] = annotations[col].replace(0, np.nan).values  # .values forces positional assignment

snp_info["w_unstd"] = 1.0
snp_info["w_std"] = 1.0 / np.sqrt(2 * p * (1 - p))

snp_info["c0"] = -2 * p       # A1A1: -2p
snp_info["c1"] = 1 - 2 * p    # A1A2: 1-2p  
snp_info["c2"] = 2 - 2 * p    # A2A2: 2-2p

snp_info.to_csv(SIM_PATH_REP / f"{SIM_VERSION}_snp_info_mph.csv", index=False)


print(f"\nWrote outputs to {SIM_PATH_REP}", flush=True)