"""
Full pipeline: SLiM tree sequence -> recapitate -> overlay mutations ->
extract true ages, frequencies, effect sizes -> phenotypes.

For a NEUTRAL trait simulation (single population, mutation-drift equilibrium).
"""

import tskit
import pyslim
import msprime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from pathlib import Path
import subprocess
import os


SIM_NE = 100000
SIM_MU = 1e-8
RNG_SEED = 42
RECOMBINATION_RATE = 1e-8

MU = 1.44e-8 # per base mutation rate 
PI_TARGET = 0.01 # fraction of the genome that's a mutational target

SIM_PATH="/Users/samm/Documents/Work/Projects/Vm/simulations/data"
SIM_VERSION="1.1"

np.random.default_rng(42)

# ---------------------------------------------------------------
# 1. Load the SLiM tree sequence
# ---------------------------------------------------------------
# Load in the SLiM treee
# We are doing this because it's easier to plop mutations onto the tree after we have simulated it forward
# rather than simulating it as you go along. 
TREE_FILE = Path(SIM_PATH) / f"{SIM_VERSION}_neutral_out.trees"
ts = tskit.load(TREE_FILE) 

# Because we simulated forward in time, some of the trees don't coalesce
# First check for multi-root trees
n_multiroot = sum(1 for t in ts.trees() if t.num_roots > 1)
print(f"Loaded: {ts.num_samples} samples, {ts.num_trees} trees, "
      f"{n_multiroot} not yet coalesced")

# ---------------------------------------------------------------
# 2. Recapitate (complete the ancient history with msprime coalescent)
#    ancestral_Ne should match your SLiM population size, since you
#    ran a constant-size single-population model.
# ---------------------------------------------------------------

## then recapacitate the trees
rts = pyslim.recapitate(
    ts,
    recombination_rate=RECOMBINATION_RATE,
    ancestral_Ne=SIM_NE,
    random_seed=RNG_SEED,
)

n_multiroot_after = sum(1 for t in rts.trees() if t.num_roots > 1)
print(f"After recapitation: {n_multiroot_after} not yet coalesced "
      f"(should be 0)")
assert n_multiroot_after == 0, "Recapitation incomplete - check ancestral_Ne"

# ---------------------------------------------------------------
# 3. Overlay neutral QTL mutations
#    mu = per-base mutation rate, pi_target = fraction of genome
#    that is mutational target (must match what you intend the
#    trait's mutational target to be)
# ---------------------------------------------------------------

## then simulate mutations over the top of the tree
## the rate, i.e. the number per generation is 
mts = msprime.sim_mutations(
    rts,
    rate=MU * PI_TARGET, # because we only want to simulate them at a rate of mutational target 
    random_seed=RNG_SEED,
    model=msprime.SLiMMutationModel(type=0),  # keeps SLiM-style metadata
    keep=True,   # keep any mutations already present (none, here)
)


print(f"After mutation overlay: {mts.num_sites} sites, "
      f"{mts.num_mutations} mutations")

DUMP_FILE = Path(SIM_PATH) / f"{SIM_VERSION}_neutral_out.recap.mut.trees"
mts.dump(DUMP_FILE)

# ---------------------------------------------------------------
# 4. Extract true ages, frequencies -> assign effect sizes
# ---------------------------------------------------------------
rng = np.random.default_rng(1)

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
print(f"Extracted {M} variants")
print(f"Age range: {ages.min():.1f} - {ages.max():.1f} generations")
print(f"Freq range: {freqs.min():.5f} - {freqs.max():.5f}")

# ---------------------------------------------------------------
# 5. Assign effect sizes (neutral trait: independent of age/freq)
# ---------------------------------------------------------------
sigma_beta = 1.0
beta = rng.normal(0, sigma_beta, size=M) # draw betas from gaussian distribution  

# ---------------------------------------------------------------
# 6. Ground-truth V_M and V_A checks
# ---------------------------------------------------------------
u_target = MU * PI_TARGET * (mts.sequence_length)   # per-gamete target rate: this is the number of mutations in the target area we expect to see per gamete
V_M_true = 2 * u_target * sigma_beta**2 # This is the total amount of additive variance put into the population from new mutations. 2 * u_target is the total number of new mutations
                                        # and then we multiple that by the expected effect size 
                                        # This is basically what we expect to see, not actually looking at any simulated sequence at all
V_A_analytic = 2 * SIM_NE * V_M_true # under mutation drift balance, with no selection on the trait - then standing additive variance is the input rate, V
                                     # multiplied by the persistance rate (think of it like a tape - the input rate is the rate of the tap and the persistence is the drain)
                                     # under pure drift, the persistence time is 2*Ne generations (small populations mean persistence time is lower and they get lost to drift faster)
                                     # So this value is what we expecte V_A to take before we do any simulations
V_A_empirical = np.sum(2 * freqs * (1 - freqs) * beta**2) # This is  the additive variation based on the sum(f * (1-f) * B) calculation. So using the data. 

print(f"\nV_M (true, analytic):   {V_M_true:.4g}")
print(f"V_A (analytic 2*Ne*V_M): {V_A_analytic:.4g}")
print(f"V_A (empirical, sum 2pq*beta^2): {V_A_empirical:.4g}")
print(f"Ratio empirical/analytic: {V_A_empirical/V_A_analytic:.3f}  "
      f"(should be close to 1 if at equilibrium)")

# ---------------------------------------------------------------
# 7. Sample individuals, build phenotypes
# ---------------------------------------------------------------
n_diploid_available = mts.num_samples // 2
n_sample = min(50000, n_diploid_available)
print(f"\nSampling {n_sample} of {n_diploid_available} available diploids")

# Build genotype matrix incrementally to avoid memory blowup
sample_idx = rng.choice(n_diploid_available, n_sample, replace=False)
hap_idx = np.sort(np.concatenate([2*sample_idx, 2*sample_idx + 1]))

g = np.zeros(n_sample)
for i, var in enumerate(mts.variants(samples=hap_idx)):
    geno = var.genotypes.reshape(-1, 2).sum(axis=1)  # diploid dosage
    g += beta[i] * geno

h2 = 0.5
V_E = np.var(g) * (1 - h2) / h2
y = g + rng.normal(0, np.sqrt(V_E), n_sample)

print(f"Empirical V_A in sample: {np.var(g):.4g}")
print(f"h2 (target): {h2}, V_E: {V_E:.4g}")

# ---------------------------------------------------------------
# 8. True per-bin variance (before touching a GRM or REML)
# ---------------------------------------------------------------
bins = [0, 25, 50, 100, 200, 500, 1000, 2000, 5000, 20000, np.inf]
print(f"\n{'Bin (gens)':<15}{'V_bin':<12}{'K_bar':<10}{'K_expected':<10}")
for lo, hi in zip(bins[:-1], bins[1:]):
    m = (ages >= lo) & (ages < hi)
    if m.sum() == 0:
        continue
    V_bin = np.sum(2 * freqs[m] * (1 - freqs[m]) * beta[m]**2) ## additive genetic variance from the bin (simply the effect size * f * (1-f))
    width = min(hi, 1e5) - lo # 
    K_bar = V_bin / (V_M_true * width) # K_bar is the proportion of the original V_M that exist per generation 
                                       # but we observe the values aggregated across the bin, so we have to divide 
                                       # by the width of the bin to get an approximate value
    K_integral_expected = 2*SIM_NE * (np.exp(-lo/(2*SIM_NE)) - np.exp(-hi/(2*SIM_NE))) 
    K_bar_expected = K_integral_expected / width   # mean of K over the bin
    K_expected = np.exp(-lo / (2 * SIM_NE))
    n_bin = m.sum()
    print(f"{lo}-{hi}: n={n_bin}, V={V_bin:.4g}, K̄={K_bar:.3f} ± {K_bar/np.sqrt(n_bin):.3f}, K_bar_expected={K_bar_expected:.3f}")

# ---------------------------------------------------------------
# 8. True per-bin variance (before touching a GRM or REML)
# ---------------------------------------------------------------

for lo, hi in zip(bins[:-1], bins[1:]):
    m = (ages >= lo) & (ages < hi)
    V_bin = np.sum(2*freqs[m]*(1-freqs[m])*beta[m]**2)
    K_int = 2*SIM_NE*(np.exp(-lo/(2*SIM_NE)) - (0 if np.isinf(hi) else np.exp(-hi/(2*SIM_NE))))
    print(f"{lo}-{hi}: n={m.sum()}, obs={V_bin:.4g}, pred={V_M_true*K_int:.4g}, "
          f"ratio={V_bin/(V_M_true*K_int):.3f}")


# ---------------------------------------------------------------
# 8. Write binned files to .bed format 
# ---------------------------------------------------------------

vcf_file = "all_variants.vcf.gz"
VCF_FILE = Path(SIM_PATH) / f"{SIM_VERSION}_neutral_out.vcf"

with gzip.open(vcf_file, 'wt') as f:
    mts.write_vcf(f)

for lo, hi in zip(bins[:-1], bins[1:]):
    m = (ages >= lo) & (ages < hi)
    variant_ids = np.where(m)[0]
    
    # Write variant IDs to file
    var_file = f"bin_{lo}_{hi}_vars.txt"
    np.savetxt(var_file, variant_ids + 1, fmt='%d')  # +1 because plink uses 1-based
    
    # Use plink to extract these variants and convert to .bed
    prefix = f"bin_{lo}_{hi}"
    cmd = f"plink2 --vcf {vcf_file} --extract {var_file} --make-bed --out {prefix}"
    subprocess.run(cmd, shell=True, check=True)
    
    print(f"Wrote {prefix}.bed/bim/fam with {m.sum()} variants")
    
    # Clean up temp file
    os.remove(var_file)


# ---------------------------------------------------------------
# 9. Write outputs for downstream R analysis
# ---------------------------------------------------------------
pd.DataFrame({
    "site_id": np.arange(M),
    "position": positions,
    "age": ages,
    "freq": freqs,
    "beta": beta,
}).to_csv("variant_info.csv", index=False)

pd.DataFrame({
    "sample_id": np.arange(n_sample),
    "y": y,
    "g": g,
}).to_csv("phenotypes.csv", index=False)

with open("out.vcf", "w") as f:
    mts.write_vcf(f, individuals=sample_idx)

print("\nWrote variant_info.csv, phenotypes.csv, out.vcf")
