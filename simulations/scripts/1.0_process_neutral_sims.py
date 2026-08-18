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


SIM_NE = 20000
SIM_MU = 1e-8
RNG_SEED = 42

MU = 1.44e-8 # per base mutation rate 
PI_TARGET = 0.01 # fraction of the genome that's a mutational target


# ---------------------------------------------------------------
# 1. Load the SLiM tree sequence
# ---------------------------------------------------------------
# Load in the SLiM treee
# We are doing this because it's easier to plop mutations onto the tree after we have simulated it forward
# rather than simulating it as you go along. 
ts = tskit.load("test.trees") 

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
    recombination_rate=SIM_MU,
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
    random_seed=2,
    model=msprime.SLiMMutationModel(type=0),  # keeps SLiM-style metadata
    keep=True,   # keep any mutations already present (none, here)
)

n_multiroot = sum(1 for t in mts.trees() if t.num_roots > 1)
print(f"Loaded: {mts.num_samples} samples, {mts.num_trees} trees, "
      f"{n_multiroot} not yet coalesced")


print(f"After mutation overlay: {mts.num_sites} sites, "
      f"{mts.num_mutations} mutations")

mts.dump("out.recap.mut.trees")

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


SIM_NE = 20000
V_M_true = 0.0288

bin_lo = np.array([0, 25, 50, 100, 200, 500, 1000, 2000, 5000, 20000])
bin_hi = np.array([25, 50, 100, 200, 500, 1000, 2000, 5000, 20000, np.inf])

# ----------------------------------------------------------------------
# Compute from your data. Requires `ages`, `freqs`, `beta` in scope.
# Set FROM_DATA = False to use the hard-coded values from your run.
# ----------------------------------------------------------------------
FROM_DATA = True

if FROM_DATA:
    V_obs, n_eff, n_bin = [], [], []
    for lo, hi in zip(bin_lo, bin_hi):
        m = (ages >= lo) & (ages < hi)
        contribs = 2 * freqs[m] * (1 - freqs[m]) * beta[m] ** 2
        V_obs.append(contribs.sum())
        # Kish effective n: variance is dominated by a few large contributors
        n_eff.append(contribs.sum() ** 2 / (contribs ** 2).sum())
        n_bin.append(m.sum())
    V_obs = np.array(V_obs)
    n_eff = np.array(n_eff)
    n_bin = np.array(n_bin)
else:
    V_obs = np.array([0.6957, 0.8032, 1.477, 2.614, 7.258, 13.4,
                      31.62, 93.76, 319.0, 686.1])
    n_bin = np.array([3107, 716, 758, 759, 989, 798, 853, 1098, 1656, 2543])
    n_eff = np.array([40, 25, 30, 35, 55, 60, 80, 130, 300, 700])

# ----------------------------------------------------------------------
# Analytic expectation: V_pred = V_M * integral of exp(-t/2Ne) over the bin
# ----------------------------------------------------------------------
def kernel_integral(lo, hi, SIM_NE):
    hi_term = 0.0 if np.isinf(hi) else np.exp(-hi / (2 * SIM_NE))
    return 2 * SIM_NE * (np.exp(-lo / (2 * SIM_NE)) - hi_term)

K_int = np.array([kernel_integral(lo, hi, SIM_NE)
                  for lo, hi in zip(bin_lo, bin_hi)])
V_pred = V_M_true * K_int

width = np.where(np.isinf(bin_hi), 5 * SIM_NE - bin_lo, bin_hi - bin_lo)
K_bar_obs = V_obs / (V_M_true * width)
K_bar_exp = K_int / width
K_bar_se = K_bar_obs / np.sqrt(n_eff)

# geometric bin midpoints for the x axis
x_mid = np.where(np.isinf(bin_hi), bin_lo * 2.0,
                 np.sqrt(bin_lo.clip(min=1) *
                         np.where(np.isinf(bin_hi), 1, bin_hi)))
x_mid[0] = 12.5

# ----------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
})

INK, ACCENT, MUTED, FLAT = "#1a1a1a", "#c1440e", "#8a8a8a", "#e8e4dc"

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

# ---------------- Panel A ----------------
ax = axes[0]
lims = [0.4, 1500]
ax.plot(lims, lims, color=MUTED, lw=1.0, ls="--", zorder=1, label="1:1")
ax.scatter(V_pred, V_obs, s=52, color=ACCENT, zorder=3,
           edgecolor="white", linewidth=0.8)

for xp, yo, lo, hi in zip(V_pred, V_obs, bin_lo, bin_hi):
    lab = f"{lo:.0f}+" if np.isinf(hi) else f"{lo:.0f}\u2013{hi:.0f}"
    ax.annotate(lab, (xp, yo), textcoords="offset points",
                xytext=(7, -9), fontsize=7, color="#6a6a6a")

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(lims); ax.set_ylim(lims); ax.set_aspect("equal")
ax.set_xlabel(r"predicted $V_{\rm bin}$   ($V_M \int_{\rm bin} e^{-t/2N_e}dt$)")
ax.set_ylabel(r"observed $V_{\rm bin}$")
ax.set_title("A   Per-bin variance matches the neutral kernel",
             loc="left", fontsize=10.5, pad=10)
ax.legend(frameon=False, fontsize=8.5, loc="upper left")
ax.xaxis.set_minor_formatter(NullFormatter())
ax.yaxis.set_minor_formatter(NullFormatter())

# ---------------- Panel B ----------------
ax = axes[1]
ax.axhspan(0.95, 1.05, color=FLAT, zorder=0)
ax.text(1.5, 1.055, r"$K > 0.95$: $V_M$ recoverable without correction",
        fontsize=8, color="#6a6a6a", va="bottom")

t_smooth = np.logspace(0, np.log10(5 * SIM_NE), 400)
ax.plot(t_smooth, np.exp(-t_smooth / (2 * SIM_NE)), color=MUTED, lw=1.0,
        ls="--", zorder=1, label=r"$K(t)=e^{-t/2N_e}$")

ax.errorbar(x_mid, K_bar_obs, yerr=K_bar_se, fmt="o", ms=6.5,
            color=ACCENT, ecolor=ACCENT, elinewidth=1.1, capsize=2.5,
            markeredgecolor="white", markeredgewidth=0.8,
            zorder=3, label=r"observed $\bar{K}$")
ax.scatter(x_mid, K_bar_exp, s=46, facecolor="none", edgecolor=INK,
           linewidth=1.1, zorder=2, label=r"expected $\bar{K}$")

ax.axhline(1.0, color=INK, lw=0.6, ls=":", zorder=1)
ax.set_xscale("log"); ax.set_ylim(0, 1.35)
ax.set_xlabel("allele age (generations)")
ax.set_ylabel(r"surviving fraction of $V_M$,  $\bar{K}$")
ax.set_title("B   Kernel is flat across the young bins",
             loc="left", fontsize=10.5, pad=10)
ax.legend(frameon=False, fontsize=8.5, loc="lower left")
ax.xaxis.set_minor_formatter(NullFormatter())

fig.suptitle(r"Neutral trait, $N_e=20{,}000$: observed vs analytic "
             r"survival kernel (no free parameters)",
             y=1.0, fontsize=11)
fig.tight_layout()
fig.savefig("neutral_validation.png", dpi=200, bbox_inches="tight")
fig.savefig("neutral_validation.pdf", bbox_inches="tight")


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