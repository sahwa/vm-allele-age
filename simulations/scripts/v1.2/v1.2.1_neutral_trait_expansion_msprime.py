#!/usr/bin/env python
"""
v1.2.1 — neutral trait under Tennessen-style European expansion.

Simulates a neutral polygenic trait on an msprime tree sequence, bins variants
by true allele age, and writes a PLINK fileset + GENIE annotation matrix.

Memory note: the genotype matrix is never held in full. Variants are streamed in
genomic chunks; each chunk is packed straight into the .bed file and its
contribution to each individual's genetic value is accumulated on the fly.
Only per-variant metadata (O(M) floats) and per-individual vectors (O(N)) persist.

Difference to version v1.2 is that we will simulate 1Gb(!) of contiguous sequence
"""

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

# =================================================================
# Demography
# =================================================================
def scale_final_size(base_demog, pop, target_final, growth_rate):
    """Rebase a population's present-day size, keeping its growth rate."""
    d = copy.deepcopy(base_demog)
    d.add_population_parameters_change(
        time=0, population=pop, initial_size=target_final, growth_rate=growth_rate
    )
    return d


# =================================================================
# PLINK .bed streaming writer
# =================================================================
def open_bed(prefix):
    f = open(f"{prefix}.bed", "wb")
    f.write(bytes([0x6C, 0x1B, 0x01]))       # magic + SNP-major
    return f


def append_bed_chunk(f, dosages, n_ind, block_size=2000):
    """
    Append variants to an open SNP-major .bed.

    dosages : (n_var, n_ind) uint8, ALT allele count 0/1/2.
    A1 = ALT, A2 = REF (matches `plink2 --vcf` default orientation).
    Encoding: 00 hom-A1, 01 missing, 10 het, 11 hom-A2
              -> dosage 2 -> 0, dosage 1 -> 2, dosage 0 -> 3
    """
    n_var = dosages.shape[0]
    if n_var == 0:
        return
    n_bytes = (n_ind + 3) // 4
    pad = n_bytes * 4 - n_ind
    lut = np.array([3, 2, 0], dtype=np.uint8)

    for s in range(0, n_var, block_size):
        block = lut[dosages[s:s + block_size]]
        if pad:
            block = np.pad(block, ((0, 0), (0, pad)))
        block = block.reshape(block.shape[0], n_bytes, 4)
        packed = (
            block[:, :, 0]
            | (block[:, :, 1] << 2)
            | (block[:, :, 2] << 4)
            | (block[:, :, 3] << 6)
        ).astype(np.uint8)
        packed.tofile(f)


def write_bim_fam(prefix, positions, ref, alt, iids, chrom=1):
    snp_ids = [f"{chrom}:{int(p)}:{r}:{a}" for p, r, a in zip(positions, ref, alt)]
    pd.DataFrame({
        "chr": chrom, "snpid": snp_ids, "cm": 0,
        "pos": positions.astype(np.int64), "a1": alt, "a2": ref,
    }).to_csv(f"{prefix}.bim", sep="\t", index=False, header=False)

    pd.DataFrame({
        "fid": 0, "iid": iids, "pid": 0, "mid": 0, "sex": 0, "pheno": -9,
    }).to_csv(f"{prefix}.fam", sep="\t", index=False, header=False)


# =================================================================
# 1. Simulate (or load cached)
# =================================================================
if TREE_FILE.is_file():
    print(f"Loading cached tree sequence: {TREE_FILE}", flush=True)
    ts = tskit.load(str(TREE_FILE))
    print(f"Successfully loaded cached tree sequence: {TREE_FILE}", flush=True)
else:
    species = stdpopsim.get_species("HomSap")
    base = species.get_demographic_model(DEMOG_MODEL).model
    demog = scale_final_size(base, POP, N_FINAL_TARGET, EUR_GROWTH_RATE)
    demog.sort_events()

    print(f"Simulating ancestry: {N_SAMPLE} diploids, {L/1e6:.0f} Mb...", flush=True)
    ts = msprime.sim_ancestry(
        samples={POP: N_SAMPLE},
        demography=demog,
        sequence_length=L,
        recombination_rate=REC,
        model=[
            msprime.DiscreteTimeWrightFisher(duration=100),
            msprime.StandardCoalescent(),
        ],
        random_seed=RNG_SEED,
    )
    print("Overlaying mutations...", flush=True)
    ts = msprime.sim_mutations(ts, rate=MU, random_seed=RNG_SEED)
    ts.dump(str(TREE_FILE))
    print(f"Saved: {TREE_FILE}", flush=True)

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

# # =================================================================
# # 2. Effect sizes (neutral trait: beta independent of age and frequency)
# # =================================================================
causal = rng_causal.random(M_RAW) < PI_CAUSAL
beta = np.zeros(M_RAW)
beta[causal] = rng_beta.normal(0, SIGMA_BETA, size=causal.sum())
print(f"{causal.sum()} causal variants ({causal.mean():.3%})", flush=True)

# # =================================================================
# # 3. Streaming pass: pack .bed, accumulate genetic values, collect metadata
# # =================================================================

tabs = ts.tables

mut_node = tabs.mutations.node
mut_time = tabs.mutations.time
ages_all = np.where(np.isnan(mut_time), tabs.nodes.time[mut_node], mut_time)

is_sample = np.zeros(ts.num_nodes, dtype=bool)
is_sample[ts.samples()] = True
sing = is_sample[mut_node]

carrier   = tabs.nodes.individual[mut_node[sing]]   # node -> individual
sing_ages = ages_all[sing]

YOUNG_MAX = 205.0        # 2T12 growth onset; older singletons are a separate regime


def count_by_bin(bins):
    """Per-individual singleton counts for a given set of bin edges."""
    idx = np.searchsorted(bins, sing_ages, side="right") - 1
    nb = len(bins) - 1
    counts = np.bincount(
        carrier * nb + idx, minlength=N_IND * nb
    ).reshape(N_IND, nb)
    labels = [
        f"{lo:.0f}+" if np.isinf(hi) else f"{lo:.0f}-{hi:.0f}"
        for lo, hi in zip(bins[:-1], bins[1:])
    ]
    return pd.DataFrame(counts, columns=labels)


a = sing_ages[sing_ages < YOUNG_MAX]

for nb in [1, 2, 3, 4, 6]:
    q = np.quantile(a, np.linspace(0, 1, nb + 1))
    q[0], q[-1] = 0.0, YOUNG_MAX
    bins = np.concatenate([q, [np.inf]])          # + one catch-all for old

    df = count_by_bin(bins)
    out = SIM_PATH_REP / f"{SIM_VERSION}_singleton_counts_nb{nb}.csv"
    df.to_csv(out, index=False)

    print(f"\nnb={nb}  edges={np.round(q, 1)}  -> {out.name}")
    print(f"{'Bin':<18}{'n_sing':>10}{'mean':>9}{'n_eff':>9}{'frac0':>8}")
    for lab in df.columns:
        c = df[lab].to_numpy(float)
        if c.sum() == 0:
            continue
        n_eff = c.sum() ** 2 / (c ** 2).sum()
        print(f"{lab:<18}{int(c.sum()):>10}{c.mean():>9.1f}"
              f"{n_eff:>9.0f}{(c == 0).mean():>8.3f}")


# # =================================================================
# # 4. Phenotype: y = g + e, scaled to the target heritability
# # =================================================================

N_COMMON = 20_000     # we want to use N common variants
SIGMA_B_SING = 1.0    # same sigma for the betas as before
H2_SING = 0.10        # variance share from singletons
H2_COMMON = 0.40      # variance share from common variants
N_COMMON = 20_000
OVERSAMPLE = 50          # tune once you know the common fraction

# --- singleton contribution (from the table pass, no genotypes) ---
beta_sing = rng_beta.normal(0, SIGMA_B_SING, size=sing.sum())
g_sing = np.bincount(carrier, weights=beta_sing, minlength=N_IND)

# --- pick common variants and pull only those columns ---
cand = rng_causal.choice(ts.num_sites, size=N_COMMON * OVERSAMPLE, replace=False)
cand.sort()

var = tskit.Variant(ts)
i=0
keep_sites, keep_gt = [], []
for site_id in cand:
    var.decode(int(site_id))
    d = var.genotypes.reshape(-1, 2).sum(axis=1)
    p = d.mean() / 2
    if 0.05 < p < 0.95:
        i+=1
        if i % 100 == 0:
            print(i)
        keep_sites.append(site_id)
        keep_gt.append(d.astype(np.uint8))
        if len(keep_sites) == N_COMMON:
            break

Gc = np.array(keep_gt)
picked = np.array(keep_sites)

p = Gc.mean(axis=1) / 2
Z = (Gc - 2 * p[:, None]) / np.sqrt(2 * p * (1 - p))[:, None]
beta_common = rng_beta.normal(0, 1, size=len(picked))
g_common = Z.T @ beta_common

g_sing_bin = np.zeros((N_IND, N_BINS))
for b in range(N_BINS):
    m = bin_idx == b
    if m.any():
        g_sing_bin[:, b] = np.bincount(
            carrier[m], weights=beta_sing[m], minlength=N_IND)

g_sing = g_sing_bin.sum(axis=1)

# --- scale each to its target share, then add noise ---
scale_sing = np.sqrt(H2_SING / g_sing.var())
g_sing_bin *= scale_sing
g_sing     *= scale_sing

g_common *= np.sqrt(H2_COMMON / g_common.var())
g = g_sing + g_common
g -= g.mean()
y = g + rng_noise.normal(0, np.sqrt(1 - H2_SING - H2_COMMON), size=N_IND)


print(f"var(g_sing)={g_sing.var():.4f}  var(g_common)={g_common.var():.4f}  "
      f"var(y)={y.var():.4f}")

g -= g.mean()
var_g = g.var()
sigma_e = np.sqrt(var_g * (1 - H2) / H2)
y = g + rng_noise.normal(0, sigma_e, size=N_IND)

print(f"\nvar(g)={var_g:.4g}  sigma_e={sigma_e:.4g}  "
      f"realised h2={var_g / y.var():.3f}", flush=True)

# # =================================================================
# # 5. Ground truth: V_M and V_A
# # =================================================================
bin_labels = [
    f"{lo:.0f}+" if np.isinf(hi) else f"{lo:.0f}-{hi:.0f}"
    for lo, hi in zip(BINS[:-1], BINS[1:])
]

print(f"\n{'Bin':<20}{'n_sing':>10}{'V_true':>12}{'share':>9}")
for b, label in enumerate(bin_labels):
    if sing_counts[:, b].sum() == 0:
        continue
    V_b = g_sing_bin[:, b].var()
    print(f"{label:<20}{sing_counts[:, b].sum():>10}{V_b:>12.5f}"
          f"{V_b / H2_SING:>9.3f}")

print(f"\nvar(g_sing)  = {g_sing.var():.4f}  (target {H2_SING})")
print(f"var(g_common)= {g_common.var():.4f}  (target {H2_COMMON})")
print(f"var(y)       = {y.var():.4f}")


# # =================================================================
# # 6. Singleton Counts
# # =================================================================

young = ages_all[sing] < 205
a = ages_all[sing][young]

for nb in [2, 3, 4]:
    q = np.quantile(a, np.linspace(0, 1, nb + 1))
    q[0], q[-1] = 0, 205
    print(nb, np.round(q, 1))

BINS = np.concatenate([q, [np.inf]])


# # =================================================================
# # 7. True per-bin genetic variance
# # =================================================================

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

iids = [f"ind{i}" for i in range(N_IND)]

pd.DataFrame({"FID": 0, "IID": iids, "PHENO": y}).to_csv(
    SIM_PATH_REP / f"{SIM_VERSION}_phenotypes.GENIE.txt",
    sep="\t", index=False,
)
pd.DataFrame({"iid": iids, "y": y, "g": g}).to_csv(
    SIM_PATH_REP / f"{SIM_VERSION}_phenotypes.csv", index=False)

print(f"\nWrote outputs to {SIM_PATH_REP}", flush=True)
