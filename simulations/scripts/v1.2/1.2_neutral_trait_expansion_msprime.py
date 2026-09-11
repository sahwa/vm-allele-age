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
H2 = 0.5                      # target narrow-sense heritability
SIGMA_BETA = 1.0              # SD of causal effect sizes
PI_CAUSAL = 0.01              # fraction of segregating sites that are causal

# Age bins (log-spaced; edges chosen to span the observed age range)
BINS = np.concatenate([[0.0], np.geomspace(10, 1e6, 9), [np.inf]])

# Streaming
CHUNK_BP = 2e6                # genomic window per chunk; tune to available RAM

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
else:
    species = stdpopsim.get_species("HomSap")
    base = species.get_demographic_model(DEMOG_MODEL).model
    demog = scale_final_size(base, POP, N_FINAL_TARGET, EUR_GROWTH_RATE)

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

# =================================================================
# 2. Effect sizes (neutral trait: beta independent of age and frequency)
# =================================================================
causal = rng_causal.random(M_RAW) < PI_CAUSAL
beta = np.zeros(M_RAW)
beta[causal] = rng_beta.normal(0, SIGMA_BETA, size=causal.sum())
print(f"{causal.sum()} causal variants ({causal.mean():.3%})", flush=True)

# =================================================================
# 3. Streaming pass: pack .bed, accumulate genetic values, collect metadata
# =================================================================
ages = np.empty(M_RAW)
freqs = np.empty(M_RAW)
positions = np.empty(M_RAW)
ref = np.empty(M_RAW, dtype="<U1")
alt = np.empty(M_RAW, dtype="<U1")
keep = np.zeros(M_RAW, dtype=bool)

g = np.zeros(N_IND)                      # genetic value, accumulated per chunk
bed_f = open_bed(GENOME_PREFIX)

print(f"\nStreaming {M_RAW} sites in {CHUNK_BP/1e6:.0f} Mb chunks...", flush=True)
offset, start = 0, 0.0
while start < L:
    end = min(start + CHUNK_BP, L)
    sub = ts.keep_intervals([[start, end]], simplify=False)
    n_c = sub.num_sites
    if n_c == 0:
        start = end
        continue
    sl = slice(offset, offset + n_c)

    G = sub.genotype_matrix()                       # (n_c, 2*N_IND) int32
    dosage = G.reshape(n_c, -1, 2).sum(axis=2).astype(np.uint8)
    freqs[sl] = G.mean(axis=1)
    del G

    tabs = sub.tables
    mt = tabs.mutations.time
    nt = tabs.nodes.time[tabs.mutations.node]
    ages[sl] = np.where(np.isnan(mt), nt, mt)
    positions[sl] = tabs.sites.position
    ref[sl] = tskit.unpack_strings(
        tabs.sites.ancestral_state, tabs.sites.ancestral_state_offset)
    alt[sl] = tskit.unpack_strings(
        tabs.mutations.derived_state, tabs.mutations.derived_state_offset)

    # Drop variants monomorphic in the sample; GENIE requires MAF > 0.
    k = (freqs[sl] > 0) & (freqs[sl] < 1)
    keep[sl] = k

    append_bed_chunk(bed_f, dosage[k], N_IND)
    g += dosage[k].T.astype(np.float64) @ beta[sl][k]

    print(f"  {start/1e6:6.0f}-{end/1e6:<6.0f} Mb  {n_c:7d} sites  "
          f"{k.sum():7d} kept", flush=True)

    del dosage, tabs, sub, k
    offset += n_c
    start = end

bed_f.close()
assert offset == M_RAW, f"site count mismatch: {offset} != {M_RAW}"

iids = [f"ind{i}" for i in range(N_IND)]
write_bim_fam(GENOME_PREFIX, positions[keep], ref[keep], alt[keep], iids)

print(f"\nAge range:  {ages.min():.1f} - {ages.max():.1f} generations")
print(f"Freq range: {freqs.min():.3g} - {freqs.max():.3g}")
print(f"Retained:   {keep.sum()} of {M_RAW} "
      f"({M_RAW - keep.sum()} monomorphic in sample)", flush=True)


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