"""
Full pipeline: SLiM tree sequence -> recapitate -> overlay mutations ->
extract true ages, frequencies, effect sizes -> phenotypes -> plink fileset
-> GENIE annotation matrix.

For a NEUTRAL trait simulation (single population, mutation-drift equilibrium).

Changes from the previous version:
  - ONE pass over mts.variants() instead of three. That traversal dominates
    runtime; everything downstream is derived from arrays built in that pass.
  - .bed/.bim/.fam written directly with numpy. No VCF, no plink subprocess.
  - g computed from the sampled dosage matrix, not from a full-population
    g_all that was then 50% discarded.
  - MAF filtering applied to ages/beta/freqs/positions simultaneously, so the
    annotation aligns to the .bim by construction. No pos_to_age.loc[bim.pos]
    round-trip, which has been a recurring source of silent misalignment.
"""

import sys
from pathlib import Path

import msprime
import numpy as np
import pandas as pd
import pyslim
import tskit


# =================================================================
# 0. Parameters
# =================================================================
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
PI_TARGET = 0.01            # fraction of the genome that is mutational target

SIM_PATH = Path(
    "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.1"
)
SIM_VERSION = "1.1"

SIM_PATH_REP = SIM_PATH / "replicates" / f"rep{REP}"
SIM_PATH_REP.mkdir(parents=True, exist_ok=True)

N_SAMPLE_TARGET = 50000     # diploids to sample for the GREML/GENIE analysis
H2 = 0.5
SIGMA_BETA = 1.0

BINS = [0, 1e2, 1e3, 1e4, 5e4, 1e5, 2e5, 5e5, np.inf]

# The mutation-overlaid tree sequence is large and only needed if you want to
# re-analyse this replicate later (e.g. tsdate). Off by default for big runs.
SAVE_MUT_TREES = False

TREE_FILE = SIM_PATH / f"{SIM_VERSION}_neutral_trait.n_100000.trees"   # shared
RTS_FILE = SIM_PATH / f"{SIM_VERSION}_neutral_trait.n_100000.recapitated.trees"
MUT_TREE_FILE = SIM_PATH_REP / f"{SIM_VERSION}_neutral_out.recap.mut.trees"
GENOME_PREFIX = SIM_PATH_REP / f"{SIM_VERSION}_neutral_out"


# =================================================================
# Helper: write a PLINK1 .bed/.bim/.fam fileset directly
# =================================================================
def write_plink_bed(dosages, prefix, positions, ref, alt, iids, chrom=1):
    """
    dosages : (n_var, n_ind) uint8, count of the ALT (derived) allele, 0/1/2
    prefix  : output path prefix

    A1 is set to ALT and A2 to REF, matching what `plink2 --vcf` does by
    default, so allele orientation is unchanged from the old VCF route.

    PLINK1 .bed: 3 magic bytes, then 2 bits per genotype, 4 per byte, with the
    first individual in the lowest-order bits.
        00 = hom A1, 01 = missing, 10 = het, 11 = hom A2
    With A1 = ALT: dosage 2 -> 00 (0), dosage 1 -> 10 (2), dosage 0 -> 11 (3)
    """
    n_var, n_ind = dosages.shape
    n_bytes = (n_ind + 3) // 4
    pad = n_bytes * 4 - n_ind

    lut = np.array([3, 2, 0], dtype=np.uint8)   # indexed by dosage

    with open(f"{prefix}.bed", "wb") as f:
        f.write(bytes([0x6C, 0x1B, 0x01]))      # magic + SNP-major
        for start in range(0, n_var, 2000):
            block = lut[dosages[start:start + 2000]]
            if pad:
                block = np.pad(block, ((0, 0), (0, pad)))
            block = block.reshape(block.shape[0], n_bytes, 4)
            packed = (block[:, :, 0]
                      | (block[:, :, 1] << 2)
                      | (block[:, :, 2] << 4)
                      | (block[:, :, 3] << 6)).astype(np.uint8)
            packed.tofile(f)

    # IDs match the old `--set-all-var-ids '@:#:$r:$a'` scheme
    snp_ids = [f"{chrom}:{int(p)}:{r}:{a}" for p, r, a in zip(positions, ref, alt)]
    pd.DataFrame({
        "chr": chrom,
        "snpid": snp_ids,
        "cm": 0,
        "pos": positions.astype(np.int64),
        "a1": alt,
        "a2": ref,
    }).to_csv(f"{prefix}.bim", sep="\t", index=False, header=False)

    pd.DataFrame({
        "fid": 0, "iid": iids, "pid": 0, "mid": 0, "sex": 0, "pheno": -9
    }).to_csv(f"{prefix}.fam", sep="\t", index=False, header=False)


# =================================================================
# 1-2. Recapitate (cached; shared across all replicates)
# =================================================================
if not RTS_FILE.exists():
    ts = tskit.load(TREE_FILE)
    rts = pyslim.recapitate(ts, recombination_rate=RECOMBINATION_RATE,
                            ancestral_Ne=SIM_NE, random_seed=RNG_SEED)
    rts.dump(RTS_FILE)
    n_multiroot_after = sum(1 for t in rts.trees() if t.num_roots > 1)
    print(f"After recapitation: {n_multiroot_after} not yet coalesced "
          f"(should be 0)", flush=True)
    assert n_multiroot_after == 0, "Recapitation incomplete - check ancestral_Ne"
else:
    rts = tskit.load(RTS_FILE)


# =================================================================
# 3. Overlay neutral QTL mutations
#    pi-scaling: mutate every site at 1% of the rate rather than designating
#    1% of sites as targets. Same expected count, no target bookkeeping.
# =================================================================
print("Overlaying mutations...", flush=True)
mts = msprime.sim_mutations(
    rts, rate=MU * PI_TARGET,
    random_seed=seed_mutations.generate_state(1)[0],
    model=msprime.JC69(), keep=True,
)
print(f"After mutation overlay: {mts.num_sites} sites, "
      f"{mts.num_mutations} mutations", flush=True)


# =================================================================
# 4. Strip multiallelic sites from mts itself, so everything downstream
#    is biallelic by construction.
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

if SAVE_MUT_TREES:
    mts.dump(MUT_TREE_FILE)


# =================================================================
# 5. Single pass over variants.
#    Collects, in one traversal: age, position, ref/alt alleles, population
#    allele frequency, and the dosage matrix for the sampled individuals.
#    This is the expensive step; nothing else re-traverses the tree sequence.
# =================================================================
n_dip_all = mts.num_samples // 2
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
# 6. Effect sizes (neutral trait: independent of age and frequency)
# =================================================================
beta = rng_beta.normal(0, SIGMA_BETA, size=M_raw)


# =================================================================
# 7. Ground-truth V_M and V_A checks
# =================================================================
u_target = MU * PI_TARGET * mts.sequence_length
V_M_true = 2 * u_target * SIGMA_BETA**2
V_A_analytic = 2 * SIM_NE * V_M_true
V_A_empirical = np.sum(2 * freqs * (1 - freqs) * beta**2)

print(f"\nV_M (true, analytic):            {V_M_true:.4g}")
print(f"V_A (analytic 2*Ne*V_M):         {V_A_analytic:.4g}")
print(f"V_A (empirical, sum 2pq*beta^2): {V_A_empirical:.4g}")
print(f"Ratio empirical/analytic: {V_A_empirical / V_A_analytic:.3f}  "
      f"(should be close to 1 if at equilibrium)", flush=True)


# =================================================================
# 8. Genetic values and phenotypes for the sampled individuals.
#    Chunked matmul: dosages.T @ beta in one go would materialise a ~14 GB
#    float array, so accumulate over blocks of variants instead.
# =================================================================
print("\nBuilding genetic values...", flush=True)
g = np.zeros(n_sample, dtype=np.float64)
for start in range(0, M_raw, 2000):
    end = min(start + 2000, M_raw)
    g += dosages[start:end].T.astype(np.float32) @ beta[start:end]

V_E = np.var(g) * (1 - H2) / H2
y = g + rng_noise.normal(0, np.sqrt(V_E), n_sample)

print(f"V_A in sample: {np.var(g):.4g},  V_E: {V_E:.4g}")
print(f"Realised h2: {np.var(g) / np.var(y):.3f} (target {H2})", flush=True)


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