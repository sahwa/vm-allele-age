"""
Full pipeline: SLiM tree sequence (stabilising selection, mutations already
tracked forward) -> extract true ages, frequencies, effect sizes -> phenotypes
-> plink fileset -> GENIE annotation matrix.

Unlike v1.1, mutations are NOT overlaid afterwards with msprime: selection
acts on them during the forward simulation, so SLiM has to track them, and
their effect sizes are read from mutation metadata rather than drawn here.

Changes from the previous version:
  - Path structure fixed: replicates/VS_{V_S}_NE_{N_E}/{REP}/
  - ONE pass over mts.variants() instead of two (ages/freqs/beta, then g_all
    separately). Also fixes the previous version silently mixing `ts` (raw)
    and `mts` (post multiallelic-filter) across sections.
  - .bed/.bim/.fam written directly with numpy. No VCF, no plink subprocess.
  - g computed from the sampled dosage matrix, not a full-population g_all
    that was then 50% discarded.
  - MAF filtering applied to ages/beta/freqs/positions simultaneously, so the
    annotation aligns to the .bim by construction.
"""

import argparse
from pathlib import Path
import itertools

import numpy as np
import pandas as pd
import tskit


# =================================================================
# 0. Parameters
# =================================================================
parser = argparse.ArgumentParser(description="SLiM selection simulation pipeline")
parser.add_argument("--rep", type=int, required=True, help="Replicate number")
parser.add_argument("--vs", type=int, required=True, help="Stabilising selection strength (V_S)")
parser.add_argument("--ne", type=int, required=True, help="Effective population size (N_e)")
parser.add_argument("--tree", required=True, help="SLiM tree file (mutations already tracked)")
args = parser.parse_args()

REP = args.rep
V_S = args.vs
N_E = args.ne
TREE_FILE = args.tree
PARAMS = f"VS_{V_S}_NE_{N_E}"

RNG_SEED = 42
ss = np.random.SeedSequence([RNG_SEED, REP])
seed_sample, seed_noise = ss.spawn(2)   # no seed for beta/mutations: SLiM already drew them

rng_sample = np.random.default_rng(seed_sample)
rng_noise = np.random.default_rng(seed_noise)

MU = 1.44e-8                # per base per generation mutation rate
PI_TARGET = 0.01            # fraction of the genome that is mutational target

SIM_PATH = Path(
    "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v2.0"
)
SIM_VERSION = "2.0"
SELECTION_TYPE = "stabilising_selection"
FILE_STEM = f"{SIM_VERSION}_{SELECTION_TYPE}_{PARAMS}"

# /well/.../data/v2.0/replicates/VS_{V_S}_NE_{N_E}/{REP}/
SIM_PATH_REP = SIM_PATH / "replicates" / PARAMS / str(REP)
SIM_PATH_REP.mkdir(parents=True, exist_ok=True)

N_SAMPLE_TARGET = 50000     # diploids to sample for the GENIE analysis
H2 = 0.5

BINS = [0, 100, 1000, 10000, 50000, np.inf]

GENOME_PREFIX = SIM_PATH_REP / f"{FILE_STEM}.out"


# =================================================================
# Helper: write a PLINK1 .bed/.bim/.fam fileset directly
# =================================================================
def write_plink_bed(dosages, prefix, positions, ref, alt, iids, chrom=1):
    """
    dosages : (n_var, n_ind) uint8, count of the ALT (derived) allele, 0/1/2
    A1 = ALT, A2 = REF, matching plink2's --vcf default orientation.
    """
    
    n_var, n_ind = dosages.shape
    assert dosages.dtype == np.uint8, f"Expected uint8, got {dosages.dtype}"
    assert np.all((dosages >= 0) & (dosages <= 2)), \
        f"Dosages out of range: {np.unique(dosages)}"
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

    snp_ids = [f"{chrom}:{int(p)}:{r}:{a}" for p, r, a in zip(positions, ref, alt)]
    pd.DataFrame({
        "chr": chrom, "snpid": snp_ids, "cm": 0,
        "pos": positions.astype(np.int64), "a1": alt, "a2": ref,
    }).to_csv(f"{prefix}.bim", sep="\t", index=False, header=False)

    pd.DataFrame({
        "fid": 0, "iid": iids, "pid": 0, "mid": 0, "sex": 0, "pheno": -9
    }).to_csv(f"{prefix}.fam", sep="\t", index=False, header=False)


# =================================================================
# 1. Load the SLiM tree sequence (mutations already tracked forward)
# =================================================================
ts = tskit.load(TREE_FILE)

n_multiroot = sum(1 for t in ts.trees() if t.num_roots > 1)
print(f"Loaded: {ts.num_samples} samples, {ts.num_trees} trees, "
      f"{n_multiroot} not yet coalesced", flush=True)

frac_multiroot = n_multiroot / ts.num_trees
assert frac_multiroot < 0.01, (
    f"{100 * frac_multiroot:.1f}% of trees uncoalesced - burn-in too short "
    f"for Ne={N_E} (tree: {TREE_FILE})"
)

total_span, weighted_tmrca = 0, 0
for tree in ts.trees():
    root_time = tree.time(tree.root) if tree.num_roots == 1 else None
    if root_time is not None:
        weighted_tmrca += root_time * tree.span
        total_span += tree.span

mean_tmrca = weighted_tmrca / total_span
print(f"Mean TMRCA: {mean_tmrca:.0f} (expect ~4*Ne = {4 * N_E} under neutrality)")


# =================================================================
# 2. Empirical Ne from pairwise coalescence times
# =================================================================
rng_ne = np.random.default_rng(42)
sample_nodes = ts.samples()   # actual sample node IDs, not a raw range

idx_pairs = rng_ne.choice(len(sample_nodes), size=(10000, 2), replace=True)
sample_pairs = [(int(sample_nodes[a]), int(sample_nodes[b]))
                 for a, b in idx_pairs if a != b]

print(f"Using {len(sample_pairs)} pairs from {len(sample_nodes)} true sample nodes")
print(f"Sample node ID range: {sample_nodes.min()} - {sample_nodes.max()}")

tmrcas, n_failed = [], 0
for tree in itertools.islice(ts.trees(), 0, None, 1000):
    for a, b in sample_pairs:
        try:
            tmrcas.append(tree.tmrca(a, b))
        except ValueError:
            n_failed += 1

tmrcas = np.array(tmrcas)
print(f"Collected {len(tmrcas)}, failed {n_failed}")

NE_EMPIRICAL = tmrcas.mean() / 2
print(f"Mean pairwise TMRCA: {tmrcas.mean():.0f}")
print(f"Ne estimate: {NE_EMPIRICAL:.0f}")


# =================================================================
# 3. Strip multiallelic sites, so everything downstream is biallelic
#    by construction. Use `mts` consistently from here on.
# =================================================================
multiallelic_site_ids = np.array(
    [site.id for site in ts.sites() if len(site.mutations) > 1]
)
print(f"Removing {len(multiallelic_site_ids)} multiallelic sites "
      f"out of {ts.num_sites}", flush=True)

if len(multiallelic_site_ids) > 0:
    tables = ts.dump_tables()
    tables.delete_sites(multiallelic_site_ids)
    tables.sort()
    mts = tables.tree_sequence()
    print(f"mts now has {mts.num_sites} biallelic sites", flush=True)
else:
    print("No multi-allelics to remove", flush=True)
    mts = ts


# =================================================================
# 4. Single pass over variants.
#    Collects age, position, ref/alt, population frequency, effect size
#    (from SLiM mutation metadata), and the dosage matrix for the sampled
#    individuals, all in one traversal.
# =================================================================
n_dip_all = mts.num_samples // 2
n_sample = min(N_SAMPLE_TARGET, n_dip_all)
sample_idx = np.sort(rng_sample.choice(n_dip_all, n_sample, replace=False))
hap_idx = np.stack([2 * sample_idx, 2 * sample_idx + 1], axis=1).ravel()

M = mts.num_sites
dosages = np.empty((M, n_sample), dtype=np.uint8)
ages = np.empty(M)
freqs = np.empty(M)
positions = np.empty(M)
beta = np.empty(M)
ref = np.empty(M, dtype="<U1")
alt = np.empty(M, dtype="<U1")

print(f"\nSingle pass over {M} variants "
      f"({n_sample} of {n_dip_all} diploids sampled)...", flush=True)

for i, var in enumerate(mts.variants()):
    gt = var.genotypes
    mut = mts.mutation(var.site.mutations[0].id)
    md = mut.metadata["mutation_list"][0]

    beta[i] = md["selection_coeff"]     # SLiM stores the drawn effect here
    ages[i] = mut.time if not np.isnan(mut.time) else mts.node(mut.node).time
    freqs[i] = gt.mean()
    positions[i] = var.site.position
    ref[i] = var.alleles[0][:1] or "A"
    alt[i] = (var.alleles[1][:1] if len(var.alleles) > 1 and var.alleles[1] else "T")
    dosages[i] = gt[hap_idx].reshape(-1, 2).sum(axis=1)

pheno_proxy = freqs * beta   # legacy quantity kept for continuity with earlier notebooks

print(f"Age range: {ages.min():.1f} - {ages.max():.1f} generations")
print(f"Freq range: {freqs.min():.5f} - {freqs.max():.5f}")
print(f"Beta range: {beta.min():.5f} - {beta.max():.5f}")


# =================================================================
# 5. Ground-truth V_M, V_A, persistence time
# =================================================================
SIGMA_BETA = 0.1   # SLiM parameter, not re-estimated from data: the empirical
                   # SD of surviving betas is lower because selection has
                   # already removed the most deleterious alleles - that's
                   # the whole point, so it must not be used here.

print(f"SD of extracted betas: {beta.std():.4f} (SLiM param: {SIGMA_BETA})")

u_target = MU * PI_TARGET * mts.sequence_length
V_M_analytic_true = 2 * u_target * SIGMA_BETA**2

V_A_empirical = np.sum(2 * freqs * (1 - freqs) * beta**2)
V_P = np.var(pheno_proxy)

PERSISTENCE_TIME = V_A_empirical / V_M_analytic_true

# s(beta) = beta^2 / (2*(V_S + V_P)): V_P dilutes selection via phenotypic
# background, not V_S alone.
S_BAR = SIGMA_BETA**2 / (2 * (V_S + V_A_empirical))
print(f"V_A empirical (= V_P proxy): {V_A_empirical:.4g}")
print(f"S_BAR corrected: {S_BAR:.6f}")
print(f"Persistence 1/S_BAR: {1 / S_BAR:.0f} generations")

REGIME_PARAMETER = S_BAR * NE_EMPIRICAL
V_A_analytic_neutral = 2 * NE_EMPIRICAL * V_M_analytic_true
V_A_analytic_selection = V_M_analytic_true / S_BAR

print(f"N_e x s_bar = {REGIME_PARAMETER:.2f}")
if REGIME_PARAMETER > 1:
    print("  -> Selection dominates; expect V_A ~ V_M / s_bar")
elif REGIME_PARAMETER < 1:
    print("  -> Drift dominates; expect V_A ~ 2*N_e*V_M")
else:
    print("  -> Intermediate regime; expect neither prediction to be exact")

print(f"V_A empirical / V_A_MSB: {V_A_empirical / V_A_analytic_selection:.3f}  (should be ~ 1.0)")
print(f"V_A empirical / V_A neutral: {V_A_empirical / V_A_analytic_neutral:.3f}  (should be << 1.0 if selection dominates)")
print(f"Persistence time observed: {PERSISTENCE_TIME:.0f} generations")
print(f"Persistence time 1/s_bar (MSB): {1 / S_BAR:.0f} generations")
print(f"Persistence time 2Ne (neutral): {2 * NE_EMPIRICAL:.0f} generations")


# =================================================================
# 6. Genetic values and phenotypes for the sampled individuals.
#    Chunked matmul to avoid materialising a large float intermediate.
# =================================================================
print("\nBuilding genetic values...", flush=True)
g = np.zeros(n_sample, dtype=np.float64)
for start in range(0, M, 2000):
    end = min(start + 2000, M)
    g += dosages[start:end].T.astype(np.float32) @ beta[start:end]

V_E = np.var(g) * (1 - H2) / H2
y = g + rng_noise.normal(0, np.sqrt(V_E), n_sample)

print(f"Sampled {n_sample} of {n_dip_all} diploids")
print(f"V_A in sample: {np.var(g):.4g}, V_E: {V_E:.4g}")
print(f"Realised h2: {np.var(g) / np.var(y):.3f} (target {H2})", flush=True)


# =================================================================
# 7. True per-bin variance, over the full population (pre-MAF-filter)
# =================================================================
bin_labels = []
for i in range(len(BINS) - 1):
    lo = int(BINS[i])
    hi = "inf" if np.isinf(BINS[i + 1]) else int(BINS[i + 1])
    bin_labels.append(f"{lo}-{hi}")

bins_assigned = pd.cut(ages, bins=BINS, labels=bin_labels, right=False)

print(f"\n{'Bin (gens)':<18}{'n':>7}{'observed':>12}")
bin_rows = []
for lo, hi in zip(BINS[:-1], BINS[1:]):
    m = (ages >= lo) & (ages < hi)
    if m.sum() == 0:
        continue
    contribs = 2 * freqs[m] * (1 - freqs[m]) * beta[m]**2
    V_bin = contribs.sum()

    label = f"{lo:.0f}+" if np.isinf(hi) else f"{lo:.0f}-{hi:.0f}"
    print(f"{label:<18}{m.sum():>7}{V_bin:>12.4g}")

    bin_rows.append({"bin_lo": lo, "bin_hi": hi, "n_variants": int(m.sum()), "V_observed": V_bin})

pd.DataFrame(bin_rows).to_csv(
    SIM_PATH_REP / f"{FILE_STEM}_bin_truth.csv", index=False)


# =================================================================
# 8. MAF filter in the sample, applied to every per-variant array at once,
#    so the annotation aligns to the .bim by construction.
# =================================================================
sample_freq = dosages.mean(axis=1) / 2
keep = (sample_freq > 0) & (sample_freq < 1)

print(f"\nAfter MAF filtering: {keep.sum()} of {M} variants retained "
      f"({M - keep.sum()} dropped as monomorphic in the sample)", flush=True)

dosages_f = dosages[keep]
positions_f = positions[keep]
ages_f = ages[keep]
ref_f, alt_f = ref[keep], alt[keep]

indv_names = [f"ind{i}" for i in range(n_sample)]

# Validate dosages before packing
assert np.all((dosages_f >= 0) & (dosages_f <= 2)), \
    f"Invalid dosage values: min={dosages_f.min()}, max={dosages_f.max()}, NaN={np.isnan(dosages_f).sum()}"
assert not np.isnan(dosages_f).any(), "NaN in dosages"

print(f"Dosages OK: shape {dosages_f.shape}, range [{dosages_f.min()}, {dosages_f.max()}]", flush=True)

write_plink_bed(dosages_f, GENOME_PREFIX, positions_f, ref_f, alt_f, indv_names)


# =================================================================
# 9. GENIE annotation matrix. Row i corresponds to row i of the .bim
#    by construction, since both come from the same `keep` mask.
# =================================================================
bins_assigned_filtered = pd.cut(ages_f, bins=BINS, labels=bin_labels, right=False)

annotations = pd.get_dummies(
    bins_assigned_filtered, prefix="bin", prefix_sep="_"
).astype(int)

assert len(annotations) == keep.sum(), "annotation rows != retained variants"
assert (annotations.sum(axis=1) == 1).all(), "each variant must fall in exactly one bin"

annotations.to_csv(
    SIM_PATH_REP / f"{FILE_STEM}_annotations_age_bins.txt",
    sep=" ", index=False, header=False,
)

pd.DataFrame({
    "column_name": annotations.columns,
    "age_bin": bin_labels,
}).to_csv(
    SIM_PATH_REP / f"{FILE_STEM}_annotations_legend.txt",
    sep=" ", index=False,
)


# =================================================================
# 10. Variant info and phenotype files
# =================================================================
pd.DataFrame({
    "site_id": np.arange(M),
    "position": positions,
    "age": ages,
    "bin": bins_assigned,
    "freq": freqs,
    "beta": beta,
    "kept": keep,
}).to_csv(SIM_PATH_REP / f"{FILE_STEM}_variant_info.csv", index=False)

pd.DataFrame({"FID": 0, "IID": indv_names, "y": y}).to_csv(
    SIM_PATH_REP / f"{FILE_STEM}_phenotypes.GENIE.txt",
    sep="\t", index=False, header=["FID", "IID", "PHENO"]
)

pd.DataFrame({
    "sample_idx": sample_idx,
    "iid": indv_names,
    "y": y,
    "g": g,
}).to_csv(SIM_PATH_REP / f"{FILE_STEM}_phenotypes.csv", index=False)

print(f"\nWrote outputs to {SIM_PATH_REP}", flush=True)