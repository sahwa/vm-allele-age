#!/usr/bin/env python3
"""
Decisive diagnostic for near-zero GENIE heritability.

Recomputes genetic values directly from the PLINK1 .bed fileset and compares
them against the g / y that the pipeline wrote out. Pure numpy - no
pandas_plink, no plink2 call.

Answers, in order:
  A. Do the annotation / phenotype files line up with the .bim / .fam at all?
     (a mismatch here zeroes h2 even when the genotypes are perfect)
  B. Does the .bed encode the same dosages the pipeline thinks it does?
  C. Does g recomputed from the .bed match the g that produced y?

Run:  python check_alignment.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- config
SIM_PATH = Path(
    "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.1"
)
SIM_VERSION = "1.1"

PREFIX = SIM_PATH / f"{SIM_VERSION}_neutral_out"
VARIANT_INFO = SIM_PATH / f"{SIM_VERSION}_variant_info.csv"
PHENO_CSV = SIM_PATH / f"{SIM_VERSION}_phenotypes.csv"
PHENO_TXT = SIM_PATH / f"{SIM_VERSION}_phenotypes.txt"
ANNOT = SIM_PATH / f"{SIM_VERSION}_annotations_age_bins.txt"

CHUNK = 2000          # variants decoded at a time
H2_TARGET = 0.5

ok = lambda b: "PASS" if b else "**FAIL**"


# ------------------------------------------------- plink1 .bed decoding
# 2 bits per genotype, SNP-major. Codes are:
#   00 -> 2 copies of A1 (.bim col 5)
#   01 -> missing
#   10 -> 1 copy of A1
#   11 -> 0 copies of A1
_CODE = {0: 2.0, 1: np.nan, 2: 1.0, 3: 0.0}
_LUT = np.empty((256, 4), dtype=np.float32)
for _b in range(256):
    for _j in range(4):
        _LUT[_b, _j] = _CODE[(_b >> (2 * _j)) & 3]


def iter_bed(bed_path, n_ind, n_var, chunk=CHUNK):
    """Yield (start_index, dosage_block) with dosage_block shape (n_var_chunk, n_ind)."""
    bytes_per_var = (n_ind + 3) // 4
    with open(bed_path, "rb") as f:
        magic = f.read(3)
        if magic[:2] != b"\x6c\x1b":
            raise ValueError(f"{bed_path} is not a PLINK1 .bed file")
        if magic[2] != 1:
            raise ValueError("individual-major .bed; re-run --make-bed to get SNP-major")

        expected = bytes_per_var * n_var
        actual = bed_path.stat().st_size - 3
        if actual != expected:
            raise ValueError(
                f".bed size mismatch: {actual} bytes for {n_var} variants x "
                f"{n_ind} individuals, expected {expected}. .bim/.fam do not "
                f"describe this .bed."
            )

        start = 0
        while start < n_var:
            nv = min(chunk, n_var - start)
            raw = np.fromfile(f, dtype=np.uint8, count=nv * bytes_per_var)
            raw = raw.reshape(nv, bytes_per_var)
            dos = _LUT[raw].reshape(nv, bytes_per_var * 4)[:, :n_ind]
            yield start, dos
            start += nv


# ---------------------------------------------------------------- load
bim = pd.read_csv(f"{PREFIX}.bim", sep=r"\s+", header=None,
                  names=["chr", "snpid", "cm", "pos", "a1", "a2"])
fam = pd.read_csv(f"{PREFIX}.fam", sep=r"\s+", header=None,
                  names=["fid", "iid", "pat", "mat", "sex", "phe"])
n_var, n_ind = len(bim), len(fam)
print(f"bim: {n_var} variants   fam: {n_ind} individuals\n")


# ============================================ A. file-shape consistency
print("=" * 62)
print("A. FILE SHAPES (a mismatch zeroes h2 with genotypes intact)")
print("=" * 62)

n_annot = sum(1 for _ in open(ANNOT))
print(f"annotation rows {n_annot} vs bim rows {n_var}   {ok(n_annot == n_var)}")
if n_annot != n_var:
    print("   -> GENIE reads annotations positionally. Every SNP after the first")
    print("      duplicated position is in the wrong bin. This alone gives ~0 h2.")

n_pheno = sum(1 for _ in open(PHENO_TXT))
print(f"phenotype rows {n_pheno} vs fam rows {n_ind}   {ok(n_pheno == n_ind)}")

# The duplicate-position trap in step 11 of the pipeline.
vinfo = pd.read_csv(VARIANT_INFO)
dup_full = int(vinfo["position"].duplicated().sum())
dup_bim = int(bim["pos"].duplicated().sum())
print(f"duplicated positions: variant_info {dup_full}, bim {dup_bim}   "
      f"{ok(dup_full == 0 and dup_bim == 0)}")
if dup_full or dup_bim:
    print("   -> pos_to_age.loc[...] returns extra rows for duplicated labels,")
    print("      inflating ages_filtered and the annotation matrix.")

annot_head = pd.read_csv(ANNOT, sep=" ", header=None, nrows=200000)
col_sums = annot_head.sum(axis=0)
empty_cols = list(np.where(col_sums.values == 0)[0])
row_sums = annot_head.sum(axis=1)
print(f"all-zero annotation columns: {empty_cols or 'none'}   {ok(not empty_cols)}")
print(f"rows assigned to no bin: {int((row_sums == 0).sum())}   "
      f"{ok((row_sums == 0).sum() == 0)}")
print(f"rows in >1 bin: {int((row_sums > 1).sum())}   {ok((row_sums > 1).sum() == 0)}")
print()


# ================================== align variant_info onto the bim rows
# Positional, not label-based, so duplicates cannot silently change length.
vinfo_pos = vinfo["position"].to_numpy()
order = np.argsort(vinfo_pos, kind="stable")
idx = order[np.searchsorted(vinfo_pos[order], bim["pos"].to_numpy())]
matched = vinfo_pos[idx] == bim["pos"].to_numpy()
print(f"bim positions found in variant_info: {matched.sum()}/{n_var}   "
      f"{ok(matched.all())}")
if not matched.all():
    raise SystemExit("Cannot align .bim to variant_info by position - stop here.")

beta = vinfo["beta"].to_numpy()[idx]
freq_expected = vinfo["freq"].to_numpy()[idx]
bin_of_var = vinfo["bin"].to_numpy()[idx]
print()


# =============================================== B/C. stream the .bed
print("=" * 62)
print("B. DECODING GENOTYPES")
print("=" * 62)

g_check = np.zeros(n_ind, dtype=np.float64)
freq_obs = np.empty(n_var, dtype=np.float64)
n_missing = 0

bins_present = pd.unique(bin_of_var[pd.notna(bin_of_var)])
g_by_bin = {b: np.zeros(n_ind, dtype=np.float64) for b in bins_present}

for start, dos in iter_bed(Path(f"{PREFIX}.bed"), n_ind, n_var):
    sl = slice(start, start + dos.shape[0])
    nan_mask = np.isnan(dos)
    if nan_mask.any():
        n_missing += int(nan_mask.sum())
        dos = np.where(nan_mask, 0.0, dos)

    freq_obs[sl] = dos.mean(axis=1) / 2.0
    b = beta[sl]
    g_check += b @ dos

    labels = bin_of_var[sl]
    for lab in np.unique(labels[pd.notna(labels)]):
        sel = labels == lab
        g_by_bin[lab] += b[sel] @ dos[sel]

    if (start // CHUNK) % 25 == 0:
        print(f"  {start + dos.shape[0]}/{n_var} variants", flush=True)

print(f"missing genotype calls: {n_missing}   {ok(n_missing == 0)}")

r_freq = np.corrcoef(freq_obs, freq_expected)[0, 1]
r_flip = np.corrcoef(freq_obs, 1 - freq_expected)[0, 1]
print(f"corr(freq_bed, freq_expected) = {r_freq:+.4f}")
print(f"corr(freq_bed, 1 - freq_expected) = {r_flip:+.4f}")
if r_flip > r_freq:
    print("   **A1 is the ancestral allele, not the derived one.**")
    print("   beta is defined on the derived allele -> dosages are 2-x.")
    print("   Re-run --make-bed with explicit REF/ALT handling, or negate beta.")
else:
    print(f"   allele orientation {ok(r_freq > 0.99)} "
          f"(sample vs population sampling noise is expected)")
print()


print("=" * 62)
print("C. GENETIC VALUES")
print("=" * 62)

pheno = pd.read_csv(PHENO_CSV)
# Join on IID rather than trusting row order.
pheno = pheno.set_index("iid").reindex(fam["iid"].to_numpy())
if pheno["y"].isna().any():
    raise SystemExit(
        "fam IIDs do not match phenotypes.csv iids - the VCF sample names and "
        "indv_names have diverged."
    )
y = pheno["y"].to_numpy()
g_saved = pheno["g"].to_numpy()

r_gg = np.corrcoef(g_check, g_saved)[0, 1]
r_gy = np.corrcoef(g_check, y)[0, 1]
r_sy = np.corrcoef(g_saved, y)[0, 1]

print(f"corr(g_from_bed, g_saved) = {r_gg:+.4f}   (expect ~1.00)   {ok(r_gg > 0.99)}")
print(f"corr(g_from_bed, y)       = {r_gy:+.4f}   (expect ~{np.sqrt(H2_TARGET):.2f})")
print(f"corr(g_saved,    y)       = {r_sy:+.4f}   (expect ~{np.sqrt(H2_TARGET):.2f})")
print()
print(f"var(g_from_bed) = {np.var(g_check):.4g}")
print(f"var(g_saved)    = {np.var(g_saved):.4g}")
print(f"implied h2 from bed genotypes = {r_gy**2:.4f}")
print()

print("per-bin corr(g_bin, y):")
for lab in sorted(g_by_bin, key=str):
    gb = g_by_bin[lab]
    if np.var(gb) == 0:
        print(f"  {str(lab):<14} var 0 - no contributing variants")
        continue
    print(f"  {str(lab):<14} r = {np.corrcoef(gb, y)[0, 1]:+.4f}   "
          f"var = {np.var(gb):.4g}")
print()


# ---------------------------------------------------------------- verdict
print("=" * 62)
print("VERDICT")
print("=" * 62)
if r_gg > 0.99 and r_gy > 0.5 and n_annot == n_var:
    print("Genotypes, phenotypes and annotations are all consistent.")
    print("The bug is in the GENIE invocation, not these files. Check that")
    print("--annot row count matches --bed, that the phenotype file's header")
    print("convention matches what your GENIE build expects, and that the")
    print("number of random vectors / jackknife blocks is sane for M and N.")
elif r_gg > 0.99 and n_annot != n_var:
    print("Genotypes and phenotypes are fine; the ANNOTATION FILE is misaligned.")
    print("Rebuild it positionally from the .bim and re-run.")
elif r_gg <= 0.99:
    print("The .bed does not reproduce g. Scrambling is upstream of GENIE,")
    print("in the VCF write or the MAF-filter reindex. Compare freq_bed against")
    print("freq_expected per variant to localise it.")
else:
    print("Mixed signals - read sections A-C above individually.")