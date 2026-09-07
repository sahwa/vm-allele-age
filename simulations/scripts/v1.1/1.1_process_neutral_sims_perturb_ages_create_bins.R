## ---------------------------------------------------------------
## Build GENIE annotation matrices under perturbed allele ages.
##
## Genotypes and .bim are unchanged - dating error doesn't touch
## genotypes. Only the annotation (which bin each SNP falls in) needs
## rebuilding per error setting. This script writes one annotation
## file per (REP, target_rmse) pair, aligned row-for-row to the
## existing MAF-filtered .bim, ready for GENIE.
##
## Run the perturbation+error model exactly as in
## perturb_ages_persistence.R so results are comparable to the
## "perfect estimator" experiment already done.
## ---------------------------------------------------------------

library(data.table)

DATA   <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.1/replicates"
N_REPS <- 100
BINS   <- c(0, 1e2, 1e3, 1e4, 5e4, 1e5, 2e5, 5e5, Inf)

bin_labels <- sapply(seq_len(length(BINS) - 1), function(i) {
    hi <- if (is.infinite(BINS[i + 1])) "inf" else format(BINS[i + 1], scientific = FALSE)
    paste0(format(BINS[i], scientific = FALSE), "-", hi)
})

## Same error model as perturb_ages_persistence.R - keep in sync if
## that script changes.
error_params <- function(maf, target_rmse,
                         k_sigma = 0.8,
                         bias_rare = 0.15,
                         bias_common = -0.05) {
    lm_min <- log10(min(maf))
    u <- (log10(maf) - lm_min) / (log10(0.5) - lm_min)
    u <- pmin(pmax(u, 0), 1)

    sigma_rel <- 1 + k_sigma * (1 - u)
    bias      <- bias_common + (bias_rare - bias_common) * (1 - u)

    scale <- sqrt(pmax(target_rmse^2 - mean(bias^2), 1e-8) / mean(sigma_rel^2))
    list(sigma = sigma_rel * scale, bias = bias)
}

## ---------------------------------------------------------------
## For one replicate/error setting: perturb ages, re-bin, write an
## annotation file aligned to the existing MAF-filtered .bim by
## position. This mirrors how the original pipeline built the
## annotation from variant_info.csv, so alignment logic matches.
## ---------------------------------------------------------------
build_annotation <- function(REP, target_rmse, seed_offset = 0) {

    rep_dir <- file.path(DATA, sprintf("rep%d", REP))
    vi  <- fread(file.path(rep_dir, "1.1_variant_info.csv"))
    bim <- fread(file.path(rep_dir, "1.1_neutral_out.bim"), header = FALSE,
                col.names = c("chr", "snpid", "cm", "pos", "a1", "a2"))

    # only variants that survived MAF filtering are in the .bim
    vi <- vi[vi$position %in% bim$pos]
    vi <- vi[match(bim$pos, vi$position)]   # force .bim row order
    # stopifnot(identical(vi$position, bim$pos))

    maf <- pmin(vi$freq, 1 - vi$freq)
    ep  <- error_params(maf, target_rmse)

    set.seed(REP * 1000 + seed_offset)
    log_age_pert <- log10(vi$age) + ep$bias + rnorm(nrow(vi), 0, ep$sigma)
    age_pert     <- 10^log_age_pert

    b <- cut(age_pert, breaks = BINS, labels = bin_labels, right = FALSE)

    ann <- as.data.table(model.matrix(~ b - 1))
    setnames(ann, paste0("bin_", bin_labels))

    stopifnot(nrow(ann) == nrow(bim))

    out_dir <- file.path(rep_dir, "dating_error")
    dir.create(out_dir, showWarnings = FALSE)
    out_file <- file.path(out_dir,
                          sprintf("1.1_annotations_age_bins.rmse%.2f.txt", target_rmse))
    fwrite(ann, out_file, sep = " ", col.names = FALSE)

    out_file
}

## ---------------------------------------------------------------
## Build annotations for the whole sweep. GENIE calls are left to a
## SLURM array (below) rather than run here, since 100 reps x 3
## rmse values x GENIE is worth parallelising properly.
## ---------------------------------------------------------------
rmse_grid <- c(0.27, 0.30, 0.33)

for (r in rmse_grid) {
    for (REP in 1:N_REPS) {
        build_annotation(REP, r)
    }
    cat("Built annotations for target_rmse =", r, "\n")
}