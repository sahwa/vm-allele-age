## ---------------------------------------------------------------
## Allele age perturbation study
##
## Question: if allele ages are estimated with realistic tsdate error
## rather than known exactly, how much persistence time (T = V_A/V_M)
## do we lose?
##
## This version perturbs ages and recomputes per-bin variance directly
## from 2pq*beta^2, i.e. it assumes a PERFECT variance estimator and
## isolates the effect of dating error alone. That is an upper bound
## on how much information survives; re-running GENIE on the perturbed
## annotations is a separate, stricter test.
##
## Error model anchored on Pope et al. 2026 Table S1 (VG algorithm):
##   bias  -0.03 to -0.00, RMSE 0.27 to 0.33, in log10 generations
## Those are averages over the whole frequency spectrum. We assume
## error is worse for rare variants, and that rare/young variants are
## biased UPWARD (Fig S12: truncated haplotypes inflate recent ages)
## while common/old variants are biased slightly downward (root-age
## regularisation compresses the deep end).
## ---------------------------------------------------------------

library(data.table)
library(ggplot2)

DATA    <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.1/replicates"
N_REPS  <- 100
BINS    <- c(0, 1e2, 1e3, 1e4, 5e4, 1e5, 2e5, 5e5, Inf)
SIM_NE  <- 100000          # T_true should come out near 2*Ne = 200,000

bin_labels <- sapply(seq_len(length(BINS) - 1), function(i) {
    hi <- if (is.infinite(BINS[i + 1])) "inf" else format(BINS[i + 1], scientific = FALSE)
    paste0(format(BINS[i], scientific = FALSE), "-", hi)
})

## ---------------------------------------------------------------
## Error model
##
## u = 0 for the rarest variants, u = 1 for the commonest.
##   sigma(u) = sigma_base * (1 + k_sigma * (1 - u))    wider at low MAF
##   bias(u)  = bias_common + (bias_rare - bias_common) * (1 - u)
##
## sigma is then rescaled so the realised overall RMSE,
## sqrt(mean(bias^2 + sigma^2)), matches `target_rmse`. This makes the
## sweep interpretable: "target_rmse = 0.27" means the whole-spectrum
## RMSE is 0.27, matching the published number, with the error
## redistributed across the frequency spectrum.
## ---------------------------------------------------------------
error_params <- function(maf, target_rmse,
                         k_sigma = 0.8,
                         bias_rare = 0.15,
                         bias_common = -0.05) {

    lm_min <- log10(min(maf))
    u <- (log10(maf) - lm_min) / (log10(0.5) - lm_min)
    u <- pmin(pmax(u, 0), 1)

    sigma_rel <- 1 + k_sigma * (1 - u)
    bias      <- bias_common + (bias_rare - bias_common) * (1 - u)

    # scale sigma so realised RMSE hits the target
    scale <- sqrt(pmax(target_rmse^2 - mean(bias^2), 1e-8) / mean(sigma_rel^2))
    list(sigma = sigma_rel * scale, bias = bias)
}

## ---------------------------------------------------------------
## Persistence time fit (natural-scale objective, obj4)
## d needs columns bin_lo, bin_hi, share
## ---------------------------------------------------------------
fit_persistence <- function(d) {
    d <- d[share > 0]
    if (nrow(d) < 3) return(NA_real_)

    pred_share <- function(T) {
        K <- T * (exp(-d$bin_lo / T) -
                  ifelse(is.infinite(d$bin_hi), 0, exp(-d$bin_hi / T)))
        K / sum(K)
    }
    obj <- function(logT) sum((pred_share(exp(logT)) - d$share)^2)
    exp(optimize(obj, c(log(100), log(5e6)))$minimum)
}

## ---------------------------------------------------------------
## Per-bin variance from a set of (age, freq, beta) triples
## ---------------------------------------------------------------
bin_variance <- function(ages, freqs, beta) {
    b <- cut(ages, breaks = BINS, labels = bin_labels, right = FALSE)
    dt <- data.table(bin = b, contrib = 2 * freqs * (1 - freqs) * beta^2)
    agg <- dt[!is.na(bin), .(V = sum(contrib)), by = bin]

    out <- data.table(bin    = factor(bin_labels, levels = bin_labels),
                      bin_lo = BINS[-length(BINS)],
                      bin_hi = BINS[-1])
    out <- merge(out, agg, by = "bin", all.x = TRUE)
    out[is.na(V), V := 0]
    out[, share := V / sum(V)]
    out[]
}

## ---------------------------------------------------------------
## One replicate, one error setting
## ---------------------------------------------------------------
run_rep <- function(REP, target_rmse, seed_offset = 0, ...) {

    vi <- fread(file.path(DATA, sprintf("rep%d", REP), "1.1_variant_info.csv"))
    vi <- vi[age > 0 & freq > 0 & freq < 1]

    maf <- pmin(vi$freq, 1 - vi$freq)
    ep  <- error_params(maf, target_rmse, ...)

    set.seed(REP * 1000 + seed_offset)
    log_age_pert <- log10(vi$age) + ep$bias + rnorm(nrow(vi), 0, ep$sigma)
    age_pert     <- 10^log_age_pert

    truth <- bin_variance(vi$age,   vi$freq, vi$beta)
    smear <- bin_variance(age_pert, vi$freq, vi$beta)

    data.table(
        REP          = REP,
        target_rmse  = target_rmse,
        T_true       = fit_persistence(truth),
        T_smeared    = fit_persistence(smear),
        realised_rmse = sqrt(mean((log_age_pert - log10(vi$age))^2)),
        realised_bias = mean(log_age_pert - log10(vi$age))
    )
}

## ---------------------------------------------------------------
## Sweep
## ---------------------------------------------------------------
rmse_grid <- c(0.27, 0.30, 0.33)

res <- rbindlist(lapply(rmse_grid, function(r)
    rbindlist(lapply(1:N_REPS, run_rep, target_rmse = r))038956
    
))

res[, rel_err := (T_smeared - T_true) / T_true]

summary_bias <- res[, .(
    med_T_true    = median(T_true,    na.rm = TRUE),
    med_T_smeared = median(T_smeared, na.rm = TRUE),
    med_rel       = median(rel_err,   na.rm = TRUE),
    q10           = quantile(rel_err, 0.1, na.rm = TRUE),
    q90           = quantile(rel_err, 0.9, na.rm = TRUE),
    realised_rmse = mean(realised_rmse),
    realised_bias = mean(realised_bias),
    type = "normal_bias"
), by = target_rmse]

print(summary_tab)

## ---------------------------------------------------------------
## Sensitivity: how much damage is bias vs variance?
## Re-run at the middle RMSE with the bias term switched off.
## ---------------------------------------------------------------
res_nobias <- rbindlist(lapply(1:N_REPS, run_rep,
                               target_rmse = 0.30,
                               bias_rare = 0, bias_common = 0))
res_nobias[, rel_err := (T_smeared - T_true) / T_true]
summary_nobias = res_nobias[, .(
    med_T_true    = median(T_true,    na.rm = TRUE),
    med_T_smeared = median(T_smeared, na.rm = TRUE),
    med_rel       = median(rel_err,   na.rm = TRUE),
    q10           = quantile(rel_err, 0.1, na.rm = TRUE),
    q90           = quantile(rel_err, 0.9, na.rm = TRUE),
    realised_rmse = mean(realised_rmse),
    realised_bias = mean(realised_bias),
    type = "nobias"
), by = target_rmse]


res_halfbias <- rbindlist(lapply(1:N_REPS, run_rep, target_rmse = 0.30,
                                 bias_rare = 0.075, bias_common = -0.025))
res_halfbias[, rel_err := (T_smeared - T_true) / T_true]
summary_halfbias = res_halfbias[, .(
    med_T_true    = median(T_true,    na.rm = TRUE),
    med_T_smeared = median(T_smeared, na.rm = TRUE),
    med_rel       = median(rel_err,   na.rm = TRUE),
    q10           = quantile(rel_err, 0.1, na.rm = TRUE),
    q90           = quantile(rel_err, 0.9, na.rm = TRUE),
    realised_rmse = mean(realised_rmse),
    realised_bias = mean(realised_bias),
    type = "half_bias"
), by = target_rmse]

res_dblbias  <- rbindlist(lapply(1:N_REPS, run_rep, target_rmse = 0.30,
                                 bias_rare = 0.30, bias_common = -0.10))
res_dblbias[, rel_err := (T_smeared - T_true) / T_true]
summary_dblbias = res_dblbias[, .(
    med_T_true    = median(T_true,    na.rm = TRUE),
    med_T_smeared = median(T_smeared, na.rm = TRUE),
    med_rel       = median(rel_err,   na.rm = TRUE),
    q10           = quantile(rel_err, 0.1, na.rm = TRUE),
    q90           = quantile(rel_err, 0.9, na.rm = TRUE),
    realised_rmse = mean(realised_rmse),
    realised_bias = mean(realised_bias),
    type = "dbl_bias"
), by = target_rmse]

rbindlist(list(
    summary_bias,
    summary_nobias,
    summary_halfbias,
    summary_dblbias
))

## ---------------------------------------------------------------
## Plot
## ---------------------------------------------------------------
p <- ggplot(res, aes(x = factor(target_rmse), y = rel_err)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_jitter(width = 0.2, alpha = 0.2, size = 1) +
    stat_summary(fun = median,
                 fun.min = ~quantile(.x, 0.1),
                 fun.max = ~quantile(.x, 0.9),
                 geom = "pointrange", colour = "#1b9e77", linewidth = 1) +
    labs(x = "Target log10 age RMSE",
         y = "Relative error in persistence time") +
    theme_light(base_size = 11)

ggsave("persistence_vs_dating_error.png", p, width = 7, height = 5, dpi = 200)