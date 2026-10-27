## ---------------------------------------------------------------------------
## date_accuracy_dt.R -- what to do with variants.csv (40.5M rows).
##
## Three deliverables, in order of how much they matter to the project:
##
##   (1) ERROR MATRICES, by MAC stratum. P(inferred age bin | true age bin).
##       This is the real version of the synthetic ERR_GRID in age_error.R.
##       Batching changes the MAC of the ascertained variants, so there is a
##       different matrix for the current design (MAC=1) and the batched one
##       (MAC~10-30) -- run the estimator under each and the difference IS the
##       dating half of the batching argument.
##
##   (2) ACCURACY CURVE vs MAC. The headline: does scatter actually fall
##       between MAC 1 and MAC 30? If it is flat, batching buys nothing on
##       dating and only the h2_sing argument survives.
##
##   (3) THE t0 AND CV_age HIT. Dating error shifts measured t0 and inflates
##       the sampling term in CV_age. Both feed the power formula directly.
## ---------------------------------------------------------------------------

library(data.table)
setDTthreads(0)

IN   <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/tsdate/variant_age_table.csv"          # or .parquet via arrow::read_parquet
OUT  <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/tsdate/date_accuracy"
dir.create(OUT, showWarnings = FALSE)

dat <- fread(IN, verbose=T)

## Same log2 grid as batch_sweep.py, so the matrices drop straight into the
## sweep without rebinning.
AGE_EDGES <- c(0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, Inf)
## MAC strata on powers of 2: a batch-singleton sits at ~N/B copies, so these
## map onto batch sizes B = N/1, N/2, N/4, ...
MAC_EDGES <- c(1, 2, 4, 8, 16, 32, 64, 128, 256, Inf)

lab <- function(e) {
  lo <- head(e, -1); hi <- tail(e, -1)
  ifelse(is.infinite(hi), paste0(lo, "+"),
         ifelse(hi - lo == 1, as.character(lo), paste0(lo, "-", hi - 1)))
}
AGE_LAB <- lab(AGE_EDGES); MAC_LAB <- lab(MAC_EDGES)

dat <- dat[age_true > 0 & age_inferred > 0]
dat[, `:=`(
  lr       = log(age_inferred / age_true),
  bin_true = findInterval(age_true,     AGE_EDGES),
  bin_inf  = findInterval(age_inferred, AGE_EDGES),
  mac_str  = findInterval(mac,          MAC_EDGES)
)]

## ---------------------------------------------------------------------------
## (2) accuracy vs MAC -- exact MAC at the low end, strata above
## ---------------------------------------------------------------------------
acc_exact <- dat[mac <= 64, .(
  n            = .N,
  age_true_med = median(age_true),
  age_inf_med  = median(age_inferred),
  bias_fold    = exp(mean(lr)),            # 1 = unbiased, 2 = twice too old
  scatter_fold = exp(sd(lr)),              # multiplicative scatter
  cor_log      = cor(log(age_true), log(age_inferred)),
  within_2x    = mean(abs(lr) < log(2))
), by = mac][order(mac)]

acc_str <- dat[, .(
  n            = .N,
  mac_med      = as.double(median(mac)),
  bias_fold    = exp(mean(lr)),
  scatter_fold = exp(sd(lr)),
  cor_log      = cor(log(age_true), log(age_inferred)),
  within_2x    = mean(abs(lr) < log(2)),
  ## the age-error model's two knobs, measured: a multiplicative bias and a
  ## lognormal CV. These are what make_age_bias()/make_age_error_matrix() want.
  cv_logage    = sd(lr)
), by = mac_str][order(mac_str)][, mac_bin := MAC_LAB[mac_str]][]

fwrite(acc_exact, file.path(OUT, "accuracy_by_mac_exact.csv"))
fwrite(acc_str,   file.path(OUT, "accuracy_by_mac_stratum.csv"))
print(acc_str)

######################

## Is the bias age-dependent as well as MAC-dependent? This is the one that
## actually distorts the decay curve: a uniform fold-bias just rescales t0,
## but a bias that varies with age bends f(lambda, t) and biases lambda.
bias_by_age <- dat[, .(n = .N, bias_fold = exp(mean(lr)), scatter = exp(sd(lr))),
                   by = .(mac_str, bin_true)][order(mac_str, bin_true)]
bias_by_age[, `:=`(mac_bin = MAC_LAB[mac_str], age_bin = AGE_LAB[bin_true])]
fwrite(bias_by_age, file.path(OUT, "bias_by_age_and_mac.csv"))

## ---------------------------------------------------------------------------
## (1) error matrices -- the main deliverable
## ---------------------------------------------------------------------------
K <- length(AGE_LAB)
err_long <- dat[, .N, by = .(mac_str, bin_true, bin_inf)]
err_long[, p := N / sum(N), by = .(mac_str, bin_true)]
fwrite(err_long, file.path(OUT, "error_matrix_long.csv"))

## ROW-STOCHASTIC: E[k_true, j_inferred] = P(assigned to j | truly in k).
## So observed counts are C_obs = C_true %*% E. Check your misclassify()
## convention and transpose if it wants the other orientation.
err_mats <- lapply(sort(unique(dat$mac_str)), function(s) {
  E <- matrix(0, K, K, dimnames = list(true = AGE_LAB, inferred = AGE_LAB))
  e <- err_long[mac_str == s]
  E[cbind(e$bin_true, e$bin_inf)] <- e$p
  empty <- rowSums(E) == 0          # bins with no variants: pass through
  E[empty, ] <- diag(K)[empty, ]
  E
})
names(err_mats) <- MAC_LAB[sort(unique(dat$mac_str))]
saveRDS(err_mats, file.path(OUT, "error_matrices.rds"))

cat("\n-- P(inferred bin | true bin), MAC = 1 (current design) --\n")
print(round(err_mats[["1"]], 3))
if ("16-31" %in% names(err_mats)) {
  cat("\n-- MAC 16-31 (B ~ 30,000) --\n")
  print(round(err_mats[["16-31"]], 3))
}

## Diagonal mass = fraction of variants landing in the right bin. If this does
## not climb substantially from MAC=1 to MAC~16-31, batching does not help the
## dating and the case rests on h2_sing alone.
diag_mass <- data.table(
  mac_bin  = names(err_mats),
  diagonal = sapply(err_mats, function(E) mean(diag(E))),
  off_by_1 = sapply(err_mats, function(E)
    mean(diag(E) + c(diag(E)[-1], 0) + c(0, diag(E)[-K]))))
fwrite(diag_mass, file.path(OUT, "diagonal_mass.csv"))
cat("\n"); print(diag_mass)

## ---------------------------------------------------------------------------
## (3) what the error does to t0 and to CV_age
## ---------------------------------------------------------------------------
## t0 measured with inferred ages vs the truth. A fold-shift here propagates
## straight into the back-extrapolation sigma2_m = sigma2_t0 * (1 + lambda*t0).
t0_tab <- dat[, .(t0_true = mean(age_true), t0_inferred = mean(age_inferred),
                  t0_fold = mean(age_inferred) / mean(age_true)),
              by = mac_str][order(mac_str)][, mac_bin := MAC_LAB[mac_str]][]
fwrite(t0_tab, file.path(OUT, "t0_shift.csv"))
cat("\n"); print(t0_tab)

## Dating error adds an independent per-variant log-age error, so a person's
## mean log age picks up extra noise var(lr)/T. That inflates V_samp and
## attenuates CV_age_eff = V_struct / sqrt(V_struct + V_samp).
V_STRUCT <- 0.111; V_SAMP <- 0.061; T_REF <- 35
cv_hit <- acc_str[, {
  extra <- cv_logage^2 / T_REF
  .(mac_bin, cv_logage,
    cv_eff_nodating = V_STRUCT / sqrt(V_STRUCT + V_SAMP),
    cv_eff_withdating = V_STRUCT / sqrt(V_STRUCT + V_SAMP + extra),
    se_inflation = sqrt(V_STRUCT + V_SAMP + extra) / sqrt(V_STRUCT + V_SAMP))
}]
fwrite(cv_hit, file.path(OUT, "cv_age_hit.csv"))
cat("\n"); print(cv_hit)

cat("\nwritten to", normalizePath(OUT), "\n")
