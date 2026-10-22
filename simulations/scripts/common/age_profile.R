## ---------------------------------------------------------------------------
## age_profile.R -- build a realistic singleton age profile from the msprime run.
## Source alongside sim_common.R / sim_diag.R.
##
## Replaces flat_profile(K, t_max), which spread singletons uniformly over
## 0-t_max generations. Real singletons are piled up at the young end: at
## n = 500k the geometric mean age is ~9 generations and only ~0.04 per person
## are older than 200. That difference drives everything downstream -- it is
## why t0 falls from ~200 to ~10 generations and the back-extrapolation
## (1 + lambda*t0) stops mattering.
##
## Input: sing_ages.npy, one entry per singleton, pooled over individuals.
## ---------------------------------------------------------------------------

AGES_FILE_DEFAULT <- file.path(
  "/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations",
  "data/v1.2/n500000_L1e+08/sing_ages.npy")

## ---------------------------------------------------------------------------
## load_age_profile()
##
##   n_ind       individuals the ages were ascertained in (500,000 for that file)
##   n_bins      equal-count bins over the young range
##   young_max   edge of the young range; 205 = 2T12 growth onset, matching the
##               production script. Everything older goes in one catch-all bin.
##   scale_total per-person singleton total to simulate. NULL keeps the observed
##               (~35 at L = 1e8); 1000 scales to genome-wide. The SHAPE of the
##               age distribution is a property of the coalescent, not of
##               sequence length, so it carries over unchanged -- only the count
##               is rescaled.
##   L_bp        sequence length the ages came from, used only for the 2*mu*L
##               consistency check.
##
## Returns prof: t, p, w, edges, lam_tot, u_per_gen, plus diagnostics.
## ---------------------------------------------------------------------------
load_age_profile <- function(file = AGES_FILE_DEFAULT,
                             n_ind = 5e5, n_bins = 7, young_max = 205,
                             scale_total = NULL, L_bp = 1e8, mu = 1.25e-8,
                             verbose = TRUE) {

  a <- RcppCNPy::npyLoad(file)
  a <- a[is.finite(a) & a >= 0]
  n_sing <- length(a)
  mean_T_obs <- n_sing / n_ind

  ## Equal-count bins over the young range, one catch-all above it. Equal counts
  ## mean roughly equal information per bin, which is what the estimator wants.
  young <- a[a < young_max]
  if (!length(young)) stop("no singletons younger than ", young_max)
  q <- unique(quantile(young, seq(0, 1, length.out = n_bins + 1), names = FALSE))
  q[1] <- 0; q[length(q)] <- young_max
  edges <- c(q, Inf)

  k   <- findInterval(a, edges, rightmost.closed = FALSE)
  cnt <- tabulate(k, nbins = length(edges) - 1L)
  keep <- cnt > 0
  ## Count-weighted mean age within each bin -- the right representative age for
  ## a bin, and the one merge_bins() and the estimator both assume.
  tbar <- vapply(seq_along(cnt), function(b)
    if (cnt[b] > 0) mean(a[k == b]) else NA_real_, numeric(1))
  tmed <- vapply(seq_along(cnt), function(b)
    if (cnt[b] > 0) median(a[k == b]) else NA_real_, numeric(1))

  lo <- head(edges, -1)[keep]; hi <- tail(edges, -1)[keep]
  w  <- ifelse(is.finite(hi), hi - lo, NA_real_)      # top bin has no width
  p  <- cnt[keep] / n_sing
  t  <- tbar[keep]

  ## Mutations per person per generation.
  ## Every mutation that arose last generation is still a singleton, so the
  ## density of singleton ages at zero IS the mutation rate: U = T * P(age < 1).
  ## This replaces the old SING_PER_BIN * K / t_max, which assumed uniform ages.
  lam_tot   <- if (is.null(scale_total)) mean_T_obs else scale_total
  u_per_gen <- lam_tot * mean(a < 1)
  u_expect  <- 2 * mu * L_bp * (lam_tot / mean_T_obs)   # what it should be

  prof <- list(t = t, p = p, w = w, edges = c(lo, tail(hi, 1)),
               lam_tot = lam_tot, u_per_gen = u_per_gen,
               t_med = tmed[keep], n_sing = cnt[keep],
               mean_T_obs = mean_T_obs, n_ind = n_ind,
               t0 = sum(p * t), gm_age = exp(mean(log(pmax(a, 1e-8)))))

  if (verbose) {
    cat(sprintf("%s\n  %d singletons in %d people -> %.1f each\n",
                basename(file), n_sing, n_ind, mean_T_obs))
    cat(sprintf("  geometric mean age %.1f gen; count-weighted mean (t0) %.1f gen\n",
                prof$gm_age, prof$t0))
    cat(sprintf("  simulating %.0f per person -> U = %.3f per generation\n",
                lam_tot, u_per_gen))
    cat(sprintf("  2*mu*L predicts %.3f  -> ratio %.2f %s\n", u_expect,
                u_per_gen / u_expect,
                if (abs(log(u_per_gen / u_expect)) < log(1.25)) "OK" else "<-- CHECK"))
    print(data.frame(bin = sprintf("%.1f-%s", lo, ifelse(is.finite(hi),
                                                          sprintf("%.1f", hi), "Inf")),
                     n = cnt[keep], frac = round(p, 4),
                     mean_age = round(t, 1), median_age = round(tmed[keep], 1),
                     per_person = round(p * lam_tot, 2), row.names = NULL))
  }
  prof
}

## ---------------------------------------------------------------------------
## The 2*mu*L check, separately, because it is the one that tests the young end
## of the profile -- where nearly all the singletons are and where the
## coalescent age check (2*Ne/n) has the least resolution.
## ---------------------------------------------------------------------------
check_u <- function(prof, L_bp = 1e8, mu = 1.25e-8) {
  u_expect <- 2 * mu * L_bp * (prof$lam_tot / prof$mean_T_obs)
  c(u_observed = prof$u_per_gen, u_expected = u_expect,
    ratio = prof$u_per_gen / u_expect)
}

## ---------------------------------------------------------------------------
## Drop-in for sim_dataset(): same model, realistic ages.
## Note what is NO LONGER a free parameter -- K and t_max both come from the
## profile now, so they leave the sweep grid.
## ---------------------------------------------------------------------------
sim_dataset_real <- function(lambda, REP, H2, N_inds, prof,
                             cv_total = SIM_DEFAULTS$cv_total,
                             cv_age   = SIM_DEFAULTS$cv_shape,
                             count_model = SIM_DEFAULTS$count_model,
                             curve = SIM_DEFAULTS$curve,
                             sigma_beta = SIM_DEFAULTS$sigma_beta) {
  set.seed(make_seed(lambda, REP, H2, N_inds, prof$lam_tot, cv_total, cv_age,
                     count_model, curve, length(prof$t)))
  C  <- sim_counts(N_inds, prof, prof$lam_tot, cv_total, cv_age, count_model)
  ph <- make_pheno(C, prof, lambda, H2, sigma_beta, curve)

  ## V_M = U * sigma2_m, with U from the age-zero density of the profile.
  V_P      <- var(ph$Y)
  h2m_true <- prof$u_per_gen * sigma_beta^2 / V_P
  t0_dat   <- sum(colMeans(C) * prof$t) / sum(colMeans(C))

  c(list(C = C, t_mid = prof$t, prof = prof,
         u_per_gen = prof$u_per_gen, V_P = V_P,
         h2m_true = h2m_true,
         h2m_t0_true = h2m_true * decay(lambda, t0_dat, curve),
         t0_dat = t0_dat), ph)
}
