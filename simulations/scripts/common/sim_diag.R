## ---------------------------------------------------------------------------
## sim_common.R -- shared simulation functions for the singleton V_M project.
## Source this at the top of every run so the generative model is identical
## everywhere.
##
## TWO KNOBS, and they feed different parameters:
##
##   cv_total  how much people vary in HOW MANY singletons they carry.
##             This is the signal for sigma2_m (V_M). Measured: ~0.38 in the
##             n=500k msprime run; 0.37-0.65 in UKB.
##
##   cv_shape  how much people vary in WHICH AGES those singletons are, as the
##             SD across people of their mean log singleton age.
##             This is the signal for s. Measured: 0.33 in the same msprime run
##             (total var 0.172, sampling var 0.061, structural 0.111).
##             cv_shape = 0 reproduces the fully "shared" model.
##
## The old independent-gamma-per-bin model is kept as count_model = "indep"
## for back-compatibility, but it is not biologically defensible: it gives
## each person a randomly scrambled age profile.
## ---------------------------------------------------------------------------

SIM_DEFAULTS <- list(
  sigma_beta = 1,
  cv_total   = 0.65,     # between-person SD/mean of total singleton count
  cv_shape   = 0.33,     # between-person SD of mean log singleton age
  curve      = "hyperbolic",
  count_model = "shared" # "shared" (+ tilt) or "indep" (old, scrambled)
)

## --- decay of per-mutation variance with allele age -------------------------
decay <- function(s, t, curve = SIM_DEFAULTS$curve) {
  switch(curve,
    exponential = exp(-s * t),
    hyperbolic  = 1 / (1 + s * t),
    stop("curve must be 'exponential' or 'hyperbolic'"))
}

## --- deterministic seeds ----------------------------------------------------
## Hash the parameters so a cell is reproducible without depending on the
## position of S in some global vector.
make_seed <- function(...) {
  key <- paste(vapply(list(...), function(x) paste(format(x, digits = 12),
                                                   collapse = "_"), ""),
               collapse = "|")
  as.integer(sum(utf8ToInt(key) * seq_len(nchar(key))) %% .Machine$integer.max)
}

## --- population mean age profile --------------------------------------------
## flat_profile  : K equal-width bins over (0, t_max], equal expected counts.
##                 This is the old toy setup -- convenient, not realistic.
## empirical_profile : bin ages and mean profile taken from a real count matrix
##                 (or from the msprime singleton_counts.csv).
flat_profile <- function(K, t_max) {
  t_mid <- (seq_len(K) - 0.5) * t_max / K
  list(t = t_mid, p = rep(1 / K, K))
}

empirical_profile <- function(C, t_mid) {
  cm <- colMeans(C)
  list(t = t_mid, p = cm / sum(cm))
}

## Parse bin labels like "0-1", "1-2", "205+" into representative ages.
## Geometric midpoint for log-spaced bins; 1.5x the edge for the open top bin.
parse_bin_ages <- function(nm) {
  vapply(nm, function(s) {
    if (grepl("\\+$", s)) return(1.5 * as.numeric(sub("\\+$", "", s)))
    e <- as.numeric(strsplit(s, "-")[[1]])
    if (e[1] <= 0) mean(e) else sqrt(e[1] * e[2])
  }, numeric(1), USE.NAMES = FALSE)
}

## --- per-person age profiles ------------------------------------------------
## Each person gets a tilt u_i ~ N(0,1) towards younger or older singletons:
##     p_ik  proportional to  pbar_k * exp(theta * u_i * z_k)
## with z the standardised log age. To first order this shifts person i's mean
## log age by theta * u_i * sd(log t), so setting theta = cv_shape / sd(log t)
## makes SD(mean log age) = cv_shape exactly as intended.
##
## cv_shape = 0 returns the common profile for everyone (the "shared" model).
## Note the population mean profile shifts very slightly for large cv_shape
## (Jensen on the exponential); negligible below ~0.5.
make_profiles <- function(N, prof, cv_shape = SIM_DEFAULTS$cv_shape) {
  pbar <- prof$p
  if (cv_shape <= 0) return(matrix(pbar, N, length(pbar), byrow = TRUE))
  lt  <- log(prof$t)
  mu  <- sum(pbar * lt)
  sdl <- sqrt(sum(pbar * (lt - mu)^2))
  if (!is.finite(sdl) || sdl <= 0) return(matrix(pbar, N, length(pbar), byrow = TRUE))
  P <- sweep(exp(outer(rnorm(N) * cv_shape / sdl, (lt - mu) / sdl)), 2, pbar, `*`)
  P / rowSums(P)
}

## --- counts -----------------------------------------------------------------
## shared: one gamma factor per person scaling their whole profile (the burden),
##         times a per-person tilt of that profile (the shape). Total count is
##         unaffected by the tilt, so the two knobs are independent.
## indep : a fresh gamma factor per person PER BIN. Note this also changes the
##         total: CV(total) comes out at roughly cv_total / sqrt(K).
sim_counts <- function(N, prof, lam_tot,
                       cv_total    = SIM_DEFAULTS$cv_total,
                       cv_shape    = SIM_DEFAULTS$cv_shape,
                       count_model = SIM_DEFAULTS$count_model) {
  K  <- length(prof$t)
  sh <- 1 / cv_total^2

  if (count_model == "indep") {
    lam_k <- lam_tot * prof$p
    C <- vapply(seq_len(K),
                function(k) rpois(N, rgamma(N, sh, rate = sh / max(lam_k[k], 1e-12))),
                numeric(N))
  } else if (count_model == "shared") {
    lam_i <- rgamma(N, sh, rate = sh / lam_tot)        # knob 1: how many
    P     <- make_profiles(N, prof, cv_shape)          # knob 2: which ages
    C <- matrix(rpois(N * K, lam_i * P), N, K)
  } else stop("count_model must be 'shared' or 'indep'")

  colnames(C) <- paste0("t", round(prof$t, 1))
  C
}

## --- phenotype --------------------------------------------------------------
## var(g_i) = sigma_beta^2 * sum_k C_ik * decay(s, t_k), then residual noise
## scaled to give the requested singleton heritability.
make_pheno <- function(C, prof, S, H2,
                       sigma_beta = SIM_DEFAULTS$sigma_beta,
                       curve      = SIM_DEFAULTS$curve) {
  s2 <- sigma_beta^2 * decay(S, prof$t, curve)
  w  <- as.vector(C %*% s2)
  g  <- sqrt(w) * rnorm(nrow(C))
  list(Y = g + rnorm(nrow(C), 0, sqrt(var(g) * (1 - H2) / H2)),
       g = g, w = w, s2_true = s2)
}

## --- one dataset ------------------------------------------------------------
## The single entry point. Everything else above is a component of this.
sim_dataset <- function(S, REP, H2, N_inds, prof, lam_tot,
                        cv_total    = SIM_DEFAULTS$cv_total,
                        cv_shape    = SIM_DEFAULTS$cv_shape,
                        count_model = SIM_DEFAULTS$count_model,
                        curve       = SIM_DEFAULTS$curve,
                        sigma_beta  = SIM_DEFAULTS$sigma_beta) {
  set.seed(make_seed(S, REP, H2, N_inds, lam_tot, cv_total, cv_shape,
                     count_model, curve, length(prof$t), max(prof$t)))
  C  <- sim_counts(N_inds, prof, lam_tot, cv_total, cv_shape, count_model)
  ph <- make_pheno(C, prof, S, H2, sigma_beta, curve)
  c(list(C = C, t_mid = prof$t, prof = prof,
         pars = list(S = S, REP = REP, H2 = H2, N_inds = N_inds,
                     lam_tot = lam_tot, cv_total = cv_total,
                     cv_shape = cv_shape, count_model = count_model,
                     curve = curve)),
    ph)
}

## --- back-compatible wrapper for the v5 sweep -------------------------------
## Same signature and return value as the old sim_diag(), so existing scripts
## keep working; LAMBDA there was per BIN, so the total is LAMBDA * K.
sim_diag <- function(S, REP, H2, CV, N_inds, K, t_max, LAMBDA,
                     type = SIM_DEFAULTS$curve,
                     cv_shape = SIM_DEFAULTS$cv_shape,
                     count_model = SIM_DEFAULTS$count_model) {
  prof <- flat_profile(K, t_max)
  d <- sim_dataset(S, REP, H2, N_inds, prof, lam_tot = LAMBDA * K,
                   cv_total = CV, cv_shape = cv_shape,
                   count_model = count_model, curve = type)
  list(C = d$C, Y = d$Y, t_mid = d$t_mid, s2_true = d$s2_true)
}

## --- measure the two knobs in any count matrix ------------------------------
## Use this on real data to calibrate, and on simulated data to check the
## simulator reproduced what you asked for. The sampling term is estimated
## directly, so no reference simulation is needed.
measure_knobs <- function(C, t_mid) {
  Tt <- rowSums(C); ok <- Tt >= 2
  lt <- log(t_mid)
  L  <- as.vector(C %*% lt) / Tt                       # per-person mean log age
  s2 <- as.vector(C %*% lt^2) / Tt - L^2               # within-person spread
  var_tot  <- var(L[ok])
  var_samp <- mean(s2[ok] / (Tt[ok] - 1))
  data.frame(mean_total = mean(Tt),
             cv_total   = sd(Tt) / mean(Tt),
             var_logage_total    = var_tot,
             var_logage_sampling = var_samp,
             cv_shape   = sqrt(max(var_tot - var_samp, 0)))
}

## --- self-check -------------------------------------------------------------
## sim_check() simulates at the requested knob settings and measures them back.
## cv_total and cv_shape should come back close to what was asked for.
sim_check <- function(N = 2e5, K = 8, t_max = 200, lam_tot = 1000,
                      cv_total = 0.65, cv_shape = 0.33) {
  prof <- flat_profile(K, t_max)
  for (cs in unique(c(0, cv_shape))) {
    C <- sim_counts(N, prof, lam_tot, cv_total, cs, "shared")
    m <- measure_knobs(C, prof$t)
    cat(sprintf("asked cv_total %.2f cv_shape %.2f  ->  got %.3f / %.3f\n",
                cv_total, cs, m$cv_total, m$cv_shape))
  }
  invisible(NULL)
}
