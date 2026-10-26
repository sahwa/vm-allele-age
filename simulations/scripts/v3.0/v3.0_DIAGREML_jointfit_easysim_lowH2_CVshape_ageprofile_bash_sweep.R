#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## v3.3 sweep: joint REML only, one SLURM array task per cell.
##
##   Rscript v3_sweep.R grid                       write cells.txt, one per cell
##   Rscript v3_sweep.R cell lambda H2 CV  run all reps for one cell
##   Rscript v3_sweep.R collect                    combine, summarise, plot
##
## Changes from v3.2:
##   - two-stage (diagGREML + lm) dropped. It was 125x worse on RMSE and did not
##     track the truth at all; the comparison figure is frozen as supplementary.
##   - count model fixed at "shared" (one per-person factor scaling a common age
##     profile). "indep" was the old bug and is no longer swept.
##   - age-profile variation is now an explicit constant, CV_AGE.
##   - simulation and fit both use the HYPERBOLIC decay, per the SLiM result.
##   - the decay rate is called lambda throughout (was S). The per-bin singleton
##     count, previously LAMBDA, is now SING_PER_BIN to free the name.
##   - h2m_t0_hat is compared against its OWN truth, sigma2_m * f(lambda, t0),
##     not against the age-0 truth. Previously it looked biased when it was not.
## ---------------------------------------------------------------------------

library(data.table)
RNGkind("L'Ecuyer-CMRG")

VERSION <- "v3.4"
BASE <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations"
source(file.path(BASE, "scripts/common/fit_sel_reml_joint.R"))
source(file.path(BASE, "scripts/common/sim_diag.R"))     # sim_diag(), decay(), ...
source(file.path(BASE, "scripts/common/age_profile.R"))     # sim_diag(), decay(), ...

OUTDIR   <- file.path(BASE, "figs", VERSION)
DATDIR   <- file.path(BASE, "data", VERSION)
CELLDIR  <- file.path(DATDIR, "cells")
CELLFILE <- file.path(DATDIR, "cells.txt")
for (d in c(OUTDIR, CELLDIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- parameters ------------------------------------------------------------
SIGMA_BETA   <- 1
SING_PER_BIN <- 1000        # singletons per person PER BIN; total is SING_PER_BIN*K
                            # (UKB 150k WGS British/Irish: ~1330 per person total)

## Load up the age profile based on the 5e5 msprime sim

## Knob 1: how much people vary in HOW MANY singletons they carry -> signal for V_M.
##   0.64 = UKB WES full sample (mean SC 17.43, SD 11.14; Wang et al.)
##   0.37 = same sample with SC outliers removed
##   0.38 = n=500k msprime run
CV_vals <- c(0.65)

## Knob 2: how much people vary in WHICH AGES those are -> signal for lambda.
## Defined as the SD across people of their mean LOG singleton age, so
## CV_AGE = 0.33 means one SD above average carries singletons exp(0.33) = 1.4x
## older than typical.
##   0    = every profile identical (the fully "shared" model; pessimistic)
##   0.33 = MEASURED in the n=500k msprime run (total var of mean log age 0.172,
##          sampling 0.061, structural 0.111 -> sqrt = 0.333)
##   0.50 = generous, roughly what the GEL NFE/non-NFE age gap implies
CV_AGE      <- 0.50
COUNT_MODEL <- "shared"

## Simulation and fit must agree. sim_diag.R uses "hyperbolic"/"exponential";
## fit_sel_reml uses "hyp"/"exp".
CURVE_SIM <- "hyperbolic"
CURVE_FIT <- "hyp"

H2_vals     <- c(0.5, 0.1, 0.034, 0.02, 0.01)  # 0.5 = reference; 0.01-0.034 realistic
lambda_vals <- c(0, 1e-4, 1e-3, 1e-2, 3e-2)
N_INDS      <- 5e5
N_REPS      <- 100

PROF <- load_age_profile(n_bins = 7, scale_total = 1000)   # ~1000/person, genome-scale
stopifnot(all.equal(PROF$n_ind, N_INDS))   # ages must match the ascertainment n


H2_REALISTIC <- c(0.01, 0.02, 0.034, 0.05)
COL_EST   <- "#2C7BB6"
COL_TRUTH <- "#D7191C"

h2_lab <- function(x) {
  v <- as.numeric(as.character(x))
  ifelse(v %in% H2_REALISTIC, paste0("h² = ", v), paste0("h² = ", v, " (ref)"))
}

cell_grid <- function()
  CJ(lambda = lambda_vals, H2 = H2_vals, CV = CV_vals, sorted = FALSE)
cell_tag <- function(lambda, H2, CV)
  sprintf("lam%g_H%g_CV%g", lambda, H2, CV)
write_atomic <- function(x, f) { tmp <- paste0(f, ".tmp"); fwrite(x, tmp); file.rename(tmp, f) }

## ---- fit one replicate -----------------------------------------------------
run_one <- function(lambda, REP, H2, CV, N_inds) {
  setDTthreads(1)
  d <- sim_dataset_real(lambda = lambda, REP = REP, H2 = H2, N_inds = N_inds,
                        prof = PROF, count_model = COUNT_MODEL, curve = CURVE_SIM)
  C <- d$C; m <- colMeans(C); t_mid <- d$t_mid
  if (any(!is.finite(m) | m <= 0)) stop("bad mean count")

  jf <- tryCatch(fit_sel_reml(y = d$Y, C = C, t = t_mid, curve = CURVE_FIT),
                 error = function(e) NULL)

  out <- data.table(lambda = lambda, REP = REP, H2 = H2, CV = CV, cv_age = CV_AGE,
                    count = COUNT_MODEL, N_inds = N_inds, K = length(t_mid),
                    V_P = d$V_P, mu_per_gen = d$u_per_gen, t0_dat = d$t0_dat,
                    h2m_true = d$h2m_true, h2m_t0_true = d$h2m_t0_true,
                    sigma2_m = NA_real_, sigma2_t0 = NA_real_, t0 = NA_real_,
                    h2m_hat = NA_real_, h2m_t0_hat = NA_real_,
                    lambda_hat = NA_real_, lambda_se = NA_real_, lrt_p = NA_real_,
                    h2_sing = NA_real_, sigma2_e = NA_real_, mu_hat = NA_real_,
                    loglik = NA_real_, conv = NA_integer_, at_bound = NA)

  if (!is.null(jf)) {
    out[, `:=`(sigma2_m = jf$sigma2_m, sigma2_t0 = jf$sigma2_t0, t0 = jf$t0,
               h2m_hat    = d$u_per_gen * jf$sigma2_m  / d$V_P,
               h2m_t0_hat = d$u_per_gen * jf$sigma2_t0 / d$V_P,
               lambda_hat = jf$s, lambda_se = jf$s_se, lrt_p = jf$lrt_p,
               h2_sing = jf$h2_sing, sigma2_e = jf$sigma2_e, mu_hat = jf$b[1],
               loglik = jf$logLik, conv = as.integer(jf$convergence),
               at_bound = jf$at_bound)]
  }
  out
}

## ---- dispatch --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else ""

if (mode == "grid") {

  g <- cell_grid()
  fwrite(g, CELLFILE, sep = " ", col.names = FALSE)
  message(nrow(g), " cells written to ", CELLFILE,
          "\ncolumns: lambda H2 CV",
          "\nfixed:   cv_age = ", CV_AGE, ", count = ", COUNT_MODEL,
          ", curve = ", CURVE_SIM,
          "\nsubmit:  sbatch --array=1-", nrow(g), " v3_submit.sh")

} else if (mode == "cell") {

  a <- suppressWarnings(as.numeric(args[-1]))
  if (length(a) != 3 || anyNA(a)) stop("usage: v3_sweep.R cell lambda H2 CV")
  lambda <- lambda_vals[which.min(abs(lambda_vals - a[1]))]   # snap to the grid value
  if (abs(lambda - a[1]) > 1e-12) stop("lambda must be one of lambda_vals")
  H2 <- a[2]; CV <- a[3]

  tag   <- cell_tag(lambda, H2, CV)
  f_rep <- file.path(CELLDIR, paste0("reps_", tag, ".csv"))
  if (file.exists(f_rep)) { message(tag, ": already done, skipping"); quit(save = "no") }

  tic <- Sys.time()
  res <- lapply(seq_len(N_REPS), function(REP)
    tryCatch(run_one(lambda, REP, H2, CV, N_INDS),
             error = function(e) { message("rep ", REP, " failed: ", conditionMessage(e)); NULL }))
  ok <- Filter(Negate(is.null), res)
  if (!length(ok)) stop("all reps failed")

  write_atomic(rbindlist(ok, use.names = TRUE), f_rep)
  message(sprintf("%s: %d/%d reps in %.1f min", tag, length(ok), N_REPS,
                  as.numeric(difftime(Sys.time(), tic, units = "mins"))))

} else if (mode == "collect") {

  library(ggplot2); library(scales)

  cells <- cell_grid()
  tags  <- cells[, cell_tag(lambda, H2, CV)]
  f_rep <- file.path(CELLDIR, paste0("reps_", tags, ".csv"))
  have  <- file.exists(f_rep)
  if (!any(have)) stop("no finished cells in ", CELLDIR)
  if (!all(have)) warning(sum(!have), " cells missing. Rerun with: sbatch --array=",
                          paste(which(!have), collapse = ","), " v3_submit.sh")

  reps  <- rbindlist(lapply(f_rep[have], fread), use.names = TRUE)
  short <- reps[, .N, by = .(lambda, H2, CV)][N < N_REPS]
  if (nrow(short)) { message("Cells with failed reps:"); print(short) }

  tag <- paste0("N", N_INDS)
  fwrite(reps, file.path(DATDIR, paste0("reps_", tag, ".csv")))

  ## Detection by the LRT, which is defined even for boundary fits. The CI rule
  ## needs a finite lambda_se, which is NA at the bound -- so it silently drops
  ## exactly the replicates where lambda was hardest to estimate.
  reps[, converged   := conv %in% 0]
  reps[, detected    := is.finite(lrt_p) & lrt_p < 0.05 & lambda_hat > 0]
  reps[, ci_lo       := lambda_hat - 1.96 * lambda_se]
  reps[, detected_ci := is.finite(ci_lo) & ci_lo > 0]

  ## ---- per-cell summary ----------------------------------------------------
  summ <- reps[, {
    lh  <- lambda_hat[is.finite(lambda_hat)]; z <- lh - lambda[1]
    okh <- is.finite(h2m_hat)
    rh  <- h2m_hat[okh]    / h2m_true[okh]        # age-0 target vs age-0 truth
    rh0 <- h2m_t0_hat[okh] / h2m_t0_true[okh]     # centroid target vs centroid truth
    .(n_total = .N, n_fit = length(lh), fail_rate = 1 - length(lh) / .N,
      frac_bound = mean(at_bound, na.rm = TRUE),
      frac_conv  = mean(converged, na.rm = TRUE),
      h2m_true   = mean(h2m_true),
      mean_lambda   = if (length(lh)) mean(lh) else NA_real_,
      median_lambda = if (length(lh)) median(lh) else NA_real_,
      bias      = if (length(z)) mean(z) else NA_real_,
      se_bias   = if (length(z) >= 2) sd(z) / sqrt(length(z)) else NA_real_,
      sd_lambda = if (length(lh) >= 2) sd(lh) else NA_real_,
      rmse      = if (length(z)) sqrt(mean(z^2)) else NA_real_,
      rel_prec  = if (length(lh) >= 2 && lambda[1] > 0) sd(lh) / lambda[1] else NA_real_,
      power     = mean(detected, na.rm = TRUE),
      power_ci  = mean(detected_ci, na.rm = TRUE),
      h2m_med_ratio    = if (length(rh))  median(rh)  else NA_real_,
      h2m_t0_med_ratio = if (length(rh0)) median(rh0) else NA_real_,
      h2m_rel_bias = if (length(rh)) mean(rh) - 1 else NA_real_,
      gmed_h2m = if (length(rh)) median(h2m_hat[okh]) else NA_real_,
      fold_lo  = if (length(rh)) quantile(h2m_hat[okh], 0.025) else NA_real_,
      fold_hi  = if (length(rh)) quantile(h2m_hat[okh], 0.975) else NA_real_)
  }, by = .(lambda, H2, CV, N_inds)]

  setorder(summ, CV, -H2, lambda)
  fwrite(summ, file.path(DATDIR, paste0("summary_", tag, ".csv")))

  ## ---- figures -------------------------------------------------------------
  for (d in list(reps, summ)) {
    d[, H2_f     := factor(H2, levels = H2_vals, labels = h2_lab(H2_vals))]
    d[, CV_f     := factor(CV, levels = CV_vals, labels = paste0("CV = ", CV_vals))]
    d[, lambda_f := factor(lambda, levels = lambda_vals)]
  }

  th <- theme_bw() +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
          axis.text = element_text(size = 14), axis.title = element_text(size = 14),
          strip.text = element_text(size = 14))
  dg2 <- position_dodge(width = 0.5)

  lambda_select <- c(0.001, 0.01)
  H2_select     <- c(0.034, 0.02, 0.01)

  r  <- reps[lambda %in% lambda_select & H2 %in% H2_select]
  su <- summ[lambda %in% lambda_select & H2 %in% H2_select]
  if (!nrow(r)) stop("nothing left after the lambda/H2 filter")

  sfx <- tag
  ## K and t_max are no longer parameters -- both come from the age profile, so
  ## t0 (the count-weighted mean age) is the summary worth printing.
  ttl <- sprintf("N = %s, %d age bins, t0 = %.1f gen, CV_age = %g",
                 format(N_INDS, big.mark = ",", scientific = FALSE),
                 length(PROF$t), PROF$t0, CV_AGE)

  ann_tab <- function(dt, est, truth) {
    a <- dt[is.finite(get(est)),
            .(bias = mean(get(est) - get(truth)),
              rmse = sqrt(mean((get(est) - get(truth))^2)), n = .N),
            by = .(lambda_f, CV_f, H2_f)]
    a[, lab := sprintf("bias %s\nRMSE %s",
                       formatC(bias, format = "g", digits = 2, flag = "+"),
                       formatC(rmse, format = "g", digits = 2))]
    a[]
  }
  ann_l   <- ann_tab(r, "lambda_hat", "lambda")
  ann_v   <- ann_tab(r, "h2m_hat",    "h2m_true")
  tru_lam <- unique(r[, .(lambda_f, lambda, CV_f, H2_f)])

  ## (i) estimates of lambda
  p1 <- ggplot(r[is.finite(lambda_hat)], aes(lambda_f, lambda_hat)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_jitter(width = 0.12, size = 1, alpha = 0.25, colour = COL_EST) +
    geom_boxplot(outlier.shape = NA, fill = NA, width = 0.6, colour = COL_EST) +
    geom_errorbar(data = tru_lam, aes(x = lambda_f, ymin = lambda, ymax = lambda),
                  width = 0.9, colour = COL_TRUTH, linewidth = 0.8, inherit.aes = FALSE) +
    geom_text(data = ann_l, aes(lambda_f, Inf, label = lab), vjust = 1.2,
              size = 2.6, lineheight = 0.9, inherit.aes = FALSE) +
    facet_grid(CV_f ~ H2_f) +
    scale_y_continuous(trans = pseudo_log_trans(sigma = 2e-3),
                       breaks = c(-0.1, -0.01, -0.001, 0, 0.001, 0.01, 0.03, 0.1)) +
    labs(x = expression("True " * lambda), y = expression(hat(lambda)),
         title = paste0("Estimates of \u03bb: ", ttl),
         subtitle = "Pseudo-log y axis; red bar = truth; points at the extremes are fits at the search bound") + th

  ## (ii) bias in lambda
  p2 <- ggplot(su, aes(lambda_f, bias)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_linerange(aes(ymin = bias - sd_lambda, ymax = bias + sd_lambda),
                   alpha = 0.3, linewidth = 2.5, colour = COL_EST) +
    geom_pointrange(aes(ymin = bias - 1.96 * se_bias, ymax = bias + 1.96 * se_bias),
                    size = 0.35, colour = COL_EST) +
    facet_grid(CV_f ~ H2_f, scales = "free_y") +
    labs(x = expression("True " * lambda), y = expression(hat(lambda) - lambda),
         title = paste0("Bias in \u03bb\u0302: ", ttl),
         subtitle = "Dark: 95% CI on the mean. Pale: SD across replicates (what one study faces).") + th

  ## (iii) power. LRT keeps boundary fits, the CI rule drops them; the gap
  ## between the lines is how much that choice matters.
  pw <- melt(su, id.vars = c("lambda_f", "CV_f", "H2_f"),
             measure.vars = c("power", "power_ci"),
             variable.name = "rule", value.name = "p")
  pw[, rule := factor(rule, c("power", "power_ci"), c("LRT p < 0.05", "95% CI excludes 0"))]
  p3 <- ggplot(pw, aes(lambda_f, p, colour = rule, group = rule)) +
    geom_hline(yintercept = 0.05, colour = "grey60", linetype = "dotted") +
    geom_hline(yintercept = 0.80, colour = COL_TRUTH, linetype = "dashed") +
    geom_line(alpha = 0.6) + geom_point(size = 2) +
    facet_grid(CV_f ~ H2_f) + scale_y_continuous(limits = c(0, 1)) +
    scale_colour_manual(values = c("LRT p < 0.05" = COL_EST,
                                   "95% CI excludes 0" = "grey35"), name = NULL) +
    labs(x = expression("True " * lambda), y = expression("P(detect " * lambda > 0 * ")"),
         title = paste0("Power to detect \u03bb > 0: ", ttl),
         subtitle = "At \u03bb = 0 this is the false-positive rate (dotted 0.05); dashed red = 80%") + th

  ## (iv) mutational heritability
  tru <- unique(r[, .(lambda_f, CV_f, H2_f, h2m_true)])
  p4 <- ggplot(r[is.finite(h2m_hat)], aes(lambda_f, h2m_hat)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_jitter(width = 0.12, size = 1, alpha = 0.25, colour = COL_EST) +
    geom_boxplot(outlier.shape = NA, fill = NA, width = 0.6, colour = COL_EST) +
    geom_point(data = tru, aes(lambda_f, h2m_true), colour = COL_TRUTH, size = 2.5,
               inherit.aes = FALSE) +
    geom_text(data = ann_v, aes(lambda_f, Inf, label = lab), vjust = 1.2,
              size = 2.6, lineheight = 0.9, inherit.aes = FALSE) +
    facet_grid(CV_f ~ H2_f) +
    scale_y_continuous(trans = pseudo_log_trans(sigma = 1e-5),
                       breaks = c(-1e-3, -1e-4, 0, 1e-4, 1e-3, 1e-2),
                       labels = function(x) format(x, scientific = TRUE)) +
    coord_cartesian(ylim = c(-3e-3, 3e-2)) +
    labs(x = expression("True " * lambda), y = expression(hat(V)[M] / V[P]),
         title = paste0("Mutational heritability: ", ttl),
         subtitle = "Red = truth; pseudo-log axis, includes zero and negative estimates; view clipped") + th

  ## (v) V_M bias: age-0 target vs centroid target, each against ITS OWN truth.
  ## sigma2_m = sigma2_t0 * (1 + lambda*t0), so error in lambda_hat is amplified
  ## by the extrapolation distance; the centroid value is not.
  rb <- melt(su, id.vars = c("lambda_f", "CV_f", "H2_f"),
             measure.vars = c("h2m_med_ratio", "h2m_t0_med_ratio"),
             variable.name = "target", value.name = "ratio")
  rb[, target := factor(target, c("h2m_med_ratio", "h2m_t0_med_ratio"),
                        c("extrapolated to age 0", "at centroid age t0"))]
  p5 <- ggplot(rb, aes(lambda_f, ratio, colour = target, group = target)) +
    geom_hline(yintercept = 1, colour = COL_TRUTH, linetype = "dashed") +
    geom_line(alpha = 0.6) + geom_point(size = 2, position = dg2) +
    facet_grid(CV_f ~ H2_f) + scale_y_log10() +
    scale_colour_manual(values = c("extrapolated to age 0" = COL_EST,
                                   "at centroid age t0" = "grey35"), name = NULL) +
    labs(x = expression("True " * lambda), y = "estimate / truth (median ratio)",
         title = paste0("Bias in mutational heritability: ", ttl),
         subtitle = "Each series against its own truth. Dashed red = unbiased. The centroid value needs no back-extrapolation.") + th

  W <- 14; Hh <- 7.5
  ggsave(file.path(OUTDIR, paste0("lambda_estimates_", sfx, ".png")), p1, width = W, height = Hh, dpi = 150)
  ggsave(file.path(OUTDIR, paste0("lambda_bias_",      sfx, ".png")), p2, width = W, height = Hh, dpi = 150)
  ggsave(file.path(OUTDIR, paste0("lambda_power_",     sfx, ".png")), p3, width = W, height = Hh, dpi = 150)
  ggsave(file.path(OUTDIR, paste0("h2m_",              sfx, ".png")), p4, width = W, height = Hh, dpi = 150)
  ggsave(file.path(OUTDIR, paste0("h2m_bias_",         sfx, ".png")), p5, width = W, height = Hh, dpi = 150)

  ## ---- headline numbers ----------------------------------------------------
  cat("\n-- false-positive rate at lambda = 0 (should be ~0.05) --\n")
  print(summ[lambda == 0, .(H2, power, power_ci, frac_bound, fail_rate)][
        order(-H2)])

  cat("\n-- lambda = 0.01 --\n")
  print(summ[lambda == 0.01, .(H2, power, bias, sd_lambda, rel_prec,
                               h2m_med_ratio, h2m_t0_med_ratio, frac_bound)][
        order(-H2)])

} else stop("usage: v3_sweep.R grid | cell lambda H2 CV | collect")
