#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## v3.3 sweep: joint REML only, one SLURM array task per cell.
##
##   Rscript v3_sweep.R grid                     write cells.txt, one line per cell
##   Rscript v3_sweep.R cell S H2 CV K t_max     run all reps for one cell
##   Rscript v3_sweep.R collect                  combine cells, summarise, plot
##
## Changes from v3.2:
##   - two-stage (diagGREML + lm) dropped. It was 125x worse on RMSE and did not
##     track the truth at all; the comparison figure is frozen as supplementary.
##   - count model fixed at "shared" (one per-person factor scaling a common age
##     profile). "indep" was the old bug and is no longer swept.
##   - age-profile variation is now an explicit constant, CV_SHAPE, rather than
##     being implicitly zero.
##   - simulation and fit both use the HYPERBOLIC decay, per the SLiM result.
## ---------------------------------------------------------------------------

library(data.table)
RNGkind("L'Ecuyer-CMRG")

VERSION <- "v3.3"
BASE <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations"
source(file.path(BASE, "scripts/common/fit_sel_reml_joint.R"))
source(file.path(BASE, "scripts/common/sim_diag.R"))     # sim_diag(), decay(), ...

OUTDIR   <- file.path(BASE, "figs", VERSION)
DATDIR   <- file.path(BASE, "data", VERSION)
CELLDIR  <- file.path(DATDIR, "cells")
CELLFILE <- file.path(DATDIR, "cells.txt")
for (d in c(OUTDIR, CELLDIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- parameters ------------------------------------------------------------
SIGMA_BETA <- 1
LAMBDA     <- 1000          # singletons per person PER BIN, so the total is LAMBDA*K
                            # (UKB 150k WGS British/Irish: ~1330 per person total)

## Knob 1: how much people vary in HOW MANY singletons they carry -> signal for V_M.
##   0.64 = UKB WES full sample (mean SC 17.43, SD 11.14; Wang et al.)
##   0.37 = same sample with SC outliers removed
##   0.38 = n=500k msprime run
CV_vals <- c(0.65)

## Knob 2: how much people vary in WHICH AGES those are -> signal for s.
## Defined as the SD across people of their mean LOG singleton age.
##   0    = every profile identical (the old "shared" model; pessimistic)
##   0.33 = MEASURED in the n=500k msprime run (total var of mean log age 0.172,
##          sampling 0.061, structural 0.111 -> sqrt = 0.333)
##   0.50 = generous, roughly what the GEL NFE/non-NFE age gap implies
## NB this is a different quantity from CV_vals above; 0.65 here would be about
## 4x the measured signal for s.
CV_SHAPE    <- 0.65
COUNT_MODEL <- "shared"

## Simulation and fit must agree. sim_common uses "hyperbolic"/"exponential";
## fit_sel_reml uses "hyp"/"exp".
CURVE_SIM <- "hyperbolic"
CURVE_FIT <- "hyp"

H2_vals <- c(0.5, 0.1, 0.034, 0.02, 0.01)   # 0.5 = reference; 0.01-0.034 realistic (Wang et al.)
Svals   <- c(0, 1e-4, 1e-3, 1e-2, 3e-2)
N_INDS  <- 5e5
K_VALS  <- c(3, 8)
T_VALS  <- c(200, 400)
N_REPS  <- 100

H2_REALISTIC <- c(0.01, 0.02, 0.034, 0.05)
COL_EST   <- "#2C7BB6"
COL_TRUTH <- "#D7191C"

h2_lab <- function(x) {
  v <- as.numeric(as.character(x))
  ifelse(v %in% H2_REALISTIC, paste0("h² = ", v), paste0("h² = ", v, " (ref)"))
}

cell_grid <- function()
  CJ(S = Svals, H2 = H2_vals, CV = CV_vals, K = K_VALS, t_max = T_VALS, sorted = FALSE)
cell_tag <- function(S, H2, CV, K, t_max)
  sprintf("S%g_H%g_CV%g_K%g_t%g", S, H2, CV, K, t_max)
write_atomic <- function(x, f) { tmp <- paste0(f, ".tmp"); fwrite(x, tmp); file.rename(tmp, f) }

## ---- fit one replicate -----------------------------------------------------
run_one <- function(S, REP, H2, CV, N_inds, K, t_max) {
  setDTthreads(1)
  d <- sim_diag(S = S, REP = REP, H2 = H2, CV = CV, N_inds = N_inds,
                K = K, t_max = t_max, LAMBDA = LAMBDA,
                type = CURVE_SIM, cv_shape = CV_SHAPE, count_model = COUNT_MODEL)
  C <- d$C; m <- colMeans(C); t_mid <- d$t_mid
  if (any(!is.finite(m) | m <= 0)) stop("bad mean count")

  ## V_M = (new mutations per person per generation) * sigma2_m.
  ## Bin width is t_max/K, so mutations per generation = LAMBDA * K / t_max.
  mu_per_gen <- LAMBDA * K / t_max
  V_P        <- var(d$Y)
  h2m_true   <- mu_per_gen * SIGMA_BETA^2 / V_P

  jf <- tryCatch(fit_sel_reml(y = d$Y, C = C, t = t_mid, curve = CURVE_FIT),
                 error = function(e) NULL)

  out <- data.table(S = S, REP = REP, H2 = H2, CV = CV, cv_shape = CV_SHAPE,
                    count = COUNT_MODEL, N_inds = N_inds, K = K, t_max = t_max,
                    V_P = V_P, mu_per_gen = mu_per_gen, h2m_true = h2m_true,
                    sigma2_m = NA_real_, sigma2_t0 = NA_real_, t0 = NA_real_,
                    h2m_hat = NA_real_, h2m_t0_hat = NA_real_,
                    s_hat = NA_real_, s_se = NA_real_, lrt_p = NA_real_,
                    h2_sing = NA_real_, sigma2_e = NA_real_, mu_hat = NA_real_,
                    loglik = NA_real_, conv = NA_integer_, at_bound = NA)

  if (!is.null(jf)) {
    out[, `:=`(sigma2_m   = jf$sigma2_m,
               sigma2_t0  = jf$sigma2_t0,
               t0         = jf$t0,
               h2m_hat    = mu_per_gen * jf$sigma2_m / V_P,
               ## sigma2_t0 needs no back-extrapolation, so it is immune to the
               ## exp(s*t0) amplification of any error in s_hat. Worth tracking.
               h2m_t0_hat = mu_per_gen * jf$sigma2_t0 / V_P,
               s_hat      = jf$s,
               s_se       = jf$s_se,
               lrt_p      = jf$lrt_p,
               h2_sing    = jf$h2_sing,
               sigma2_e   = jf$sigma2_e,
               mu_hat     = jf$b[1],
               loglik     = jf$logLik,
               conv       = as.integer(jf$convergence),
               at_bound   = jf$at_bound)]
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
          "\ncolumns: S H2 CV K t_max",
          "\nfixed:   cv_shape = ", CV_SHAPE, ", count = ", COUNT_MODEL,
          ", curve = ", CURVE_SIM,
          "\nsubmit:  sbatch --array=1-", nrow(g), " v3_submit.sh")

} else if (mode == "cell") {

  a <- suppressWarnings(as.numeric(args[-1]))
  if (length(a) != 5 || anyNA(a)) stop("usage: v3_sweep.R cell S H2 CV K t_max")
  S <- Svals[which.min(abs(Svals - a[1]))]        # snap to the grid value
  if (abs(S - a[1]) > 1e-12) stop("S must be one of Svals")
  H2 <- a[2]; CV <- a[3]; K <- as.integer(a[4]); t_max <- a[5]

  tag   <- cell_tag(S, H2, CV, K, t_max)
  f_rep <- file.path(CELLDIR, paste0("reps_", tag, ".csv"))
  if (file.exists(f_rep)) { message(tag, ": already done, skipping"); quit(save = "no") }

  t0  <- Sys.time()
  res <- lapply(seq_len(N_REPS), function(REP)
    tryCatch(run_one(S, REP, H2, CV, N_INDS, K, t_max),
             error = function(e) { message("rep ", REP, " failed: ", conditionMessage(e)); NULL }))
  ok <- Filter(Negate(is.null), res)
  if (!length(ok)) stop("all reps failed")

  write_atomic(rbindlist(ok, use.names = TRUE), f_rep)
  message(sprintf("%s: %d/%d reps in %.1f min", tag, length(ok), N_REPS,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))

} else if (mode == "collect") {

  library(ggplot2); library(scales)

  cells <- cell_grid()
  tags  <- cells[, cell_tag(S, H2, CV, K, t_max)]
  f_rep <- file.path(CELLDIR, paste0("reps_", tags, ".csv"))
  have  <- file.exists(f_rep)
  if (!any(have)) stop("no finished cells in ", CELLDIR)
  if (!all(have)) warning(sum(!have), " cells missing. Rerun with: sbatch --array=",
                          paste(which(!have), collapse = ","), " v3_submit.sh")

  reps  <- rbindlist(lapply(f_rep[have], fread), use.names = TRUE)
  short <- reps[, .N, by = .(S, H2, CV, K, t_max)][N < N_REPS]
  if (nrow(short)) { message("Cells with failed reps:"); print(short) }

  tag <- paste0("N", N_INDS)
  fwrite(reps, file.path(DATDIR, paste0("reps_", tag, ".csv")))

  ## Detection by the LRT, which is defined even for boundary fits. The old
  ## CI rule needed a finite s_se, which is NA at the bound -- so it silently
  ## dropped exactly the replicates where s was hardest to estimate.
  reps[, converged := conv %in% 0]
  reps[, detected  := is.finite(lrt_p) & lrt_p < 0.05 & s_hat > 0]
  reps[, ci_lo     := s_hat - 1.96 * s_se]
  reps[, detected_ci := is.finite(ci_lo) & ci_lo > 0]

  ## ---- per-cell summary ----------------------------------------------------
  summ <- reps[, {
    sh  <- s_hat[is.finite(s_hat)]; z <- sh - S
    okh <- is.finite(h2m_hat)
    rh  <- h2m_hat[okh] / h2m_true[okh]
    rh0 <- h2m_t0_hat[okh] / h2m_true[okh]
    .(n_total = .N, n_fit = length(sh), fail_rate = 1 - length(sh) / .N,
      frac_bound = mean(at_bound, na.rm = TRUE),
      frac_conv  = mean(converged, na.rm = TRUE),
      h2m_true   = mean(h2m_true),
      mean_s = if (length(sh)) mean(sh) else NA_real_,
      median_s = if (length(sh)) median(sh) else NA_real_,
      bias    = if (length(z)) mean(z) else NA_real_,
      se_bias = if (length(z) >= 2) sd(z) / sqrt(length(z)) else NA_real_,
      sd_s    = if (length(sh) >= 2) sd(sh) else NA_real_,
      rmse    = if (length(z)) sqrt(mean(z^2)) else NA_real_,
      rel_prec = if (length(sh) >= 2 && S > 0) sd(sh) / S else NA_real_,
      power      = mean(detected, na.rm = TRUE),
      power_ci   = mean(detected_ci, na.rm = TRUE),
      h2m_med_ratio    = if (length(rh))  median(rh)  else NA_real_,
      h2m_t0_med_ratio = if (length(rh0)) median(rh0) else NA_real_,
      h2m_rel_bias = if (length(rh)) mean(rh) - 1 else NA_real_,
      gmed_h2m = if (length(rh)) median(h2m_hat[okh]) else NA_real_,
      fold_lo  = if (length(rh)) quantile(h2m_hat[okh], 0.025) else NA_real_,
      fold_hi  = if (length(rh)) quantile(h2m_hat[okh], 0.975) else NA_real_)
  }, by = .(S, H2, CV, N_inds, K, t_max)]

  setorder(summ, K, t_max, CV, -H2, S)
  fwrite(summ, file.path(DATDIR, paste0("summary_", tag, ".csv")))

  ## ---- figures -------------------------------------------------------------
  for (d in list(reps, summ)) {
    d[, H2_f := factor(H2, levels = H2_vals, labels = h2_lab(H2_vals))]
    d[, CV_f := factor(CV, levels = CV_vals, labels = paste0("CV = ", CV_vals))]
    d[, S_f  := factor(S,  levels = Svals)]
  }

  th <- theme_bw() +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
          axis.text = element_text(size = 14), axis.title = element_text(size = 14),
          strip.text = element_text(size = 14))
  dg2 <- position_dodge(width = 0.5)

  S_select  <- c(0.001, 0.01)
  H2_select <- c(0.034, 0.02, 0.01)

  for (kk in K_VALS) for (tt in T_VALS) {
    r <- reps[K == kk & t_max == tt & S %in% S_select & H2 %in% H2_select]
    s <- summ[K == kk & t_max == tt & S %in% S_select & H2 %in% H2_select]
    if (!nrow(r)) next
    sfx <- paste0(tag, "_K", kk, "_t", tt)
    ttl <- sprintf("N = %s, %d bins, ages 0-%d gen, cv_shape = %g",
                   format(N_INDS, big.mark = ",", scientific = FALSE), kk, tt, CV_SHAPE)

    ann_tab <- function(dt, est, truth) {
      a <- dt[is.finite(get(est)),
              .(bias = mean(get(est) - get(truth)),
                rmse = sqrt(mean((get(est) - get(truth))^2)), n = .N),
              by = .(S_f, CV_f, H2_f)]
      a[, lab := sprintf("bias %s\nRMSE %s",
                         formatC(bias, format = "g", digits = 2, flag = "+"),
                         formatC(rmse, format = "g", digits = 2))]
      a[]
    }
    ann_s <- ann_tab(r, "s_hat",   "S")
    ann_v <- ann_tab(r, "h2m_hat", "h2m_true")
    tru_s <- unique(r[, .(S_f, S, CV_f, H2_f)])

    ## (i) estimates of s
    p1 <- ggplot(r[is.finite(s_hat)], aes(S_f, s_hat)) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_jitter(width = 0.12, size = 1, alpha = 0.25, colour = COL_EST) +
      geom_boxplot(outlier.shape = NA, fill = NA, width = 0.6, colour = COL_EST) +
      geom_errorbar(data = tru_s, aes(x = S_f, ymin = S, ymax = S), width = 0.9,
                    colour = COL_TRUTH, linewidth = 0.8, inherit.aes = FALSE) +
      geom_text(data = ann_s, aes(S_f, Inf, label = lab), vjust = 1.2,
                size = 2.6, lineheight = 0.9, inherit.aes = FALSE) +
      facet_grid(CV_f ~ H2_f) +
      scale_y_continuous(trans = pseudo_log_trans(sigma = 2e-3),
                         breaks = c(-0.1, -0.01, -0.001, 0, 0.001, 0.01, 0.03, 0.1)) +
      labs(x = "True s", y = expression(hat(s)), title = paste0("Estimates of s: ", ttl),
           subtitle = "Pseudo-log y axis; red bar = truth; points at the extremes are fits at the search bound") + th

    ## (ii) bias in s
    p2 <- ggplot(s, aes(S_f, bias)) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_linerange(aes(ymin = bias - sd_s, ymax = bias + sd_s),
                     alpha = 0.3, linewidth = 2.5, colour = COL_EST) +
      geom_pointrange(aes(ymin = bias - 1.96 * se_bias, ymax = bias + 1.96 * se_bias),
                      size = 0.35, colour = COL_EST) +
      facet_grid(CV_f ~ H2_f, scales = "free_y") +
      labs(x = "True s", y = expression(hat(s) - s), title = paste0("Bias in ŝ: ", ttl),
           subtitle = "Dark: 95% CI on the mean. Pale: SD across replicates (what one study faces).") + th

    ## (iii) power. Two detection rules: LRT keeps boundary fits, the CI rule
    ## drops them -- the gap between the lines is how much that choice matters.
    pw <- melt(s, id.vars = c("S_f", "CV_f", "H2_f"),
               measure.vars = c("power", "power_ci"),
               variable.name = "rule", value.name = "p")
    pw[, rule := factor(rule, c("power", "power_ci"), c("LRT p < 0.05", "95% CI excludes 0"))]
    p3 <- ggplot(pw, aes(S_f, p, colour = rule, group = rule)) +
      geom_hline(yintercept = 0.05, colour = "grey60", linetype = "dotted") +
      geom_hline(yintercept = 0.80, colour = COL_TRUTH, linetype = "dashed") +
      geom_line(alpha = 0.6) + geom_point(size = 2) +
      facet_grid(CV_f ~ H2_f) + scale_y_continuous(limits = c(0, 1)) +
      scale_colour_manual(values = c("LRT p < 0.05" = COL_EST,
                                     "95% CI excludes 0" = "grey35"), name = NULL) +
      labs(x = "True s", y = "P(detect s > 0)", title = paste0("Power to detect s > 0: ", ttl),
           subtitle = "At s = 0 this is the false-positive rate (dotted 0.05); dashed red = 80%") + th

    ## (iv) mutational heritability
    tru <- unique(r[, .(S_f, CV_f, H2_f, h2m_true)])
    p4 <- ggplot(r[is.finite(h2m_hat)], aes(S_f, h2m_hat)) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_jitter(width = 0.12, size = 1, alpha = 0.25, colour = COL_EST) +
      geom_boxplot(outlier.shape = NA, fill = NA, width = 0.6, colour = COL_EST) +
      geom_point(data = tru, aes(S_f, h2m_true), colour = COL_TRUTH, size = 2.5,
                 inherit.aes = FALSE) +
      geom_text(data = ann_v, aes(S_f, Inf, label = lab), vjust = 1.2,
                size = 2.6, lineheight = 0.9, inherit.aes = FALSE) +
      facet_grid(CV_f ~ H2_f) +
      scale_y_continuous(trans = pseudo_log_trans(sigma = 1e-5),
                         breaks = c(-1e-3, -1e-4, 0, 1e-4, 1e-3, 1e-2),
                         labels = function(x) format(x, scientific = TRUE)) +
      coord_cartesian(ylim = c(-3e-3, 3e-2)) +
      labs(x = "True s", y = expression(hat(V)[M] / V[P]),
           title = paste0("Mutational heritability: ", ttl),
           subtitle = "Red = truth; pseudo-log axis, includes zero and negative estimates; view clipped") + th

    ## (v) V_M bias, extrapolated to age 0 vs evaluated at the centroid t0.
    ## sigma2_m = sigma2_t0 * (1 + s*t0), so any error in s_hat is amplified by
    ## the extrapolation distance. sigma2_t0 should be the more robust of the two.
    rb <- melt(s, id.vars = c("S_f", "CV_f", "H2_f"),
               measure.vars = c("h2m_med_ratio", "h2m_t0_med_ratio"),
               variable.name = "target", value.name = "ratio")
    rb[, target := factor(target, c("h2m_med_ratio", "h2m_t0_med_ratio"),
                          c("extrapolated to age 0", "at centroid age t0"))]
    p5 <- ggplot(rb, aes(S_f, ratio, colour = target, group = target)) +
      geom_hline(yintercept = 1, colour = COL_TRUTH, linetype = "dashed") +
      geom_line(alpha = 0.6) + geom_point(size = 2, position = dg2) +
      facet_grid(CV_f ~ H2_f) + scale_y_log10() +
      scale_colour_manual(values = c("extrapolated to age 0" = COL_EST,
                                     "at centroid age t0" = "grey35"), name = NULL) +
      labs(x = "True s", y = "estimate / truth (median ratio)",
           title = paste0("Bias in mutational heritability: ", ttl),
           subtitle = "Dashed red = unbiased. The centroid value needs no back-extrapolation, so it is not amplified by error in s.") + th

    W <- 14; Hh <- 7.5
    ggsave(file.path(OUTDIR, paste0("s_estimates_", sfx, ".png")), p1, width = W, height = Hh, dpi = 150)
    ggsave(file.path(OUTDIR, paste0("s_bias_",      sfx, ".png")), p2, width = W, height = Hh, dpi = 150)
    ggsave(file.path(OUTDIR, paste0("s_power_",     sfx, ".png")), p3, width = W, height = Hh, dpi = 150)
    ggsave(file.path(OUTDIR, paste0("h2m_",         sfx, ".png")), p4, width = W, height = Hh, dpi = 150)
    ggsave(file.path(OUTDIR, paste0("h2m_bias_",    sfx, ".png")), p5, width = W, height = Hh, dpi = 150)
  }

  ## ---- headline numbers ----------------------------------------------------
  cat("\n-- false-positive rate at s = 0 (should be ~0.05) --\n")
  print(summ[S == 0, .(K, t_max, H2, power, power_ci, frac_bound, fail_rate)][order(K, t_max, -H2)])

  cat("\n-- s = 0.01 --\n")
  print(summ[S == 0.01, .(K, t_max, H2, power, bias, sd_s, rel_prec,
                          h2m_med_ratio, h2m_t0_med_ratio, frac_bound)][order(K, t_max, -H2)])

} else stop("usage: v3_sweep.R grid | cell S H2 CV K t_max | collect")
