#!/usr/bin/env Rscript
## v5.1 sweep: two-stage (diagGREML + lm) vs joint REML, one SLURM array task per cell
##   Rscript v5_sweep.R grid                    write cells.txt, one line per cell
##   Rscript v5_sweep.R cell S H2 CV K t_max    run all reps for one cell
##   Rscript v5_sweep.R collect                 combine cells, summarise, plot

library(data.table)
RNGkind("L'Ecuyer-CMRG")   # what furrr's seed = TRUE was using underneath set.seed()

VERSION <- "v5.1"
BASE <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations"
source(file.path(BASE, "programs/diagGREML.R"))

OUTDIR   <- file.path(BASE, "figs", VERSION)
DATDIR   <- file.path(BASE, "data", VERSION)
CELLDIR  <- file.path(DATDIR, "cells")
CELLFILE <- file.path(DATDIR, "cells.txt")
for (d in c(OUTDIR, CELLDIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- parameters ------------------------------------------------------------
SIGMA_BETA <- 1
LAMBDA     <- 1000                 # singletons per person per bin
                                   # (UKB 150k WGS British/Irish: ~1330 per person total)

H2_vals <- c(0.5, 0.1, 0.034, 0.02, 0.01)   # 0.5 = reference; 0.01-0.034 realistic (Wang et al.)
Svals   <- c(0, 1e-4, 1e-3, 1e-2, 3e-2)
CV_vals <- c(0.15, 0.30, 0.65)              # 0.65 ~ UKB (mean SC 17.4, SD 11.1); 0.37 if outliers removed
N_INDS  <- 5e5
K_VALS  <- c(3, 8)
T_VALS  <- c(200, 400)
N_REPS  <- 100

H2_REALISTIC <- c(0.01, 0.02, 0.034, 0.05)

METHODS   <- c("two-stage", "joint")
COL_M     <- c("two-stage" = "grey25", "joint" = "#2C7BB6")
COL_TRUTH <- "#D7191C"

h2_lab <- function(x) {
  v <- as.numeric(as.character(x))
  ifelse(v %in% H2_REALISTIC, paste0("h\u00b2 = ", v), paste0("h\u00b2 = ", v, " (ref)"))
}

cell_grid <- function() CJ(S = Svals, H2 = H2_vals, CV = CV_vals, K = K_VALS, t_max = T_VALS)
cell_tag  <- function(S, H2, CV, K, t_max) sprintf("S%g_H%g_CV%g_K%g_t%g", S, H2, CV, K, t_max)
write_atomic <- function(x, f) { tmp <- paste0(f, ".tmp"); fwrite(x, tmp); file.rename(tmp, f) }

## ---- joint REML ------------------------------------------------------------
## y_i ~ N(Z_i b, v_i),  v_i = s2m * sum_k x_ik exp(-s t_k) + s2e
## C: N x K matrix of singleton counts; t: bin ages; Z: fixed-effect design matrix
fit_sel_reml <- function(y, C, t, Z = matrix(1, length(y), 1), reml = TRUE) {
  storage.mode(C) <- "double"                  # avoid re-coercing integer counts at every evaluation

  var_i <- function(par)                       # par = (log s2m, s, log s2e)
    exp(par[1]) * as.vector(C %*% exp(-par[2] * t)) + exp(par[3])

  gls <- function(v) {                         # fixed effects given the variances
    ZtViZ <- crossprod(Z, Z / v)
    list(ZtViZ = ZtViZ, b = solve(ZtViZ, crossprod(Z, y / v)))
  }

  nll <- function(par) {
    v <- var_i(par)
    if (!all(is.finite(v)) || any(v <= 0)) return(Inf)   # reject overflowing trial steps
    g <- gls(v)
    r <- y - as.vector(Z %*% g$b)
    out <- sum(log(v)) + sum(r^2 / v)
    if (reml) out <- out + as.numeric(determinant(g$ZtViZ)$modulus)
    0.5 * out
  }

  start <- c(log(0.02 * var(y) / mean(rowSums(C))), 0, log(0.98 * var(y)))
  opt <- optim(start, nll, method = "BFGS", hessian = TRUE,
               control = list(parscale = c(1, 1e-2, 1), reltol = 1e-12, maxit = 500))

  se <- tryCatch(suppressWarnings(sqrt(diag(solve(opt$hessian)))),
                 error = function(e) rep(NA_real_, 3))
  list(s = opt$par[2], s_se = se[2],
       sigma2_m = exp(opt$par[1]), log_sigma2_m_se = se[1],
       sigma2_e = exp(opt$par[3]),
       b = as.vector(gls(var_i(opt$par))$b),
       logLik = -opt$value, convergence = opt$convergence)
}

## ---- simulate one dataset --------------------------------------------------
sim_diag <- function(S, REP, H2, CV, N_inds, K, t_max) {
  set.seed(REP + 1e3 * which(Svals == S) + round(CV * 1e4) + K * 1e5 + t_max * 1e6)
  t_mid <- (seq_len(K) - 0.5) * t_max / K
  shape <- 1 / CV^2
  C <- sapply(seq_len(K), function(k) rpois(N_inds, rgamma(N_inds, shape, rate = shape / LAMBDA)))
  colnames(C) <- paste0("t", round(t_mid))
  s2 <- SIGMA_BETA^2 * exp(-S * t_mid)
  g  <- rowSums(sqrt(sweep(C, 2, s2, `*`)) * matrix(rnorm(N_inds * K), N_inds, K))
  Y  <- g + rnorm(N_inds, 0, sqrt(var(g) * (1 - H2) / H2))
  list(C = C, Y = Y, t_mid = t_mid, s2_true = s2)
}

## ---- fit one replicate with both methods -----------------------------------
run_one <- function(S, REP, H2, CV, N_inds, K, t_max) {
  setDTthreads(1)
  d <- sim_diag(S, REP, H2, CV, N_inds, K, t_max)
  C <- d$C; m <- colMeans(C); t_mid <- d$t_mid
  if (any(!is.finite(m) | m <= 0)) stop("bad mean count")

  ## mutational heritability: V_M = (new mutations per person per generation) * sigma2_m
  ## bin width is t_max/K, so mutations per person per generation = LAMBDA * K / t_max
  mu_per_gen <- LAMBDA * K / t_max
  V_P        <- var(d$Y)
  h2m_true   <- mu_per_gen * SIGMA_BETA^2 / V_P

  ## (a) two-stage: diagGREML per bin, then lm(log sigma2_b ~ age)
  A   <- setNames(lapply(seq_len(K), function(k) C[, k] / m[k]), colnames(C))
  fit <- suppressMessages(fit_diagGREML(y = d$Y, X = matrix(1, nrow(C), 1), A = A,
                                        constraint = FALSE, magic0316 = TRUE))
  idx <- match(names(A), fit$Vlistnames)
  if (anyNA(idx)) stop("component name mismatch")

  bins <- data.table(S = S, REP = REP, H2 = H2, CV = CV, N_inds = N_inds,
                     K = K, t_max = t_max,
                     bin = seq_len(K), mid = t_mid, m = m, s2_true = d$s2_true,
                     est = as.numeric(fit$varcmp)[idx],
                     se  = sqrt(diag(fit$Hi))[idx])
  bins[, `:=`(sigma2_b = est / m, se_b = se / m)]

  pos <- bins[is.finite(sigma2_b) & sigma2_b > 0]
  mm  <- if (nrow(pos) >= 3) tryCatch(lm(log(sigma2_b) ~ mid, data = pos),
                                      error = function(e) NULL) else NULL

  ## (b) joint: one-stage REML for sigma2_m and s
  jf <- tryCatch(fit_sel_reml(y = d$Y, C = C, t = t_mid), error = function(e) NULL)

  out <- data.table(S = S, REP = REP, H2 = H2, CV = CV, N_inds = N_inds,
                    K = K, t_max = t_max,
                    n_pos = nrow(pos), n_neg = K - nrow(pos),
                    V_P = V_P, mu_per_gen = mu_per_gen, h2m_true = h2m_true,
                    sigma2_m_sep = NA_real_, sigma2_m_se_log_sep = NA_real_,
                    h2m_hat_sep = NA_real_, s_hat_sep = NA_real_, s_se_sep = NA_real_,
                    sigma2_m_joint = NA_real_, sigma2_m_se_log_joint = NA_real_,
                    h2m_hat_joint = NA_real_, s_hat_joint = NA_real_, s_se_joint = NA_real_,
                    sigma2_e_joint = NA_real_, mu_joint = NA_real_,
                    loglik_joint = NA_real_, conv_joint = NA_integer_)

  if (!is.null(mm)) {
    cf <- summary(mm)$coefficients
    out[, `:=`(sigma2_m_sep        = exp(cf[1, 1]),
               sigma2_m_se_log_sep = cf[1, 2],
               h2m_hat_sep         = mu_per_gen * exp(cf[1, 1]) / V_P,
               s_hat_sep           = -cf[2, 1],
               s_se_sep            = cf[2, 2])]
  }
  if (!is.null(jf)) {
    out[, `:=`(sigma2_m_joint        = jf$sigma2_m,
               sigma2_m_se_log_joint = jf$log_sigma2_m_se,
               h2m_hat_joint         = mu_per_gen * jf$sigma2_m / V_P,
               s_hat_joint           = jf$s,
               s_se_joint            = jf$s_se,
               sigma2_e_joint        = jf$sigma2_e,
               mu_joint              = jf$b[1],
               loglik_joint          = jf$logLik,
               conv_joint            = as.integer(jf$convergence))]
  }
  list(rep = out, bins = bins)
}

## ---- dispatch --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else ""

if (mode == "grid") {

  fwrite(cell_grid(), CELLFILE, sep = " ", col.names = FALSE)
  message(nrow(cell_grid()), " cells written to ", CELLFILE)

} else if (mode == "cell") {

  a <- suppressWarnings(as.numeric(args[-1]))
  if (length(a) != 5 || anyNA(a)) stop("usage: v5_sweep.R cell S H2 CV K t_max")
  S <- Svals[which.min(abs(Svals - a[1]))]        # snap to grid: the seed uses which(Svals == S)
  if (abs(S - a[1]) > 1e-12) stop("S must be one of Svals")
  H2 <- a[2]; CV <- a[3]; K <- as.integer(a[4]); t_max <- a[5]

  tag   <- cell_tag(S, H2, CV, K, t_max)
  f_rep <- file.path(CELLDIR, paste0("reps_", tag, ".csv"))
  f_bin <- file.path(CELLDIR, paste0("bins_", tag, ".csv"))
  if (file.exists(f_rep)) { message(tag, ": already done, skipping"); quit(save = "no") }

  t0  <- Sys.time()
  res <- lapply(seq_len(N_REPS), function(REP)
    tryCatch(run_one(S, REP, H2, CV, N_INDS, K, t_max),
             error = function(e) { message("rep ", REP, " failed: ", conditionMessage(e)); NULL }))
  ok <- Filter(Negate(is.null), res)
  if (!length(ok)) stop("all reps failed")

  write_atomic(rbindlist(lapply(ok, `[[`, "bins"), use.names = TRUE), f_bin)
  write_atomic(rbindlist(lapply(ok, `[[`, "rep"),  use.names = TRUE), f_rep)  # written last = "done" marker
  message(sprintf("%s: %d/%d reps in %.1f min", tag, length(ok), N_REPS,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))

} else if (mode == "collect") {

  library(ggplot2); library(scales)

  cells <- cell_grid()
  tags  <- cells[, cell_tag(S, H2, CV, K, t_max)]
  f_rep <- file.path(CELLDIR, paste0("reps_", tags, ".csv"))
  f_bin <- file.path(CELLDIR, paste0("bins_", tags, ".csv"))
  have  <- file.exists(f_rep)
  if (!any(have)) stop("no finished cells in ", CELLDIR)
  if (!all(have)) warning(sum(!have), " cells missing. Rerun with: sbatch --array=",
                          paste(which(!have), collapse = ","), " v5_submit.sh")

  reps <- rbindlist(lapply(f_rep[have], fread), use.names = TRUE)
  bins <- rbindlist(lapply(f_bin[have], fread), use.names = TRUE)
  short <- reps[, .N, by = .(S, H2, CV, K, t_max)][N < N_REPS]
  if (nrow(short)) { message("Cells with failed reps:"); print(short) }

  tag <- paste0("N", N_INDS)
  fwrite(reps, file.path(DATDIR, paste0("reps_", tag, ".csv")))
  fwrite(bins, file.path(DATDIR, paste0("bins_", tag, ".csv")))

  ## ---- long format: one row per replicate per method -----------------------
  ## joint fits count as failed unless optim converged and gave a finite SE
  ok_j   <- reps[, conv_joint %in% 0 & is.finite(s_se_joint)]
  common <- reps[, .(S, REP, H2, CV, N_inds, K, t_max, n_neg, h2m_true)]
  long <- rbind(
    cbind(common, reps[, .(method = "two-stage",
                           s_hat = s_hat_sep, s_se = s_se_sep, h2m_hat = h2m_hat_sep)]),
    cbind(common, reps[, .(method = "joint",
                           s_hat   = ifelse(ok_j, s_hat_joint,   NA_real_),
                           s_se    = ifelse(ok_j, s_se_joint,    NA_real_),
                           h2m_hat = ifelse(ok_j, h2m_hat_joint, NA_real_))]))
  long[, method := factor(method, levels = METHODS)]
  long[, ci_lo := s_hat - 1.96 * s_se]
  long[, detected := is.finite(s_hat) & is.finite(s_se) & ci_lo > 0]
  fwrite(long, file.path(DATDIR, paste0("reps_long_", tag, ".csv")))

  ## ---- per-cell, per-method summary ----------------------------------------
  summ <- long[, {
    sh <- s_hat[is.finite(s_hat)]; z <- sh - S
    hm <- h2m_hat[is.finite(h2m_hat) & h2m_hat > 0]
    .(n_total = .N, n_fit = length(sh), fail_rate = 1 - length(sh) / .N,
      mean_n_neg = mean(n_neg, na.rm = TRUE),
      h2m_true   = h2m_true[1],
      mean_s = if (length(sh)) mean(sh) else NA_real_,
      median_s = if (length(sh)) median(sh) else NA_real_,
      bias    = if (length(z)) mean(z) else NA_real_,
      se_bias = if (length(z) >= 2) sd(z) / sqrt(length(z)) else NA_real_,
      sd_s    = if (length(sh) >= 2) sd(sh) else NA_real_,
      rmse    = if (length(z)) sqrt(mean(z^2)) else NA_real_,
      medae   = if (length(z)) median(abs(z)) else NA_real_,
      rel_prec = if (length(sh) >= 2 && S > 0) sd(sh) / S else NA_real_,
      power   = mean(detected, na.rm = TRUE),
      gm_h2m   = if (length(hm)) exp(mean(log(hm))) else NA_real_,
      gmed_h2m = if (length(hm)) median(hm) else NA_real_,
      fold_lo  = if (length(hm)) quantile(hm, 0.025) else NA_real_,
      fold_hi  = if (length(hm)) quantile(hm, 0.975) else NA_real_)
  }, by = .(method, S, H2, CV, N_inds, K, t_max)]
  setorder(summ, method, K, t_max, CV, -H2, S)
  fwrite(summ, file.path(DATDIR, paste0("summary_", tag, ".csv")))

  ## ---- figures: one set per (K, t_max), facet CV x H2, colour = method ------
  for (d in list(long, summ)) {
    d[, H2_f := factor(H2, levels = H2_vals, labels = h2_lab(H2_vals))]
    d[, CV_f := factor(CV, levels = CV_vals, labels = paste0("CV = ", CV_vals))]
    d[, S_f  := factor(S,  levels = Svals)]
  }

  col_scale <- scale_colour_manual(values = COL_M, name = NULL)
  th  <- theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  dg  <- position_dodge(width = 0.75)
  jd  <- position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75)
  dg2 <- position_dodge(width = 0.5)

  for (kk in K_VALS) for (tt in T_VALS) {
    r <- long[K == kk & t_max == tt]; s <- summ[K == kk & t_max == tt]
    if (!nrow(r)) next
    sfx <- paste0(tag, "_K", kk, "_t", tt)
    ttl <- sprintf("N = %s, %d bins, ages 0-%d gen",
                   format(N_INDS, big.mark = ",", scientific = FALSE), kk, tt)

    ## (i) estimates of s - pseudo-log y axis (handles 0 and negatives)
    p1 <- ggplot(r[is.finite(s_hat)], aes(S_f, s_hat, colour = method)) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_boxplot(outlier.shape = NA, fill = NA, width = 0.6, position = dg) +
      geom_point(position = jd, size = 0.3, alpha = 0.25) +
      geom_point(data = unique(r[, .(S_f, S, CV_f, H2_f)]), aes(S_f, S),
                 colour = COL_TRUTH, size = 2.5, inherit.aes = FALSE) +
      facet_grid(CV_f ~ H2_f) + col_scale +
      scale_y_continuous(trans = pseudo_log_trans(sigma = 1e-4),
                         breaks = c(-1e-2, -1e-3, 0, 1e-4, 1e-3, 1e-2, 3e-2, 1e-1)) +
      labs(x = "True s", y = expression(hat(s)), title = paste0("Estimates of s: ", ttl),
           subtitle = "Pseudo-log y axis (linear near zero, log beyond); red = truth") + th

    ## (ii) bias in s
    p2 <- ggplot(s, aes(S_f, bias, colour = method)) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_linerange(aes(ymin = bias - sd_s, ymax = bias + sd_s),
                     alpha = 0.3, linewidth = 2.5, position = dg2) +
      geom_pointrange(aes(ymin = bias - 1.96 * se_bias, ymax = bias + 1.96 * se_bias),
                      size = 0.35, position = dg2) +
      facet_grid(CV_f ~ H2_f, scales = "free_y") + col_scale +
      labs(x = "True s", y = expression(hat(s) - s), title = paste0("Bias in \u015d: ", ttl),
           subtitle = "Dark: 95% CI on the mean. Pale: SD across replicates (what one study faces).") + th

    ## (iii) power / false-positive rate
    p3 <- ggplot(s, aes(S_f, power, colour = method, group = method)) +
      geom_hline(yintercept = 0.05, colour = "grey60", linetype = "dotted") +
      geom_hline(yintercept = 0.80, colour = COL_TRUTH, linetype = "dashed") +
      geom_line(alpha = 0.6) + geom_point(size = 2) +
      facet_grid(CV_f ~ H2_f) + col_scale + scale_y_continuous(limits = c(0, 1)) +
      labs(x = "True s", y = "P(95% CI excludes 0)", title = paste0("Power to detect s > 0: ", ttl),
           subtitle = "At s = 0 this is the false-positive rate (dotted 0.05); dashed red = 80%") + th

    ## (iv) mutational heritability V_M / V_P
    tru <- unique(r[, .(S_f, CV_f, H2_f, h2m_true)])
    p4 <- ggplot(r[is.finite(h2m_hat) & h2m_hat > 0], aes(S_f, h2m_hat, colour = method)) +
      geom_boxplot(outlier.shape = NA, fill = NA, width = 0.6, position = dg) +
      geom_point(position = jd, size = 0.3, alpha = 0.25) +
      geom_point(data = tru, aes(S_f, h2m_true), colour = COL_TRUTH, size = 2.5,
                 inherit.aes = FALSE) +
      facet_grid(CV_f ~ H2_f) + col_scale +
      scale_y_log10(labels = label_number(accuracy = 0.0001)) +
      labs(x = "True s", y = expression(hat(V)[M] / V[P]),
           title = paste0("Mutational heritability: ", ttl),
           subtitle = "Red = truth; log axis. Literature V_M/V_P is typically ~1e-3.") + th

    W <- 14; Hh <- 7.5
    ggsave(file.path(OUTDIR, paste0("s_estimates_", sfx, ".png")), p1, width = W, height = Hh, dpi = 150)
    ggsave(file.path(OUTDIR, paste0("s_bias_",      sfx, ".png")), p2, width = W, height = Hh, dpi = 150)
    ggsave(file.path(OUTDIR, paste0("s_power_",     sfx, ".png")), p3, width = W, height = Hh, dpi = 150)
    ggsave(file.path(OUTDIR, paste0("h2m_",         sfx, ".png")), p4, width = W, height = Hh, dpi = 150)
  }

  ## ---- headline numbers ---------------------------------------------------
  print(summ[S == 0.01, .(method, K, t_max, CV, H2, power, bias, sd_s, rel_prec, fail_rate)][
        order(K, t_max, CV, -H2, method)])

} else stop("usage: v5_sweep.R grid | cell S H2 CV K t_max | collect")
