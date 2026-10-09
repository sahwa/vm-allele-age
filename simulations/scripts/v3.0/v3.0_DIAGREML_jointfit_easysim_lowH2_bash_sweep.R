#!/usr/bin/env Rscript
## v5.1 sweep: two-stage (diagGREML + lm) vs joint REML, one SLURM array task per cell
##   Rscript v5_sweep.R grid                    write cells.txt, one line per cell
##   Rscript v5_sweep.R cell S H2 CV K t_max    run all reps for one cell
##   Rscript v5_sweep.R collect                 combine cells, summarise, plot

library(data.table)
RNGkind("L'Ecuyer-CMRG")   # what furrr's seed = TRUE was using underneath set.seed()

VERSION <- "v3.2"
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
## Joint REML by profile likelihood over s.
## y_i ~ N(Z_i b, v_i),  v_i = a * w_i(s) + s2e,
## w_i(s) = sum_k x_ik exp(-s (t_k - t0)), scaled to mean 1; t0 = count-weighted mean age.
## a can be negative. sigma2_m is the per-mutation variance extrapolated to age 0.
fit_sel_reml <- function(y, C, t, Z = matrix(1, length(y), 1), s_max = 0.1) {
  storage.mode(C) <- "double"
  cm <- colMeans(C)
  t0 <- sum(cm * t) / sum(cm)
  tc <- t - t0
  vy <- var(y)
  theta <- c(0.05 * vy, 0.95 * vy)              # (a, s2e); warm-started across values of s

  ## REML fit of (a, s2e) for one fixed s, by Fisher scoring with step-halving
  inner <- function(s, th) {
    w <- as.vector(C %*% exp(-s * tc)); wbar <- mean(w); w <- w / wbar
    obj <- function(th) {
      v <- th[1] * w + th[2]
      if (th[2] <= 0 || min(v) <= 0) return(NULL)
      ZtViZ <- crossprod(Z, Z / v)
      b <- solve(ZtViZ, crossprod(Z, y / v))
      r <- y - as.vector(Z %*% b)
      list(v = v, r = r, b = as.vector(b), ZtViZ = ZtViZ,
           nll = 0.5 * (sum(log(v)) + sum(r^2 / v) + as.numeric(determinant(ZtViZ)$modulus)))
    }
    cur <- obj(th)
    if (is.null(cur)) { th <- c(0, vy); cur <- obj(th) }
    converged <- FALSE
    for (it in 1:50) {
      v <- cur$v
      h <- rowSums((Z %*% solve(cur$ZtViZ)) * Z)            # REML leverage term
      q <- cur$r^2 / v^2 - (1 - h / v) / v
      score <- 0.5 * c(sum(q * w), sum(q))
      info  <- 0.5 * matrix(c(sum(w^2 / v^2), sum(w / v^2),
                              sum(w / v^2),   sum(1 / v^2)), 2)
      step <- solve(info, score)
      new <- NULL
      for (half in 1:25) {
        new <- obj(th + step)
        if (!is.null(new) && new$nll <= cur$nll + 1e-8) break
        new <- NULL; step <- step / 2
      }
      if (is.null(new)) break
      th <- th + step; cur <- new
      if (max(abs(step)) < 1e-8 * vy) { converged <- TRUE; break }
    }
    list(nll = cur$nll, th = th, b = cur$b, wbar = wbar, converged = converged)
  }

  prof <- function(s) { f <- inner(s, theta); theta <<- f$th; f$nll }

  ## coarse search on a grid that is dense near zero, then refine around the best point
  g    <- s_max * 10^seq(-3, 0, length.out = 13)
  grid <- c(-rev(g), 0, g)
  pg   <- vapply(grid, prof, 0)
  j    <- which.min(pg)
  s_hat <- optimize(prof, c(grid[max(j - 1, 1)], grid[min(j + 1, length(grid))]),
                    tol = 1e-6)$minimum
  f <- inner(s_hat, theta)
  at_bound <- abs(s_hat) > s_max * (1 - 1e-3)

  ## SE of s from the curvature of the profile; likelihood-ratio test of s = 0
  h    <- 2e-4
  curv <- (prof(s_hat + h) - 2 * f$nll + prof(s_hat - h)) / h^2
  s_se <- if (!at_bound && is.finite(curv) && curv > 0) 1 / sqrt(curv) else NA_real_
  lrt  <- max(0, 2 * (pg[grid == 0] - f$nll))

  a_t0 <- f$th[1] / f$wbar                       # per-mutation variance at age t0
  list(s = s_hat, s_se = s_se,
       sigma2_m = a_t0 * exp(s_hat * t0), log_sigma2_m_se = NA_real_,
       sigma2_e = f$th[2], b = f$b,
       logLik = -f$nll, convergence = if (f$converged) 0L else 1L,
       lrt = lrt, lrt_p = pchisq(lrt, 1, lower.tail = FALSE),
       at_bound = at_bound, sigma2_t0 = a_t0, t0 = t0)
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
    sh  <- s_hat[is.finite(s_hat)]; z <- sh - S
    okh <- is.finite(h2m_hat)
    zh  <- h2m_hat[okh] - h2m_true[okh]          # per-rep error in h2m
    rh  <- h2m_hat[okh] / h2m_true[okh]          # per-rep ratio to the truth
    hm  <- h2m_hat[okh & h2m_hat > 0]
    .(n_total = .N, n_fit = length(sh), fail_rate = 1 - length(sh) / .N,
      mean_n_neg = mean(n_neg, na.rm = TRUE),
      h2m_true   = mean(h2m_true),
      mean_s = if (length(sh)) mean(sh) else NA_real_,
      median_s = if (length(sh)) median(sh) else NA_real_,
      bias    = if (length(z)) mean(z) else NA_real_,
      se_bias = if (length(z) >= 2) sd(z) / sqrt(length(z)) else NA_real_,
      sd_s    = if (length(sh) >= 2) sd(sh) else NA_real_,
      rmse    = if (length(z)) sqrt(mean(z^2)) else NA_real_,
      medae   = if (length(z)) median(abs(z)) else NA_real_,
      rel_prec = if (length(sh) >= 2 && S > 0) sd(sh) / S else NA_real_,
      power   = mean(detected, na.rm = TRUE),
      h2m_bias      = if (length(zh)) mean(zh) else NA_real_,
      h2m_se_bias   = if (length(zh) >= 2) sd(zh) / sqrt(length(zh)) else NA_real_,
      h2m_rel_bias  = if (length(rh)) mean(rh) - 1 else NA_real_,
      h2m_med_ratio = if (length(rh)) median(rh) else NA_real_,
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
  th  <- theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "bottom", axis.text = element_text(size=16), axis.title = element_text(size=16), strip.text = element_text(size=16))
  dg  <- position_dodge(width = 0.75)
  jd  <- position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75)
  dg2 <- position_dodge(width = 0.5)

  for (kk in K_VALS) for (tt in T_VALS) {
  r <- long[K == kk & t_max == tt]; s <- summ[K == kk & t_max == tt]

  CV_select <- 0.65
  S_select  <- c(0.001, 0.01)
  H2_select <- c(0.034, 0.02, 0.01)

  r <- r[CV %in% CV_select & S %in% S_select & H2 %in% H2_select]
  s <- s[CV %in% CV_select & S %in% S_select & H2 %in% H2_select]
    if (!nrow(r)) next
    sfx <- paste0(tag, "_K", kk, "_t", tt)
    ttl <- sprintf("N = %s, %d bins, ages 0-%d gen",
                   format(N_INDS, big.mark = ",", scientific = FALSE), kk, tt)

    ## (i) estimates of s - pseudo-log y axis (handles 0 and negatives)
    tru_s <- unique(r[, .(S_f, S, CV_f, H2_f)])
    ann_tab <- function(r, est, truth) {
    a <- r[is.finite(get(est)),
            .(bias = mean(get(est) - get(truth)),
              rmse = sqrt(mean((get(est) - get(truth))^2)),
              n = .N),
            by = .(method, S_f, CV_f, H2_f)]
      a[, lab := sprintf("bias %s\nRMSE %s",
                        formatC(bias, format = "g", digits = 2, flag = "+"),
                        formatC(rmse, format = "g", digits = 2))]
      a[]
    }

    ann_s <- ann_tab(r, "s_hat",   "S")          # for the s figure
    ann_v <- ann_tab(r, "h2m_hat", "h2m_true")   # for the V_M figure

    p1 <- ggplot(r[is.finite(s_hat)], aes(S_f, s_hat, colour = method)) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_point(position = jd, size = 1, alpha = 0.25) +
      geom_boxplot(outlier.shape = NA, fill = NA, width = 0.6, position = dg) +
      geom_errorbar(data = tru_s, aes(x = S_f, ymin = S, ymax = S), width = 0.9,
                    colour = COL_TRUTH, linewidth = 0.8, inherit.aes = FALSE) +
      facet_grid(CV_f ~ H2_f) + col_scale +
      geom_text(data = ann_s, aes(S_f, Inf, label = lab, colour = method),
          position = position_dodge(width = 0.9), vjust = 1.2,
          size = 2.6, lineheight = 0.9, show.legend = FALSE) +
      scale_y_continuous(trans = pseudo_log_trans(sigma = 2e-3),
                        breaks = c(-0.1, -0.01, -0.001, 0, 0.001, 0.01, 0.03, 0.1),
                        labels = c("-0.1", "-0.01", "-0.001", "0", "0.001", "0.01", "0.03", "0.1")) +
      labs(x = "True s", y = expression(hat(s)), title = paste0("Estimates of s: ", ttl),
          subtitle = "Pseudo-log y axis; red bar = truth; points at \u00b10.1 are joint fits at the search bound") + th

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
    p4 <- ggplot(r[is.finite(h2m_hat)], aes(S_f, h2m_hat, colour = method)) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_boxplot(outlier.shape = NA, fill = NA, width = 0.6, position = dg) +
      geom_point(position = jd, size = 1, alpha = 0.25) +
      geom_point(data = tru, aes(S_f, h2m_true), colour = COL_TRUTH, size = 2.5,
                inherit.aes = FALSE) +
      facet_grid(CV_f ~ H2_f) + col_scale +
      scale_y_continuous(trans = pseudo_log_trans(sigma = 1e-5),
                        breaks = c(-1e-3, -1e-4, 0, 1e-4, 1e-3, 1e-2),
                        labels = function(x) format(x, scientific = TRUE)) +
      coord_cartesian(ylim = c(-3e-3, 3e-2)) +
      geom_text(data = ann_v, aes(S_f, Inf, label = lab, colour = method),
          position = position_dodge(width = 0.9), vjust = 1.2,
          size = 2.6, lineheight = 0.9, show.legend = FALSE) +
      labs(x = "True s", y = expression(hat(V)[M] / V[P]),
          title = paste0("Mutational heritability: ", ttl),
          subtitle = "Red = truth; pseudo-log axis, includes zero and negative estimates; view clipped") + th

    hb <- r[is.finite(h2m_hat) & h2m_hat > 0,
            .(m = mean(log(h2m_hat / h2m_true)),
              se = sd(log(h2m_hat / h2m_true)) / sqrt(.N)),
            by = .(method, S_f, CV_f, H2_f)]

    ## (v) mutational heritability V_M / V_P bias
    p5 <- ggplot(hb, aes(S_f, exp(m), colour = method)) +
      geom_hline(yintercept = 1, colour = COL_TRUTH, linetype = "dashed") +
      geom_pointrange(aes(ymin = exp(m - 1.96 * se), ymax = exp(m + 1.96 * se)),
                      size = 1, position = dg2) +
      facet_grid(CV_f ~ H2_f) + col_scale +
      scale_y_log10() +
      labs(x = "True s", y = expression(hat(V)[M] / V[M] ~ "(estimate / truth)"),
           title = paste0("Bias in mutational heritability: ", ttl),
           subtitle = "Geometric mean ratio with 95% CI; dashed red = unbiased (ratio of 1); log axis") + th



    W <- 14; Hh <- 7.5
    ggsave(file.path(OUTDIR, paste0("s_estimates_", sfx, ".png")), p1, width = W, height = Hh, dpi = 150)
    ggsave(file.path(OUTDIR, paste0("s_bias_",      sfx, ".png")), p2, width = W, height = Hh, dpi = 150)
    ggsave(file.path(OUTDIR, paste0("s_power_",     sfx, ".png")), p3, width = W, height = Hh, dpi = 150)
    ggsave(file.path(OUTDIR, paste0("h2m_",         sfx, ".png")), p4, width = W, height = Hh, dpi = 150)
    ggsave(file.path(OUTDIR, paste0("h2m_bias_",    sfx, ".png")), p5, width = W, height = Hh, dpi = 150)
  }

  ## ---- headline numbers ---------------------------------------------------
  print(summ[S == 0.01, .(method, K, t_max, CV, H2, power, bias, sd_s, rel_prec, fail_rate)][order(K, t_max, CV, -H2, method)])

} else stop("usage: v5_sweep.R grid | cell S H2 CV K t_max | collect")
