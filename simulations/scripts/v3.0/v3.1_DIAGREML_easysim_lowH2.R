library(future); library(furrr); library(purrr)
library(ggplot2); library(data.table); library(scales)
library(progressr); handlers(global = TRUE); handlers("cli")
f = glue::glue

source("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/diagGREML.R")

OUTDIR <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v5.0"
DATDIR <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v5.0"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATDIR, recursive = TRUE, showWarnings = FALSE)

## ---- parameters ------------------------------------------------------------
SIGMA_BETA <- 1; K <- 8; t_max <- 200; LAMBDA <- 50
t_mid <- (seq_len(K) - 0.5) * t_max / K

## 0.5 and 1 are reference points showing the estimator is sound;
## 0.01-0.05 is the realistic range (Wang et al., UKB WES: 1.9-3.4% for significant traits)
H2_vals <- c(1, 0.5, 0.1, 0.05, 0.034, 0.02, 0.01)
Svals   <- c(0, 1e-4, 1e-3, 1e-2)
CV_vals <- c(0.15, 0.30, 0.65)        # 0.65 ~ UKB observed (mean SC 17.4, SD 11.1)
N_vals  <- c(5e5)
N_REPS  <- 100

H2_REALISTIC <- c(0.01, 0.02, 0.034, 0.05)

COL_EST   <- "grey25"       # estimates
COL_TRUTH <- "#D7191C"      # truth, used for nothing else
COL_SD    <- "#2C7BB6"      # per-replicate SD band

h2_lab <- function(x) {
  v <- as.numeric(as.character(x))
  ifelse(v %in% H2_REALISTIC, paste0("h\u00b2 = ", v), paste0("h\u00b2 = ", v, " (ref)"))
}

## ---- simulate one dataset --------------------------------------------------
## seed excludes H2: same genetics across heritability levels, only residual differs
sim_diag <- function(S, REP, H2, CV, N_inds) {
  set.seed(REP + 1e3 * which(Svals == S) + round(CV * 1e4) + N_inds)
  shape <- 1 / CV^2
  C <- sapply(seq_len(K), function(k) rpois(N_inds, rgamma(N_inds, shape, rate = shape / LAMBDA)))
  colnames(C) <- paste0("t", round(t_mid))
  s2 <- SIGMA_BETA^2 * exp(-S * t_mid)
  g  <- rowSums(sqrt(sweep(C, 2, s2, `*`)) * matrix(rnorm(N_inds * K), N_inds, K))
  Y  <- g + rnorm(N_inds, 0, sqrt(var(g) * (1 - H2) / H2))
  list(C = C, Y = Y)
}

## ---- fit one replicate -----------------------------------------------------
run_one <- function(S, REP, H2, CV, N_inds) {
  setDTthreads(1)
  d <- sim_diag(S, REP, H2, CV, N_inds)
  C <- d$C; m <- colMeans(C)
  if (any(!is.finite(m) | m <= 0)) stop("bad mean count")

  A   <- setNames(lapply(seq_len(K), function(k) C[, k] / m[k]), colnames(C))
  fit <- suppressMessages(fit_diagGREML(y = d$Y, X = matrix(1, nrow(C), 1), A = A,
                                        constraint = FALSE, magic0316 = TRUE))
  # message(f("Completed diagGREML for S={S}, H2={H2}, CV={CV}, N={N_inds}, Rep {REP}"))
  idx <- match(names(A), fit$Vlistnames)
  if (anyNA(idx)) stop("component name mismatch")

  bins <- data.table(S = S, REP = REP, H2 = H2, CV = CV, N_inds = N_inds,
                     bin = seq_len(K), mid = t_mid, m = m,
                     s2_true = SIGMA_BETA^2 * exp(-S * t_mid),
                     est = as.numeric(fit$varcmp)[idx],
                     se  = sqrt(diag(fit$Hi))[idx])
  bins[, `:=`(sigma2_b = est / m, se_b = se / m)]

  pos <- bins[is.finite(sigma2_b) & sigma2_b > 0]
  mm  <- if (nrow(pos) >= 3) tryCatch(lm(log(sigma2_b) ~ mid, data = pos),
                                      error = function(e) NULL) else NULL

  out <- data.table(S = S, REP = REP, H2 = H2, CV = CV, N_inds = N_inds,
                    n_pos = nrow(pos), n_neg = K - nrow(pos),
                    sigma2_m = NA_real_, sigma2_m_se_log = NA_real_,
                    s_hat = NA_real_, s_se = NA_real_)
  if (!is.null(mm)) {
    cf <- summary(mm)$coefficients
    out[, `:=`(sigma2_m        = exp(cf[1, 1]),
               sigma2_m_se_log = cf[1, 2],       # SE on the log scale
               s_hat           = -cf[2, 1],
               s_se            = cf[2, 2])]
  }
  list(rep = out, bins = bins)
}

## ---- run the grid for one N ------------------------------------------------
run_for_N <- function(N_inds, n_workers = NULL) {

  if (is.null(n_workers))
    n_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = parallel::detectCores()))
  if (!is.finite(n_workers) || n_workers < 1) n_workers <- 1L

  params <- CJ(S = Svals, REP = seq_len(N_REPS), H2 = H2_vals, CV = CV_vals)[, N_inds := N_inds]
  message(sprintf("N = %s: %d fits on %d workers", format(N_inds, big.mark = ","),
                  nrow(params), n_workers))

  plan(multicore, workers = as.integer(n_workers))
  plan(multicore, workers = as.integer(n_workers))

  res <- with_progress({
    p <- progressor(steps = nrow(params))
    future_pmap(
      as.list(params),
      function(S, REP, H2, CV, N_inds) {
        out <- safely(run_one)(S, REP, H2, CV, N_inds)
        p()
        out
      },
      .options = furrr_options(seed = TRUE, scheduling = Inf,
                               packages = "data.table")
    )
  })

  plan(sequential)
  errs <- keep(res, ~ !is.null(.x$error))
  if (length(errs)) warning(length(errs), "/", length(res), " hard failures at N = ", N_inds,
                            ". First: ", conditionMessage(errs[[1]]$error))

  ok   <- compact(map(res, "result"))
  reps <- rbindlist(map(ok, "rep"),  use.names = TRUE)
  bins <- rbindlist(map(ok, "bins"), use.names = TRUE)

  ## ---- raw output, so every figure can be redrawn without rerunning -------
  tag <- paste0("N", N_inds)
  fwrite(reps, file.path(DATDIR, paste0("reps_", tag, ".csv")))
  fwrite(bins, file.path(DATDIR, paste0("bins_", tag, ".csv")))

  ## ---- per-cell summary ---------------------------------------------------
  reps[, `:=`(ci_lo = s_hat - 1.96 * s_se, ci_hi = s_hat + 1.96 * s_se)]
  reps[, detected := is.finite(s_hat) & is.finite(s_se) & ci_lo > 0]

  summ <- reps[, {
    sh <- s_hat[is.finite(s_hat)]; z <- sh - S
    sm <- sigma2_m[is.finite(sigma2_m)]
    .(n_total   = .N, n_fit = length(sh), fail_rate = 1 - length(sh) / .N,
      mean_n_neg = mean(n_neg, na.rm = TRUE),
      mean_s     = if (length(sh)) mean(sh) else NA_real_,
      median_s   = if (length(sh)) median(sh) else NA_real_,
      bias       = if (length(z)) mean(z) else NA_real_,
      se_bias    = if (length(z) >= 2) sd(z) / sqrt(length(z)) else NA_real_,
      sd_s       = if (length(sh) >= 2) sd(sh) else NA_real_,
      rmse       = if (length(z)) sqrt(mean(z^2)) else NA_real_,
      medae      = if (length(z)) median(abs(z)) else NA_real_,
      rel_prec   = if (length(sh) >= 2 && S > 0) sd(sh) / S else NA_real_,
      power      = mean(detected, na.rm = TRUE),   # at S = 0 this is the false-positive rate
      ## sigma2_m: geometric mean and fold-range, since it's symmetric on the log scale
      gm_sigma2_m  = if (length(sm)) exp(mean(log(sm))) else NA_real_,
      fold_lo      = if (length(sm)) quantile(sm, 0.025) else NA_real_,
      fold_hi      = if (length(sm)) quantile(sm, 0.975) else NA_real_)
  }, by = .(S, H2, CV, N_inds)]
  setorder(summ, CV, -H2, S)
  fwrite(summ, file.path(DATDIR, paste0("summary_", tag, ".csv")))

  ## ---- plotting helpers ---------------------------------------------------
  for (d in list(reps, summ)) {
    d[, H2_f := factor(H2, levels = H2_vals, labels = h2_lab(H2_vals))]
    d[, CV_f := factor(CV, levels = CV_vals, labels = paste0("CV = ", CV_vals))]
    d[, S_f  := factor(S,  levels = Svals)]
  }
  ttl <- paste0("N = ", format(N_inds, big.mark = ",", scientific = FALSE))
  sub_real <- "Panels marked (ref) are above the realistic range for singleton h\u00b2"
  gr <- facet_grid(CV_f ~ H2_f)

  ## (i) estimates of s across replicates
  p1 <- ggplot(reps[is.finite(s_hat)], aes(S_f, s_hat)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_boxplot(outlier.shape = NA, colour = COL_EST, fill = NA, width = 0.6) +
    geom_jitter(width = 0.12, size = 0.3, alpha = 0.25, colour = COL_EST) +
    geom_point(data = unique(reps[, .(S_f, S, CV_f, H2_f)]), aes(S_f, S),
               colour = COL_TRUTH, size = 2.5) +
    gr + coord_cartesian(ylim = c(-0.02, 0.03)) +
    labs(x = "True s", y = expression(hat(s)), title = paste0("Estimates of s, ", ttl),
         subtitle = paste0(sub_real, "; view clipped, box statistics use all data")) +
    theme_bw() + theme(panel.grid.minor = element_blank())

  ## (ii) bias in s
  p2 <- ggplot(summ, aes(S_f, bias)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_linerange(aes(ymin = bias - sd_s, ymax = bias + sd_s),
                   colour = COL_SD, alpha = 0.35, linewidth = 2.5) +
    geom_pointrange(aes(ymin = bias - 1.96 * se_bias, ymax = bias + 1.96 * se_bias),
                    colour = COL_EST, size = 0.35) +
    gr + coord_cartesian(ylim = c(-0.015, 0.015)) +
    labs(x = "True s", y = expression(hat(s) - s), title = paste0("Bias in \u015d, ", ttl),
         subtitle = "Dark: 95% CI on the mean. Pale blue: SD across replicates (what one study faces).") +
    theme_bw() + theme(panel.grid.minor = element_blank())

  ## (iii) power / false-positive rate
  p3 <- ggplot(summ, aes(S_f, power, group = 1)) +
    geom_hline(yintercept = 0.05, colour = "grey60", linetype = "dotted") +
    geom_hline(yintercept = 0.80, colour = COL_TRUTH, linetype = "dashed") +
    geom_line(colour = COL_EST, alpha = 0.6) +
    geom_point(colour = COL_EST, size = 2) +
    gr + scale_y_continuous(limits = c(0, 1)) +
    labs(x = "True s", y = "P(95% CI excludes 0)",
         title = paste0("Power to detect s > 0, ", ttl),
         subtitle = "At s = 0 this is the false-positive rate (dotted 0.05); dashed red = 80% power") +
    theme_bw() + theme(panel.grid.minor = element_blank())

  ## (iv) sigma2_m recovery, log scale (truth = SIGMA_BETA^2)
  p4 <- ggplot(reps[is.finite(sigma2_m) & sigma2_m > 0], aes(S_f, sigma2_m)) +
    geom_hline(yintercept = SIGMA_BETA^2, colour = COL_TRUTH, linetype = "dashed", linewidth = 0.8) +
    geom_boxplot(outlier.shape = NA, colour = COL_EST, fill = NA, width = 0.6) +
    geom_jitter(width = 0.12, size = 0.3, alpha = 0.25, colour = COL_EST) +
    gr + scale_y_log10(breaks = c(0.1, 0.3, 1, 3, 10), labels = label_number(accuracy = 0.1)) +
    coord_cartesian(ylim = c(0.1, 10)) +
    labs(x = "True s", y = expression(hat(sigma)[m]^2),
         title = paste0("Estimates of ", "\u03c3\u00b2", "\u2098, ", ttl),
         subtitle = paste0("Truth = ", SIGMA_BETA^2,
                           " (dashed red); log axis, view clipped to [0.1, 10]")) +
    theme_bw() + theme(panel.grid.minor = element_blank())

  W <- 14; Hh <- 7
  ggsave(file.path(OUTDIR, paste0("s_estimates_", tag, ".png")), p1, width = W, height = Hh, dpi = 150)
  ggsave(file.path(OUTDIR, paste0("s_bias_",      tag, ".png")), p2, width = W, height = Hh, dpi = 150)
  ggsave(file.path(OUTDIR, paste0("s_power_",     tag, ".png")), p3, width = W, height = Hh, dpi = 150)
  ggsave(file.path(OUTDIR, paste0("sigma2m_",     tag, ".png")), p4, width = W, height = Hh, dpi = 150)

  list(reps = reps, bins = bins, summary = summ, errors = errs,
       p1 = p1, p2 = p2, p3 = p3, p4 = p4)
}

## ---- run -------------------------------------------------------------------
all_results <- map(N_vals, run_for_N)
names(all_results) <- paste0("N", N_vals)

## combined files across all N, for cross-N figures later
fwrite(rbindlist(map(all_results, "reps")),    file.path(DATDIR, "reps_allN.csv"))
fwrite(rbindlist(map(all_results, "summary")), file.path(DATDIR, "summary_allN.csv"))

## ---- the headline numbers --------------------------------------------------
all_summ <- rbindlist(map(all_results, "summary"))
all_summ[S == 0.01, .(N_inds, CV, H2, power, bias, sd_s, rel_prec)][order(N_inds, CV, -H2)]
