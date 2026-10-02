library(future); library(furrr); library(purrr)
library(ggplot2); library(data.table); library(scales)

source("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/diagGREML.R")

OUTDIR <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v3.0"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## ---- parameters ------------------------------------------------------------
SIGMA_BETA <- 1; K <- 8; t_max <- 200; LAMBDA <- 50
t_mid <- (seq_len(K) - 0.5) * t_max / K

H2_vals <- c(0.01, 0.034, 0.1)
Svals   <- c(0, 1e-4, 1e-2)
CV_vals <- seq(0.15, 0.15)
N_vals  <- c(2e5, 3e5, 5e5)
N_REPS  <- 50

## ---- simulate one dataset --------------------------------------------------
sim_diag <- function(S, REP, H2, CV, N_inds) {
  set.seed(REP + 1000 * which(Svals == S) + round(CV * 1e4) + N_inds)
  shape <- 1 / CV^2
  C <- sapply(seq_len(K), function(k) rpois(N_inds, rgamma(N_inds, shape, rate = shape / LAMBDA)))
  colnames(C) <- paste0("t", round(t_mid))
  s2 <- SIGMA_BETA^2 * exp(-S * t_mid)
  g  <- rowSums(sqrt(sweep(C, 2, s2, `*`)) * matrix(rnorm(N_inds * K), N_inds, K))
  Y  <- g + rnorm(N_inds, 0, sqrt(var(g) * (1 - H2) / H2))
  list(C = C, Y = Y, s2_true = s2)
}

## ---- fit one dataset: REML, then log-linear and NLS second stages ----------
run_diagGREML <- function(S, REP, H2, CV, N_inds) {
  setDTthreads(1)
  d <- sim_diag(S, REP, H2, CV, N_inds)
  C <- d$C; m <- colMeans(C)
  if (any(!is.finite(m) | m <= 0)) stop("bad mean count")

  A <- setNames(lapply(seq_len(ncol(C)), function(k) C[, k] / m[k]), colnames(C))
  fit <- suppressMessages(fit_diagGREML(y = d$Y, X = matrix(1, nrow(C), 1), A = A,
                                        constraint = FALSE, magic0316 = TRUE))
  idx <- match(names(A), fit$Vlistnames)
  if (anyNA(idx)) stop("component name mismatch")

  res <- data.table(mid = t_mid, m = m,
                    est = as.numeric(fit$varcmp)[idx],
                    se  = sqrt(diag(fit$Hi))[idx])
  res[, `:=`(sigma2_b = est / m, se_b = se / m)]
  n_neg <- sum(!is.finite(res$est) | res$est <= 0)

  ## log-linear fit (drops non-positive components)
  pos   <- res[is.finite(sigma2_b) & sigma2_b > 0]
  m_log <- if (nrow(pos) >= 2) tryCatch(lm(log(sigma2_b) ~ mid, data = pos),
                                        error = function(e) NULL) else NULL

  ## NLS fit on the untransformed scale (keeps negatives)
  start_a <- if (!is.null(m_log)) exp(coef(m_log)[[1]]) else SIGMA_BETA^2
  start_s <- if (!is.null(m_log)) max(-coef(m_log)[[2]], 1e-8) else 1e-3
  if (!is.finite(start_a) || start_a <= 0) start_a <- SIGMA_BETA^2
  if (!is.finite(start_s) || start_s < 0)  start_s <- 1e-3

  m_nls <- tryCatch({
    nd <- res[is.finite(sigma2_b) & is.finite(se_b) & se_b > 0]
    if (nrow(nd) < 2) stop("too few points")
    nls(sigma2_b ~ a * exp(-s_par * mid), data = nd, weights = 1 / se_b^2,
        start = list(a = start_a, s_par = start_s),
        algorithm = "port", lower = c(a = 1e-12, s_par = 0))
  }, error = function(e) NULL)

  data.table(S = S, REP = REP, H2 = H2, CV = CV, N_inds = N_inds, n_neg = n_neg,
             sigma2_m_log = if (is.null(m_log)) NA_real_ else exp(coef(m_log)[[1]]),
             s_hat_log    = if (is.null(m_log)) NA_real_ else -coef(m_log)[[2]],
             sigma2_m_nls = if (is.null(m_nls)) NA_real_ else coef(m_nls)[["a"]],
             s_hat_nls    = if (is.null(m_nls)) NA_real_ else coef(m_nls)[["s_par"]])
}

## ---- run the full grid for one N, make plots and tables --------------------
run_for_N <- function(N_inds, n_workers = NULL) {

  if (is.null(n_workers))
    n_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = parallel::detectCores()))
  if (!is.finite(n_workers) || n_workers < 1) n_workers <- 1L

  params <- CJ(S = Svals, REP = seq_len(N_REPS), H2 = H2_vals, CV = CV_vals)[, N_inds := N_inds]

  plan(multicore, workers = as.integer(n_workers))
  res <- future_pmap(as.list(params), safely(run_diagGREML),
                     .options = furrr_options(seed = TRUE, scheduling = Inf,
                                              packages = "data.table"))

  errs <- keep(res, ~ !is.null(.x$error))
  if (length(errs)) warning(length(errs), "/", length(res), " failed at N = ", N_inds,
                            ". First: ", conditionMessage(errs[[1]]$error))

  ok <- compact(map(res, "result"))
  if (!length(ok)) stop("all runs failed at N = ", N_inds)
  res_total <- rbindlist(ok, use.names = TRUE, fill = TRUE)
  res_total[, S_f := factor(S, levels = Svals)]

  s_long <- melt(res_total, id.vars = c("S", "S_f", "H2", "CV", "N_inds"),
                 measure.vars = c("s_hat_log", "s_hat_nls"),
                 variable.name = "method", value.name = "s_hat")
  s_long[, method := fifelse(method == "s_hat_log", "log", "nls")]

  vm_long <- melt(res_total, id.vars = c("S", "S_f", "H2", "CV", "N_inds"),
                  measure.vars = c("sigma2_m_log", "sigma2_m_nls"),
                  variable.name = "method", value.name = "sigma2_m")
  vm_long[, method := fifelse(method == "sigma2_m_log", "log", "nls")]

  ttl <- paste0("N = ", format(N_inds, big.mark = ",", scientific = FALSE))

  ## p1: s_hat recovery
  p1 <- ggplot(s_long, aes(S_f, s_hat, colour = method)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(0.7)) +
    ggbeeswarm::geom_quasirandom(dodge.width = 0.7, size = 0.6, alpha = 0.5) +
    geom_point(data = unique(s_long[, .(S_f, S)]), aes(S_f, S),
              colour = "red", size = 3, inherit.aes = FALSE) +
    facet_grid(CV ~ H2, labeller = label_both, scales = "free_y") +
    labs(x = "True s", y = expression(hat(s)), title = ttl) + theme_bw()
   ## p2: sigma2_m recovery (log axis, clipped view so outliers don't dominate)
  p2 <- ggplot(vm_long[is.finite(sigma2_m) & sigma2_m > 0],
               aes(S_f, sigma2_m, colour = method)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(0.7)) +
    geom_hline(yintercept = SIGMA_BETA^2, colour = "red", linetype = "dashed", linewidth = 1) +
    facet_grid(CV ~ H2, labeller = label_both) +
    scale_y_log10(breaks = breaks_log(6), labels = label_number(accuracy = 0.01)) +
    coord_cartesian(ylim = c(0.1, 10)) +
    labs(x = "True s", y = expression(hat(sigma)[m]^2), title = ttl) + theme_bw()

  ## summary table: bias, SD, RMSE per cell
  summary_tbl <- s_long[, {
    z <- (s_hat - S)[is.finite(s_hat - S)]
    .(mean_s_hat = mean(s_hat, na.rm = TRUE),
      bias     = if (length(z)) mean(z) else NA_real_,
      rel_bias = if (length(z) && S != 0) mean(z) / S else NA_real_,
      sd       = sd(s_hat, na.rm = TRUE),
      rmse     = if (length(z)) sqrt(mean(z^2)) else NA_real_,
      se_bias  = if (length(z) >= 2) sd(z) / sqrt(length(z)) else NA_real_,
      n_valid  = length(z), n_total = .N)
  }, by = .(S, H2, CV, method)]
  setorder(summary_tbl, method, H2, CV, S)

  ## p3: bias in s_hat
  p3 <- ggplot(summary_tbl, aes(factor(S), bias, colour = method)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_linerange(aes(ymin = bias - sd, ymax = bias + sd),
                  position = position_dodge(0.4), alpha = 0.25, linewidth = 2) +
    geom_pointrange(aes(ymin = bias - 1.96 * se_bias, ymax = bias + 1.96 * se_bias),
                    position = position_dodge(0.4)) +
    facet_grid(CV ~ H2, labeller = label_both, scales = "free_y") +
    labs(x = "True s", y = expression(hat(s) - s),
        title = paste0("Bias in \u015d, ", ttl),
        subtitle = "Dark bar: 95% CI on the mean. Pale bar: SD across replicates.") +
    theme_bw()

  ## p4: RMSE of s_hat
  p4 <- ggplot(summary_tbl, aes(factor(S), rmse, colour = method, group = method)) +
    geom_point(position = position_dodge(0.4), size = 2) +
    geom_line(position = position_dodge(0.4), alpha = 0.5) +
    facet_grid(CV ~ H2, labeller = label_both, scales = "free_y") +
    labs(x = "True s", y = expression(RMSE(hat(s))),
         title = paste0("RMSE of \u015d, ", ttl)) + theme_bw()

  tag <- paste0("N", N_inds)
  ggsave(file.path(OUTDIR, paste0("s_hat_grid_",   tag, ".png")), p1, width = 10, height = 8)
  ggsave(file.path(OUTDIR, paste0("sigma2m_grid_", tag, ".png")), p2, width = 10, height = 8)
  ggsave(file.path(OUTDIR, paste0("s_hat_bias_",   tag, ".png")), p3, width = 10, height = 8)
  ggsave(file.path(OUTDIR, paste0("s_hat_rmse_",   tag, ".png")), p4, width = 10, height = 8)
  fwrite(summary_tbl, file.path(OUTDIR, paste0("s_hat_summary_", tag, ".csv")))

  list(res_total = res_total, summary_tbl = summary_tbl, errors = errs,
       p1 = p1, p2 = p2, p3 = p3, p4 = p4)
}

## ---- run all sample sizes --------------------------------------------------
all_results <- map(N_vals, run_for_N)
names(all_results) <- paste0("N", N_vals)
plan(sequential)

## ---- one combined table across all N --------------------------------------
all_summary <- rbindlist(map(all_results, "summary_tbl"), idcol = "N")
fwrite(all_summary, file.path(OUTDIR, "s_hat_summary_allN.csv"))

## the comparison that answers "does more N fix low heritability?"
all_summary[H2 == 0.1 & method == "nls", .(N, CV, S, bias, rmse)]
