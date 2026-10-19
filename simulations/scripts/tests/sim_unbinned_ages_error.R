#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Bias sweeps for the joint singleton estimator, as a SLURM array.
##
##   Rscript sim_ubinned_ages_error.R grid              # write the cell table
##   Rscript sim_ubinned_ages_error.R cell <ID>         # run one cell
##   Rscript sim_ubinned_ages_error.R collect           # combine, summarise, plot
##
## Two experiments, deliberately not crossed:
##   A) BINNING    sweep s x K, no age error.  Plotted against the dimensionless
##                 group that governs the Jensen bias, not against K, so one
##                 curve covers every parameter setting.
##   B) AGE ERROR  ungrouped fit, noise and bias components separated.
##
## Settings are tuned for BIAS DETECTION, not realism: h2 and CV are pushed high
## so the sampling SE is small enough that structural bias is visible.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

COMMON_DIR <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/common"
source(file.path(COMMON_DIR, "fit_sel_reml_joint.R"))
source(file.path(COMMON_DIR, "age_error.R"))
source(file.path(COMMON_DIR, "sim_fast.R"))

DATA <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/tests/sim_unbinned_ages_err"
FIGS <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/tests/sim_unbinned_ages_err"
CELLS <- file.path(DATA, "cells")
OUT   <- file.path(DATA, "out")
for (d in c(DATA, FIGS, CELLS, OUT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

AGES_FILE <- file.path(
  "/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations",
  "data/v1.2/n500000_L1e+08/sing_ages.npy"
)

## --- parameters -------------------------------------------------------------
N_REPS     <- 50
N_inds     <- 1e5          # structural bias is O(1) in N; buy precision with reps
H2         <- 0.8          # high on purpose: shrinks SE so bias is visible
CV_SHARED  <- 0.65         # ditto. SE(a)/a ~ sqrt(2/N)/(CV * h2_sing)
N_SING     <- 1000
SIGMA_BETA <- 1
CURVE      <- "hyp"
AGE_CAP_Q  <- 0.999        # 0.9999 gave a 1092-generation grid for 0.01% of data
N_AGE_BINS <- 60           # log-spaced: matches both the skew and the error kernel

## s must be large enough to be identifiable at all, or s_hat is noise and the
## v_ratio is meaningless. The low-x end of the collapse curve is reached with
## large K instead, not with small s.
S_VALS  <- c(1e-3, 3e-3, 1e-2, 3e-2)
K_VALS  <- c(2, 4, 8, 20, 40)        # plus NA = ungrouped reference
S_ERR   <- 1e-2                      # age-error sweep: needs signal to attenuate

## --- age distribution -------------------------------------------------------
## Shape is a property of the coalescent and demography, not sequence length, so
## the L = 1e8 distribution carries over; only the per-person count is rescaled.
## Log spacing because the decay is flat within a generation and dating error is
## multiplicative, so equal-ratio bins carry roughly equal information.
## CHECK: that AGES is in generations, not years.
build_age_pmf <- function(ages, n_bins = N_AGE_BINS, cap_q = AGE_CAP_Q) {
  ages  <- ages[is.finite(ages) & ages > 0]
  cap   <- quantile(ages, cap_q)
  edges <- unique(c(0, exp(seq(0, log(cap), length.out = n_bins))))
  mid   <- (head(edges, -1) + tail(edges, -1)) / 2
  cnt   <- tabulate(findInterval(pmin(ages, cap), edges, all.inside = TRUE),
                    nbins = length(mid))
  keep  <- cnt > 0
  list(t = mid[keep], p = cnt[keep] / sum(cnt))
}

get_pmf <- function() build_age_pmf(RcppCNPy::npyLoad(AGES_FILE))

## --- decay curve ------------------------------------------------------------
fdecay <- function(s, t, curve = CURVE)
  if (curve == "exp") exp(-s * t) else 1 / (1 + s * t)

## --- merge to K coarse bins -------------------------------------------------
## Returns the count-weighted within-bin age VARIANCE as well: for a skewed age
## distribution this, not the nominal width, is what sets the Jensen bias.
merge_bins <- function(C_fine, t_fine, K) {
  j  <- split(seq_along(t_fine), cut(seq_along(t_fine), breaks = K, labels = FALSE))
  cm <- colMeans(C_fine)
  C  <- vapply(j, function(jj) rowSums(C_fine[, jj, drop = FALSE]), numeric(nrow(C_fine)))
  tt <- vapply(j, function(jj) sum(cm[jj] * t_fine[jj]) / sum(cm[jj]), numeric(1))
  vv <- vapply(seq_along(j), function(b) {
    jj <- j[[b]]; sum(cm[jj] * (t_fine[jj] - tt[b])^2) / sum(cm[jj])
  }, numeric(1))
  m  <- vapply(j, function(jj) sum(cm[jj]), numeric(1))
  list(C = C, t = as.vector(tt), vbar = sum(m * vv) / sum(m))
}

## Dimensionless group governing the within-bin averaging bias.
## mean(f)/f(t0) ~ 1 + V*f''(t0)/(2 f(t0)) and sinh(x)/x ~ 1 + x^2/6, so
## x = sqrt(3 V f''/f). For exp with uniform ages this is exactly s*w/2.
jensen_x <- function(s, t0, vbar, curve = CURVE) {
  d2_over_f <- if (curve == "exp") s^2 else 2 * s^2 / (1 + s * t0)^2
  sqrt(3 * vbar * d2_over_f)
}

## --- fit helpers ------------------------------------------------------------
fit1 <- function(y, C, t, s_max) {
  tryCatch(fit_sel_reml(y = y, C = C, t = t, curve = CURVE, s_max = s_max),
           error = function(e) NULL)
}
grab <- function(ft) data.table(
  s_hat    = if (is.null(ft)) NA_real_ else ft$s,
  s_se     = if (is.null(ft)) NA_real_ else ft$s_se,
  sigma2_m = if (is.null(ft)) NA_real_ else ft$sigma2_m,
  h2_sing  = if (is.null(ft)) NA_real_ else ft$h2_sing,
  lrt_p    = if (is.null(ft)) NA_real_ else ft$lrt_p,
  at_bound = if (is.null(ft)) NA       else ft$at_bound
)

## ===========================================================================
## GRID
## ===========================================================================
grid_file <- file.path(DATA, "cells.csv")

do_grid <- function() {
  bins <- CJ(sweep = "bin", REP = seq_len(N_REPS), S = S_VALS, err = NA_character_)
  errs <- CJ(sweep = "err", REP = seq_len(N_REPS), S = S_ERR,  err = names(ERR_GRID))
  cells <- rbind(bins, errs)
  cells[, cell := .I]
  setcolorder(cells, "cell")
  fwrite(cells, grid_file)
  pmf <- get_pmf()
  cat(sprintf("%d cells written to %s\n", nrow(cells), grid_file))
  cat(sprintf("age grid: %d bins, %.1f-%.0f gen, mean %.1f\n",
              length(pmf$t), min(pmf$t), max(pmf$t), sum(pmf$t * pmf$p)))
  cat(sprintf("submit with: sbatch --array=1-%d\n", nrow(cells)))
}

## ===========================================================================
## CELL
## ===========================================================================
## One cell = one replicate of one (sweep, parameter) combination. For "bin"
## that is all K on a single simulated dataset, so the K fits stay paired within
## replicate and common random numbers cancel the shared sampling noise.
do_cell <- function(id) {
  cells <- fread(grid_file)
  cl <- cells[cell == id]
  if (!nrow(cl)) stop("no such cell: ", id)
  pmf <- get_pmf()

  out <- if (cl$sweep == "bin") {
    d  <- sim_counts(cl$REP, pmf, cl$S)              # P = NULL: no age error
    Y  <- make_pheno(d)
    sm <- max(0.1, 10 * cl$S)                        # keep fits off the upper bound
    rbindlist(lapply(c(NA_integer_, K_VALS), function(K) {
      m <- if (is.na(K)) list(C = d$C, t = d$t, vbar = 0) else merge_bins(d$C, d$t, K)
      t0 <- sum(m$t * colMeans(m$C)) / sum(colMeans(m$C))
      cbind(data.table(cell = id, sweep = "bin", REP = cl$REP, S = cl$S,
                       K = K, err = NA_character_, x = jensen_x(cl$S, t0, m$vbar)),
            grab(fit1(Y, m$C, m$t, sm)))
    }))
  } else {
    P <- if (cl$err == "none") NULL else
      build_err_matrices(pmf$t, t_ref = sum(pmf$t * pmf$p))[[cl$err]]
    d <- sim_counts(cl$REP, pmf, cl$S, P = P)
    cbind(data.table(cell = id, sweep = "err", REP = cl$REP, S = cl$S,
                     K = NA_integer_, err = cl$err, x = NA_real_),
          grab(fit1(make_pheno(d), d$C, d$t, 0.1)))
  }

  fwrite(out, file.path(CELLS, sprintf("cell_%05d.csv", id)))
  cat(sprintf("cell %d (%s) done, %d rows\n", id, cl$sweep, nrow(out)))
}

## ===========================================================================
## COLLECT
## ===========================================================================
do_collect <- function() {
  library(ggplot2); library(patchwork)
  suppressWarnings(try(source(file.path(COMMON_DIR, "themes/zoo_theme.R")), silent = TRUE))
  if (!exists("theme_zoo")) theme_zoo <- function() theme_minimal()

  f <- list.files(CELLS, pattern = "^cell_.*csv$", full.names = TRUE)
  expected <- nrow(fread(grid_file))
  cat(sprintf("%d of %d cells present\n", length(f), expected))
  res <- rbindlist(lapply(f, fread), fill = TRUE)
  fwrite(res, file.path(OUT, "all_cells.csv"))

  ## --- binning: paired ratio to the ungrouped fit of the SAME replicate -----
  rb  <- res[sweep == "bin"]
  ref <- rb[is.na(K), .(REP, S, s_ref = s_hat, v_ref = sigma2_m)]
  rb  <- merge(rb[!is.na(K)], ref, by = c("REP", "S"))
  rb[, `:=`(s_ratio = s_hat / s_ref, v_ratio = sigma2_m / v_ref)]

  summ_bin <- rb[, .(n = .N, n_bound = sum(at_bound, na.rm = TRUE),
                     x = mean(x, na.rm = TRUE),
                     s_ratio = mean(s_ratio, na.rm = TRUE),
                     v_ratio = mean(v_ratio, na.rm = TRUE),
                     v_ratio_med = median(v_ratio, na.rm = TRUE)),
                 keyby = .(S, K)]
  print(summ_bin); fwrite(summ_bin, file.path(OUT, "summ_binning.csv"))

  ## --- age error ------------------------------------------------------------
  re <- res[sweep == "err"]
  summ_err <- re[, .(n = .N, n_bound = sum(at_bound, na.rm = TRUE),
                     s_mean = mean(s_hat, na.rm = TRUE),
                     s_sd   = sd(s_hat, na.rm = TRUE),
                     attenuation = mean(s_hat, na.rm = TRUE) / S_ERR,
                     v_mean = mean(sigma2_m, na.rm = TRUE),
                     v_bias_pct = 100 * (mean(sigma2_m, na.rm = TRUE) / SIGMA_BETA^2 - 1)),
                 keyby = err]
  print(summ_err); fwrite(summ_err, file.path(OUT, "summ_age_error.csv"))

  ## --- plots ----------------------------------------------------------------
  ## The collapse: every (s, K) on one curve, with sinh(x)/x overlaid.
  xs <- 10^seq(-2.5, 0.8, length.out = 200)
  p1 <- ggplot(summ_bin, aes(x, v_ratio)) +
    geom_line(data = data.frame(x = xs, v_ratio = sinh(xs) / xs),
              colour = "red", linetype = "dashed") +
    geom_point(aes(shape = factor(K)), size = 2.5) +
    scale_x_log10() + scale_y_log10() +
    labs(x = "x  (dimensionless bin width)",
         y = expression(hat(sigma)[m]^2 ~ "/ ungrouped"), shape = "K") + theme_zoo()

  p2 <- ggplot(summ_bin, aes(x, s_ratio, shape = factor(K))) +
    geom_hline(yintercept = 1, colour = "red", linetype = "dashed") +
    geom_point(size = 2.5) + scale_x_log10() +
    labs(x = "x", y = expression(hat(s) ~ "/ ungrouped"), shape = "K") + theme_zoo()

  re[, err := factor(err, levels = names(ERR_GRID))]
  p3 <- ggplot(re, aes(err, s_hat)) +
    ggbeeswarm::geom_quasirandom(size = 0.8) +
    geom_hline(yintercept = S_ERR, colour = "red", linetype = "dashed") +
    labs(x = NULL, y = expression(hat(s))) + theme_zoo() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(FIGS, "binning_collapse.png"), p1 / p2, width = 7, height = 8, dpi = 200)
  ggsave(file.path(FIGS, "age_error_decomposition.png"), p3, width = 7, height = 4.5, dpi = 200)
  cat("figures written to", FIGS, "\n")
}

## ===========================================================================
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "grid"
switch(mode,
       grid    = do_grid(),
       cell    = do_cell(as.integer(args[2])),
       collect = do_collect(),
       stop("mode must be one of: grid, cell, collect"))
