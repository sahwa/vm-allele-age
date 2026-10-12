library(parallel)
library(data.table)

COMMON_DIR="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/common"
source(file.path(COMMON_DIR, "fit_sel_reml_joint.R"))

N_CORES    <- 12
N_REPS     <- 30
N_inds     <- 5e5
H2         <- 0.5
S          <- 3e-2
t_max      <- 400
CV         <- 0.65
LAMBDA_TOT <- 1000                             # singletons per person, all ages
SIGMA_BETA <- 1
K_FINE     <- 120                              # resolution the data are generated at
K_vals     <- c(2, 3, 4, 6, 8, 12, 20, 40, 120)
stopifnot(K_FINE %% K_vals == 0)

## simulate on the fine grid: each fine bin has its own age and effect variance
sim_fine <- function(REP) {
  set.seed(REP)
  t_fine <- (seq_len(K_FINE) - 0.5) * t_max / K_FINE
  s2     <- SIGMA_BETA^2 * exp(-S * t_fine)
  shape  <- 1 / CV^2
  lam    <- LAMBDA_TOT / K_FINE
  g <- numeric(N_inds)
  C <- matrix(0L, N_inds, K_FINE)
  for (f in seq_len(K_FINE)) {
    C[, f] <- rpois(N_inds, rgamma(N_inds, shape, rate = shape / lam))
    g <- g + sqrt(C[, f] * s2[f]) * rnorm(N_inds)
  }
  Y <- g + rnorm(N_inds, 0, sqrt(var(g) * (1 - H2) / H2))
  list(C = C, Y = Y, t_fine = t_fine)
}

## merge neighbouring fine bins into K coarse ones, with count-weighted mean ages
merge_bins <- function(d, K) {
  grp <- rep(seq_len(K), each = K_FINE / K)
  C   <- sapply(split(seq_len(K_FINE), grp), function(j) rowSums(d$C[, j, drop = FALSE]))
  cm  <- colMeans(d$C)
  t   <- as.vector(tapply(cm * d$t_fine, grp, sum) / tapply(cm, grp, sum))
  list(C = C, t = t)
}

res <- mclapply(seq_len(N_REPS), function(REP) {
  print(REP)
  d <- sim_fine(REP)
  rbindlist(lapply(K_vals, function(K) {
    m  <- merge_bins(d, K)
    ft <- tryCatch(fit_sel_reml(y = d$Y, C = m$C, t = m$t, curve = "exp"),
                   error = function(e) NULL)
    data.table(REP = REP, K = K,
               s_hat    = if (is.null(ft)) NA_real_ else ft$s,
               s_se     = if (is.null(ft)) NA_real_ else ft$s_se,
               sigma2_m = if (is.null(ft)) NA_real_ else ft$sigma2_m,
               at_bound = if (is.null(ft)) NA else ft$at_bound)
  }))
}, mc.cores = N_CORES)
res <- rbindlist(res)

summ <- res[, .(n_fit             = sum(is.finite(s_hat)),
                s_median          = median(s_hat, na.rm = TRUE),
                s_bias_pct        = 100 * (median(s_hat, na.rm = TRUE) / S - 1),
                sigma2_m_median   = median(sigma2_m, na.rm = TRUE),
                sigma2_m_bias_pct = 100 * (median(sigma2_m, na.rm = TRUE) / SIGMA_BETA^2 - 1)),
            keyby = K]
print(summ)
