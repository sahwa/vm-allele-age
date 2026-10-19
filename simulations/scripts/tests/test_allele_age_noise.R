library(parallel)
library(data.table)
library(patchwork)

COMMON_DIR="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/common"
source(file.path(COMMON_DIR, "fit_sel_reml_joint.R"))
source("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/common/themes/zoo_theme.R")

FIGS="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/tests/test_allele_age_noise"

N_CORES    <- parallel::detectCores() - 4
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
merge_bins <- function(C_fine, t_fine, K) {
  grp <- rep(seq_len(K), each = ncol(C_fine) / K)
  C   <- sapply(split(seq_len(ncol(C_fine)), grp),
                function(j) rowSums(C_fine[, j, drop = FALSE]))
  cm  <- colMeans(C_fine)
  t   <- as.vector(tapply(cm * t_fine, grp, sum) / tapply(cm, grp, sum))
  list(C = C, t = t)
}


## create a matrix of miscopying probabilities
make_prob_matrix <- function(K, sigma = 1) {
  P <- exp(-outer(1:K, 1:K, "-")^2 / (2 * sigma^2))
  P / rowSums(P)
}

misclassify <- function(C, P, tol = 1e-12) {
  Cobs <- matrix(0L, nrow(C), ncol(C))
  for (f in seq_len(ncol(C))) {
    p   <- P[f, ]
    sup <- which(p > tol)                  # observed bins this true bin can reach
    p   <- p[sup] / sum(p[sup])
    rem <- C[, f]                          # each person's singletons still to place
    p_left <- 1
    for (k in seq_along(sup)) {
      x <- if (k == length(sup)) rem else rbinom(length(rem), rem, min(1, p[k] / p_left))
      Cobs[, sup[k]] <- Cobs[, sup[k]] + x
      rem    <- rem - x
      p_left <- p_left - p[k]
    }
  }
  Cobs
}

pmatrix_fine = make_prob_matrix(K=K_FINE)

res <- mclapply(seq_len(N_REPS), function(REP) {
  print(REP)
  d     <- sim_fine(REP)
  C_obs <- misclassify(d$C, pmatrix_fine)          # once per rep
  rbindlist(lapply(K_vals, function(K) {
    m  <- merge_bins(C_obs, d$t_fine, K)           # or d$C for the error-free comparison
    ft <- fit_sel_reml(y = d$Y, C = m$C, t = m$t, curve = "exp")
      data.table(
        REP = REP, K = K,
        s_hat    = if (is.null(ft)) NA_real_ else ft$s,
        s_se     = if (is.null(ft)) NA_real_ else ft$s_se,
        sigma2_m = if (is.null(ft)) NA_real_ else ft$sigma2_m,
        at_bound = if (is.null(ft)) NA else ft$at_bound,
        y_var = var(d$Y)
      )
  }))
}, mc.cores = N_CORES)
res <- rbindlist(res)

fwrite(res, "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/tests/test_allele_age_noise/all_reps_K_res.txt")
res = fread("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/tests/test_allele_age_noise/all_reps_K_res.txt")

summ <- res[, .(n_fit             = sum(is.finite(s_hat)),
                s_median          = median(s_hat, na.rm = TRUE),
                s_bias_pct        = 100 * (median(s_hat, na.rm = TRUE) / S - 1),
                sigma2_m_median   = median(sigma2_m, na.rm = TRUE),
                sigma2_m_bias_pct = 100 * (median(sigma2_m, na.rm = TRUE) / SIGMA_BETA^2 - 1)),
            keyby = K]
print(summ)


res[, K := factor(K, levels = unique(sort(K)))]

p_s = res %>%
  ggplot(aes(x=K, y=s_hat)) +
  ggbeeswarm::geom_quasirandom() +
  geom_hline(yintercept = S, colour = 'red', linetype = 'dashed') +
  theme_zoo()

p_s_file = file.path(FIGS, "S_est_per_K.png")
ggsave(p_s_file, p_s)

p_v_m = res %>%
  ggplot(aes(x=K, y=sigma2_m)) +
  ggbeeswarm::geom_quasirandom() +
  geom_hline(yintercept = 1, colour = 'red', linetype = 'dashed') +
  scale_y_log10()

p_combined = p_s + p_v_m

p_combined_filed = file.path(FIGS, "v_m_s_est_per_K.png")
ggsave(p_combined_filed, p_combined)
