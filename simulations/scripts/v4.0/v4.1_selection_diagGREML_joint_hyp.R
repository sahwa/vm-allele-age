library(data.table)
library(purrr)
library(ggplot2)
library(patchwork)
f = glue::glue

COMMON="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/common"
DATA="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v4.0"
FIGS="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v4.0"
OUT = file.path(DATA, "collated")
dir.create(OUT, showWarnings = FALSE)
source("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/diagGREML.R")
source(file.path(COMMON, "fit_sel_reml_joint.R"))

S_vals = c(0, 0.002, 0.005, 0.01)
N_REPS = 25
K = 50
n_sample = 2500
H2 = 0.5

get_counts = function(S, K, H2, REP) {
  set.seed(1e4 * which(S_vals == S) + REP)           # reproducible phenotype noise
  sing <- fread(file.path(DATA, f("sing_K{S}_seed{REP}.tsv")))   # carrier, age, beta
  stopifnot(max(sing$carrier) <= n_sample)

  # ---- Age bins: equal-count, and each age falls in exactly one bin ---------
  breaks <- unique(quantile(sing$age, probs = seq(0, 1, length.out = K + 1), type = 1))
  sing[, bin := cut(age, breaks, include.lowest = TRUE, labels = FALSE)]
  K_eff <- length(breaks) - 1                        # can be < K if many ages are tied

  # ---- Per-bin summaries, including the oracle (true mean beta^2) ----------
  bins <- sing[, .(mid = mean(age), oracle_b2 = mean(beta^2), n_sing = .N), keyby = bin]

  # ---- Genetic value and phenotype for everyone in the sample --------------
  g <- numeric(n_sample)                             # zero for people with no singletons
  agg <- rowsum(sing$beta, sing$carrier)             # sum of true betas per carrier
  g[as.integer(rownames(agg))] <- agg[, 1]
  Y <- g + rnorm(n_sample, sd = sqrt(var(g) * (1 - H2) / H2))

  # ---- Counts: every person x every bin, rows aligned with Y ---------------
  C <- table(factor(sing$carrier, levels = seq_len(n_sample)),
             factor(sing$bin,     levels = seq_len(K_eff)))
  C <- matrix(C, nrow = n_sample)
  m <- colMeans(C)
  A <- setNames(lapply(seq_len(K_eff), function(k) C[, k] / m[k]),
                paste0("bin", seq_len(K_eff)))
  list(A = A, C = C, Y = Y, sing = sing, t = bins$mid, bins = bins, m = m,
       K_eff = K_eff, y_var = var(Y))
}

## ---- two-stage: diagGREML per bin, then log-linear and NLS on the ages ------
run_diagGREML = function(d, K_S, REP, H2) {
  bins <- copy(d$bins)
  fit <- suppressMessages(fit_diagGREML(y = d$Y, A = d$A, X = matrix(1, n_sample, 1),
                                        constraint = FALSE, magic0316 = TRUE))
  idx <- match(names(d$A), fit$Vlistnames)
  bins[, `:=`(m = d$m,
              est = as.numeric(fit$varcmp)[idx],
              se  = sqrt(diag(fit$Hi))[idx])]
  bins[, `:=`(sigma2_b = est / m, se_b = se / m)]    # per-mutation variance and its SE

  m_orc <- lm(log(oracle_b2) ~ mid, data = bins)

  pos   <- bins[is.finite(est) & est > 0]
  m_log <- if (nrow(pos) >= 3) lm(log(sigma2_b) ~ mid, data = pos) else NULL
  m_nls <- if (is.null(m_log)) NULL else tryCatch(
    nls(sigma2_b ~ a * exp(-s_par * mid), data = bins, weights = 1 / se_b^2,
        start = list(a = exp(coef(m_log)[[1]]), s_par = max(-coef(m_log)[[2]], 1e-4))),
    error = function(e) NULL)

  #### V_M ####
  new_per_person  <- d$sing[age == 0, .N] / n_sample  # empirical: should be close to 2 * u_target
  oracle_sigma2_m <- exp(coef(m_orc)[[1]])
  sigma2_m_log    <- if (is.null(m_log)) NA_real_ else exp(coef(m_log)[[1]])

  list(
    summary = data.table(
      K_S = K_S, REP = REP, H2 = H2, K_eff = d$K_eff,
      s_pred          = K_S,                           # predicted slope: K_SEL * sigma^2, sigma^2 = 1
      new_per_person  = new_per_person,                # check against 2 * u_target
      oracle_sigma2_m = oracle_sigma2_m,
      oracle_s_hat    = -coef(m_orc)[[2]],
      V_M_oracle      = oracle_sigma2_m * new_per_person,
      n_neg        = sum(bins$est <= 0),
      sigma2_m_log = sigma2_m_log,
      s_hat_log    = if (is.null(m_log)) NA_real_ else -coef(m_log)[[2]],
      V_M_reml     = sigma2_m_log * new_per_person,
      sigma2_m_nls = if (is.null(m_nls)) NA_real_ else coef(m_nls)[["a"]],
      s_hat_nls    = if (is.null(m_nls)) NA_real_ else coef(m_nls)[["s_par"]]
    ),
    bins = bins[, .(bin, mid, n_sing, oracle_b2, sigma2_b, se_b)]
  )
}

## ---- one (K_SEL, seed): all fits, returned as one row plus a bins table -----
run_once = function(S, REP) {
  d  <- get_counts(S, K, H2, REP)
  dg <- run_diagGREML(d, K_S = S, REP = REP, H2 = H2)
  fit_exp <- tryCatch(fit_sel_reml(y = d$Y, C = d$C, t = d$t, curve = "exp"),
                      error = function(e) NULL)
  fit_hyp <- tryCatch(fit_sel_reml(y = d$Y, C = d$C, t = d$t, curve = "hyp"),
                      error = function(e) NULL)

  gj <- function(m, x) if (is.null(m)) NA_real_ else as.numeric(m[[x]])
  joint <- data.table(
    s_hat_joint_exp    = gj(fit_exp, "s"),        s_se_joint_exp    = gj(fit_exp, "s_se"),
    sigma2_m_joint_exp = gj(fit_exp, "sigma2_m"), sigma2_t0_exp     = gj(fit_exp, "sigma2_t0"),
    ll_joint_exp       = gj(fit_exp, "logLik"),   bound_joint_exp   = gj(fit_exp, "at_bound"),
    s_hat_joint_hyp    = gj(fit_hyp, "s"),        s_se_joint_hyp    = gj(fit_hyp, "s_se"),
    sigma2_m_joint_hyp = gj(fit_hyp, "sigma2_m"), sigma2_t0_hyp     = gj(fit_hyp, "sigma2_t0"),
    ll_joint_hyp       = gj(fit_hyp, "logLik"),   bound_joint_hyp   = gj(fit_hyp, "at_bound")
  )
  joint[, `:=`(V_M_joint_exp = sigma2_m_joint_exp * dg$summary$new_per_person,
               V_M_joint_hyp = sigma2_m_joint_hyp * dg$summary$new_per_person,
               hyp_wins      = ll_joint_hyp > ll_joint_exp)]

  list(rep  = cbind(dg$summary, joint, y_var = d$y_var),
       bins = cbind(K_S = S, REP = REP, dg$bins))
}

## ---- run everything ---------------------------------------------------------
params <- CJ(S = S_vals, REP = seq_len(N_REPS))

res <- pmap(params, function(S, REP) tryCatch(run_once(S, REP), error = function(e) {
  message(f("K_S = {S}, REP = {REP} failed: {conditionMessage(e)}")); NULL }))
ok <- compact(res)

## ---- collate ----------------------------------------------------------------
reps <- rbindlist(map(ok, "rep"),  use.names = TRUE, fill = TRUE)   # one row per (K_S, REP)
bins <- rbindlist(map(ok, "bins"), use.names = TRUE, fill = TRUE)   # one row per (K_S, REP, bin)



## long: one row per (K_S, REP, method), same columns for every method
long <- rbind(
  reps[, .(K_S, REP, method = "oracle",     s_hat = oracle_s_hat,    sigma2_m = oracle_sigma2_m,    V_M = V_M_oracle,                    V_P = y_var)],
  reps[, .(K_S, REP, method = "two-stage",  s_hat = s_hat_log,       sigma2_m = sigma2_m_log,       V_M = V_M_reml,                      V_P = y_var)],
  reps[, .(K_S, REP, method = "nls",        s_hat = s_hat_nls,       sigma2_m = sigma2_m_nls,       V_M = sigma2_m_nls * new_per_person, V_P = y_var)],
  reps[, .(K_S, REP, method = "joint exp",  s_hat = s_hat_joint_exp, sigma2_m = sigma2_m_joint_exp, V_M = V_M_joint_exp,                 V_P = y_var)],
  reps[, .(K_S, REP, method = "joint hyp",  s_hat = s_hat_joint_hyp, sigma2_m = sigma2_m_joint_hyp, V_M = V_M_joint_hyp,                 V_P = y_var)])
long[, method := factor(method, levels = c("oracle", "two-stage", "nls", "joint exp", "joint hyp"))]
long[, h2m := V_M / V_P]                                            # mutational heritability

## summary: one row per (K_S, method)
summ <- long[, .(n_fit           = sum(is.finite(s_hat)),
                 s_median        = median(s_hat, na.rm = TRUE),
                 s_mad           = mad(s_hat, na.rm = TRUE),
                 sigma2_m_median = median(sigma2_m, na.rm = TRUE),
                 sigma2_m_mad    = mad(sigma2_m, na.rm = TRUE),
                 V_M_median      = median(V_M, na.rm = TRUE),
                 h2m_median      = median(h2m, na.rm = TRUE),
                 h2m_mad         = mad(h2m, na.rm = TRUE)),
             keyby = .(K_S, method)]

## how often the data prefer the hyperbolic curve
curve_pick <- reps[, .(pct_hyp_wins = 100 * mean(hyp_wins, na.rm = TRUE),
                       n_bound_exp  = sum(bound_joint_exp > 0, na.rm = TRUE),
                       n_bound_hyp  = sum(bound_joint_hyp > 0, na.rm = TRUE)), keyby = K_S]
print(curve_pick)

### plot V_M / V_P

est   <- long[method != "oracle"]
truth <- long[method == "oracle", .(h2m = median(h2m, na.rm = TRUE)), by = K_S]
med   <- est[, .(h2m = median(h2m, na.rm = TRUE)), by = .(K_S, method)]
cols  <- c("two-stage" = "grey25", "nls" = "#E69F00",
           "joint exp" = "#2C7BB6", "joint hyp" = "#3E3A6E")
ylims <- c(-4, 16) * median(truth$h2m)               # same window as before, relative to the oracle

p1_h2m <- ggplot(est[is.finite(h2m)], aes(method, h2m, colour = method)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_hline(data = truth, aes(yintercept = h2m), colour = "#D7191C", linetype = "dashed") +
  ggbeeswarm::geom_quasirandom(size = 1.2, alpha = 0.6, width = 0.3) +
  geom_crossbar(data = med, aes(ymin = h2m, ymax = h2m), width = 0.6, linewidth = 0.4) +
  facet_wrap(~ K_S, nrow = 1,
             labeller = labeller(K_S = function(x) paste0("K_SEL = ", x))) +
  scale_colour_manual(values = cols, guide = "none") +
  coord_cartesian(ylim = ylims) +
  labs(x = NULL, y = expression(hat(V)[M] / V[P]),
       title = "Mutational heritability: estimates against the SLiM oracle",
       subtitle = "Dashed red = oracle (truth); bars = medians; view clipped to -4x to 16x the oracle") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(FIGS, "h2m_single_joint_oracle.png"), p1_h2m, width = 11, height = 4.5, dpi = 150)



## error of each estimate against the oracle from the same run
orc <- long[method == "oracle", .(K_S, REP, h2m_true = h2m, s_true = s_hat)]
err <- merge(long[method != "oracle"], orc, by = c("K_S", "REP"))
err <- melt(err[, .(K_S, REP, method, h2m = h2m - h2m_true, s = s_hat - s_true)],
            id.vars = c("K_S", "REP", "method"),
            variable.name = "quantity", value.name = "err")

perf <- err[is.finite(err),
            .(n = .N, bias = mean(err), se = sd(err) / sqrt(.N), rmse = sqrt(mean(err^2))),
            by = .(K_S, method, quantity)]

perf_long <- melt(perf, id.vars = c("K_S", "method", "quantity", "n", "se"),
                  measure.vars = c("bias", "rmse"), variable.name = "metric")
perf_long[metric == "rmse", se := NA_real_]          # the interval applies to bias only

dg <- position_dodge(width = 0.5)

p_perf <- ggplot(perf_long, aes(method, value, colour = method, shape = metric)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_linerange(aes(ymin = value - 1.96 * se, ymax = value + 1.96 * se),
                 position = dg, na.rm = TRUE) +
  geom_point(size = 2.5, position = dg) +
  geom_text(data = perf, aes(method, -Inf, label = n), inherit.aes = FALSE,
            vjust = -0.6, size = 3, colour = "grey40") +
  facet_grid(quantity ~ K_S, scales = "free_y",
             labeller = labeller(K_S = function(x) paste0("K_SEL = ", x),
                                 quantity = c(h2m = "V_M / V_P", s = "s"))) +
  scale_colour_manual(values = cols, guide = "none") +
  scale_shape_manual(values = c(bias = 16, rmse = 17),
                     labels = c(bias = "Bias (95% CI)", rmse = "RMSE"), name = NULL) +
  labs(x = NULL, y = "Error against the oracle",
       title = "Bias and RMSE by method",
       subtitle = "Grey numbers = fits contributing to each point") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(FIGS, "bias_rmse_by_method.png"), p_perf, width = 11, height = 6, dpi = 150)


fwrite(reps, file.path(OUT, "reps_wide.csv"))
fwrite(long, file.path(OUT, "reps_long.csv"))
fwrite(bins, file.path(OUT, "bins.csv"))
fwrite(summ, file.path(OUT, "summary.csv"))

print(summ)
