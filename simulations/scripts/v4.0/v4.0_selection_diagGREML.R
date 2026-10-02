library(data.table)
library(purrr)
library(future)
library(furrr)
library(patchwork)
f = glue::glue

source("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/diagGREML.R")
DATA="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v4.0"
FIGS="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v4.0"

n_sample = 2500

run_diagGREML <- run_diagGREML <- function(K_S, REP, n_sample = 2500, K = 8, H2 = 0.5, seed = 1) {
  set.seed(seed)
  counts_file = file.path(DATA, f("sing_K{K_S}_seed{REP}.tsv"))
  sing <- fread(counts_file)                         # one row per singleton: carrier, age, beta
  stopifnot(max(sing$carrier) <= n_sample)           # n_sample must match N_SAMPLE in the SLiM run

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

  # ---- REML -----------------------------------------------------------------
  fit <- fit_diagGREML(y = Y, A = A, X = matrix(1, n_sample, 1),
                                        constraint = FALSE, magic0316 = TRUE)
  idx <- match(names(A), fit$Vlistnames)
  bins[, `:=`(m = m,
              est = as.numeric(fit$varcmp)[idx],
              se  = sqrt(diag(fit$Hi))[idx])]
  bins[, `:=`(sigma2_b = est / m, se_b = se / m)]    # per-mutation variance and its SE

  m_orc <- lm(log(oracle_b2) ~ mid, data = bins)

  # ---- Second stage: log-linear (drops negatives) and NLS (keeps them) -----
  pos   <- bins[est > 0]
  m_log <- lm(log(sigma2_b) ~ mid, data = pos)

  m_nls <- tryCatch(
    nls(sigma2_b ~ a * exp(-s_par * mid), data = bins, weights = 1 / se_b^2,
        start = list(a = exp(coef(m_log)[[1]]), s_par = max(-coef(m_log)[[2]], 1e-4))),
    error = function(e) NULL)

  #### V_M ####
  new_per_person <- sing[age == 0, .N] / n_sample     # empirical: should be close to 2 * u_target

  oracle_sigma2_m <- exp(coef(m_orc)[[1]])
  sigma2_m_log    <- exp(coef(m_log)[[1]])

  V_M_oracle <- oracle_sigma2_m * new_per_person
  V_M_reml   <- sigma2_m_log    * new_per_person

  list(
    summary = data.table(
      K_S = K_S, REP = REP, H2 = H2, K_eff = K_eff,
      s_pred          = K_S,                           # predicted slope: K_SEL * sigma^2, sigma^2 = 1
      new_per_person  = new_per_person,                 # check against 2 * u_target
      oracle_sigma2_m = oracle_sigma2_m,
      oracle_s_hat    = -coef(m_orc)[[2]],
      V_M_oracle      = V_M_oracle,
      n_neg        = sum(bins$est <= 0),
      sigma2_m_log = sigma2_m_log, s_hat_log = -coef(m_log)[[2]],
      V_M_reml     = V_M_reml,
      sigma2_m_nls = if (is.null(m_nls)) NA_real_ else coef(m_nls)[["a"]],
      s_hat_nls    = if (is.null(m_nls)) NA_real_ else coef(m_nls)[["s_par"]]
    ),
    bins = bins[, .(bin, mid, n_sing, oracle_b2, sigma2_b, se_b)]
  )
}

MU = 1.28e-8; PI_TARGET = 1; L=1e8; SIGMA_BETA = 1
V_M_true <- 2 * MU * PI_TARGET * L * SIGMA_BETA^2   # from the SLiM parameters, same for every K_SEL


K_S_vals = c(0, 0.002, 0.005, 0.01)
REP_Vals = 2:25
n_sample = 2500
params = CJ(K_S = K_S_vals, REP = REP_Vals)
res <- pmap(params, possibly(run_diagGREML, otherwise = NULL))
summ <- rbindlist(map(compact(res), "summary"))

##### V_M #####
agg_vm <- summ[, .(
  true      = V_M_true,
  npp       = mean(new_per_person),
  oracle    = mean(V_M_oracle), oracle_se = sd(V_M_oracle) / sqrt(.N),
  reml      = mean(V_M_reml),   reml_se   = sd(V_M_reml)   / sqrt(.N)
), by = K_S]

plot_vm <- rbind(
  agg_vm[, .(K_S, what = "True (2 * u_target * sigma^2)", est = true,   se = 0)],
  agg_vm[, .(K_S, what = "Oracle (true betas)",           est = oracle, se = oracle_se)],
  agg_vm[, .(K_S, what = "REML (log fit)",                est = reml,   se = reml_se)]
)
plot_vm[, what := factor(what, levels = unique(what))]

p_vm_path = file.path(FIGS, "P_VM.png")
p_vm = ggplot(plot_vm, aes(factor(K_S), est, colour = what)) +
  geom_hline(aes(yintercept = V_M_true), colour = "grey60", linetype = "dashed") +
  geom_pointrange(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se),
                  position = position_dodge(width = 0.4)) +
  labs(x = "K_SEL", y = expression(V[M]), colour = NULL,
       title = expression(paste("Predicted, true (oracle) and estimated ", V[M])),
       subtitle = "Mean over replicates, 95% CI; dashed line = true value from SLiM parameters") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(p_vm_path, p_vm)






# summ[, K_S := factor(K_S, levels = K_S_vals)]

agg <- summ[, .(
  pred      = s_pred[1],
  oracle    = mean(oracle_s_hat), oracle_se = sd(oracle_s_hat) / sqrt(.N),
  reml      = mean(s_hat_log),    reml_se   = sd(s_hat_log)    / sqrt(.N)
), by = K_S]

plot_dt <- rbind(
  agg[, .(K_S, what = "Oracle (true betas)",          est = oracle, se = oracle_se)],
  agg[, .(K_S, what = "REML (log fit)",               est = reml,   se = reml_se)]
)



plot_dt[, what := factor(what, levels = unique(what))]

p4 = ggplot(plot_dt, aes(as.numeric(K_S), est, colour = what)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se),
                  position = position_dodge(width = 0.0006)) +
  geom_line(aes(group = what), position = position_dodge(width = 0.0006), alpha = 0.5) +
  labs(x = "K_SEL", y = "slope (s)", colour = NULL,
       title = "True and estimated slope",
       subtitle = "Mean over replicates, 95% CI") +
  theme_bw() + theme(legend.position = "bottom")
p4_path = file.path(FIGS, "Predicted_true_oracle_estimated_slope.png")
ggsave(p4_path, p4)



p1 = summ %>%
  ggplot(aes(x=K_S, y=s_hat_log)) +
  ggbeeswarm::geom_quasirandom()

p1_path = file.path(FIGS, "P1_K_S_S_hat.png")
ggsave(p1_path, p1)

p2 = summ %>%
  ggplot(aes(x=K_S, y=sigma2_m_log)) +
  ggbeeswarm::geom_quasirandom()

p2_path = file.path(FIGS, "P1_K_S_sigma2_m_log.png")
ggsave(p2_path, p2)

p3 = summ %>%
  ggplot(aes(x=oracle_s_hat, y=s_hat_log)) +
  ggbeeswarm::geom_quasirandom() +
  facet_wrap(~K_S, scales='free')

p3_path = file.path(FIGS, "P1_oracle_s_hat_s_hat_log.png")
ggsave(p3_path, p3)



summ[, .(pred = s_pred[1],
         oracle = mean(oracle_s_hat), oracle_sd = sd(oracle_s_hat),
         reml_log = mean(s_hat_log), reml_nls = mean(s_hat_nls, na.rm = TRUE)),
     by = K_S]

ggplot(summ, aes(factor(K_S), oracle_s_hat)) +
  geom_boxplot() +
  geom_point(aes(y = s_pred), colour = "red", size = 3) +
  labs(x = "K_SEL", y = "oracle slope", title = "Oracle slope vs prediction (red)")

