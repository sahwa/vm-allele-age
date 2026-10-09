library(data.table)
library(purrr)
library(ggplot2)
f = glue::glue

DATA="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v4.0"
FIGS="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v4.0"
OUT = file.path(DATA, "collated")
dir.create(OUT, showWarnings = FALSE)
source("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/diagGREML.R")

S_vals = c(0, 0.002, 0.005, 0.01)
N_REP = 25
K = 8
n_sample = 2500
H2 = 0.1

## ---- fit_sel_reml: keep your existing function here, unchanged ---------------

## ---- pooled dataset for one K_SEL: seeds stacked as extra people -------------
build_data = function(S) {
  set.seed(1e4 * which(S_vals == S))                 # reproducible phenotype noise
  files <- file.path(DATA, f("sing_K{S}_seed{1:N_REP}.tsv"))
  sing  <- rbindlist(map(files, fread), idcol = "rep")            # carrier, age, beta
  stopifnot(max(sing$carrier) <= n_sample)
  sing[, person := carrier + (rep - 1L) * n_sample]  # each seed adds n_sample new people
  n_tot <- n_sample * N_REP

  breaks <- unique(quantile(sing$age, probs = seq(0, 1, length.out = K + 1), type = 1))
  sing[, bin := cut(age, breaks, include.lowest = TRUE, labels = FALSE)]
  K_eff <- length(breaks) - 1

  g   <- numeric(n_tot)
  agg <- rowsum(sing$beta, sing$person)
  g[as.integer(rownames(agg))] <- agg[, 1]
  Y <- g + rnorm(n_tot, sd = sqrt(var(g) * (1 - H2) / H2))

  C <- table(factor(sing$person, levels = seq_len(n_tot)),
             factor(sing$bin,    levels = seq_len(K_eff)))
  C <- matrix(C, nrow = n_tot)
  list(sing = sing, C = C, Y = Y, K_eff = K_eff,
       rep_of = rep(seq_len(N_REP), each = n_sample))
}

## ---- all four methods, using only the seeds in `keep` ------------------------
fit_all = function(d, keep = seq_len(N_REP)) {
  rows <- d$rep_of %in% keep
  sing <- d$sing[rep %in% keep]
  C <- d$C[rows, , drop = FALSE]; Y <- d$Y[rows]; n <- sum(rows)

  bins <- sing[, .(mid = mean(age), oracle_b2 = mean(beta^2)), keyby = bin]
  m    <- colMeans(C)
  A    <- setNames(lapply(seq_len(d$K_eff), function(k) C[, k] / m[k]),
                   paste0("bin", seq_len(d$K_eff)))
  npp  <- sing[age == 0, .N] / n                     # new mutations per person

  ## oracle: true mean beta^2 per bin
  m_orc <- lm(log(oracle_b2) ~ mid, data = bins)

  ## two-stage (log-linear) and NLS, both from the per-bin REML estimates
  fit <- suppressMessages(fit_diagGREML(y = Y, A = A, X = matrix(1, n, 1),
                                        constraint = FALSE, magic0316 = TRUE))
  idx <- match(names(A), fit$Vlistnames)
  bins[, `:=`(sigma2_b = as.numeric(fit$varcmp)[idx] / m,
              se_b     = sqrt(diag(fit$Hi))[idx] / m)]
  pos   <- bins[is.finite(sigma2_b) & sigma2_b > 0]
  m_log <- if (nrow(pos) >= 3) lm(log(sigma2_b) ~ mid, data = pos) else NULL
  m_nls <- if (is.null(m_log)) NULL else tryCatch(
    nls(sigma2_b ~ a * exp(-s_par * mid), data = bins, weights = 1 / se_b^2,
        start = list(a = exp(coef(m_log)[[1]]), s_par = max(-coef(m_log)[[2]], 1e-4))),
    error = function(e) NULL)

  ## joint
  jf <- tryCatch(fit_sel_reml(y = Y, C = C, t = bins$mid), error = function(e) NULL)

  s2m <- c(exp(coef(m_orc)[[1]]),
           if (is.null(m_log)) NA_real_ else exp(coef(m_log)[[1]]),
           if (is.null(m_nls)) NA_real_ else coef(m_nls)[["a"]],
           if (is.null(jf))    NA_real_ else jf$sigma2_m)
  s   <- c(-coef(m_orc)[[2]],
           if (is.null(m_log)) NA_real_ else -coef(m_log)[[2]],
           if (is.null(m_nls)) NA_real_ else coef(m_nls)[["s_par"]],
           if (is.null(jf))    NA_real_ else jf$s)
  data.table(method = c("oracle", "two-stage", "nls", "joint"),
             s_hat = s, sigma2_m = s2m, V_M = s2m * npp, n_people = n)
}

## ---- full-data estimate plus leave-one-seed-out jackknife SE -----------------
jk_se = function(x) {
  x <- x[is.finite(x)]; n <- length(x)
  if (n < 2) NA_real_ else sqrt((n - 1) / n * sum((x - mean(x))^2))
}

run_K = function(S) {
  d    <- build_data(S)
  full <- fit_all(d)
  jk   <- rbindlist(lapply(seq_len(N_REP), function(i)
            fit_all(d, keep = setdiff(seq_len(N_REP), i))))
  se   <- jk[, .(V_M_se = jk_se(V_M), sigma2_m_se = jk_se(sigma2_m), s_se = jk_se(s_hat),
                 n_jk = sum(is.finite(V_M))), by = method]
  cbind(K_S = S, merge(full, se, by = "method"))
}

res <- rbindlist(lapply(S_vals, function(S) tryCatch(run_K(S), error = function(e) {
  message(f("K_S = {S} failed: {conditionMessage(e)}")); NULL })))
res[, method := factor(method, levels = c("oracle", "two-stage", "nls", "joint"))]
setorder(res, K_S, method)
fwrite(res, file.path(OUT, "pooled_estimates.csv"))
print(res)

## ---- plot V_M with 95% jackknife intervals ----------------------------------
cols <- c(oracle = "#D7191C", "two-stage" = "grey25", nls = "#E69F00", joint = "#2C7BB6")

p_VM <- ggplot(res[is.finite(V_M)], aes(method, V_M, colour = method)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_hline(data = res[method == "oracle"], aes(yintercept = V_M),
             colour = "#D7191C", linetype = "dashed") +
  geom_errorbar(aes(ymin = V_M - 1.96 * V_M_se, ymax = V_M + 1.96 * V_M_se),
                width = 0.2, na.rm = TRUE) +
  geom_point(size = 2.5) +
  facet_wrap(~ K_S, nrow = 1,
             labeller = labeller(K_S = function(x) paste0("K_SEL = ", x))) +
  scale_colour_manual(values = cols, guide = "none") +
  labs(x = NULL, y = expression(hat(V)[M]),
       title = "Mutational variance from pooled seeds",
       subtitle = "Points = estimate from all seeds; bars = 95% jackknife interval; dashed red = oracle") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())

ggsave(file.path(FIGS, "V_M_pooled_jackknife.png"), p_VM, width = 10, height = 4.5, dpi = 150)
