library(data.table)
library(purrr)
library(future)
library(furrr)
library(patchwork)
f = glue::glue

DATA="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v4.0"
FIGS="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v4.0"
OUT = file.path(DATA, "collated")
dir.create(OUT, showWarnings = FALSE)
source("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/diagGREML.R")

S_vals = c(0, 0.002, 0.005, 0.01)
N_REPS = 25
K = 8
n_sample = 2500
H2 = 0.1

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
  list(A = A, C = C, Y = Y, sing = sing, t = bins$mid, bins = bins, m = m, K_eff = K_eff, y_var = var(Y))
}

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

  # ---- Second stage: log-linear (drops negatives) and NLS (keeps them) -----
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

## ---- one (K_SEL, seed): both fits, returned as one row plus a bins table ----
run_once = function(S, REP) {
  d  <- get_counts(S, K, H2, REP)
  dg <- run_diagGREML(d, K_S = S, REP = REP, H2 = H2)
  jf <- tryCatch(fit_sel_reml(y = d$Y, C = d$C, t = d$t), error = function(e) NULL)

  gj <- function(x) if (is.null(jf)) NA_real_ else as.numeric(jf[[x]])
  joint <- data.table(s_hat_joint = gj("s"), s_se_joint = gj("s_se"),
                      sigma2_m_joint = gj("sigma2_m"),
                      sigma2_t0_joint = gj("sigma2_t0"), t0_joint = gj("t0"),
                      sigma2_e_joint = gj("sigma2_e"), lrt_p_joint = gj("lrt_p"),
                      at_bound_joint = gj("at_bound"), conv_joint = gj("convergence"))
  joint[, V_M_joint := sigma2_m_joint * dg$summary$new_per_person]

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
  reps[, .(K_S, REP, method = "oracle",    s_hat = oracle_s_hat, sigma2_m = oracle_sigma2_m, V_M = V_M_oracle,                   V_P = y_var)],
  reps[, .(K_S, REP, method = "two-stage", s_hat = s_hat_log,    sigma2_m = sigma2_m_log,    V_M = V_M_reml,                     V_P = y_var)],
  reps[, .(K_S, REP, method = "nls",       s_hat = s_hat_nls,    sigma2_m = sigma2_m_nls,    V_M = sigma2_m_nls * new_per_person, V_P = y_var)],
  reps[, .(K_S, REP, method = "joint",     s_hat = s_hat_joint,  sigma2_m = sigma2_m_joint,  V_M = V_M_joint,                    V_P = y_var)])
long[, method := factor(method, levels = c("oracle", "two-stage", "nls", "joint"))]
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

### plot V_M / V_P

est   <- long[method != "oracle"]
truth <- long[method == "oracle", .(h2m = median(h2m, na.rm = TRUE)), by = K_S]
med   <- est[, .(h2m = median(h2m, na.rm = TRUE)), by = .(K_S, method)]
cols  <- c("two-stage" = "grey25", "nls" = "#E69F00", "joint" = "#2C7BB6")
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
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())

ggsave(file.path(FIGS, "h2m_single_joint_oracle.png"), p1_h2m, width = 10, height = 4.5, dpi = 150)


#### bias and RMSE ####

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
        legend.position = "bottom")

ggsave(file.path(FIGS, "bias_rmse_by_method.png"), p_perf, width = 10, height = 6, dpi = 150)


fwrite(reps, file.path(OUT, "reps_wide.csv"))
fwrite(long, file.path(OUT, "reps_long.csv"))
fwrite(bins, file.path(OUT, "bins.csv"))
fwrite(summ, file.path(OUT, "summary.csv"))

print(summ)
