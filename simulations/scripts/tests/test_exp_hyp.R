f = glue::glue
library(data.table)
library(purrr)
library(parallel)

COMMON="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/common"

source(file.path(COMMON, "fit_sel_reml_joint.R"))
source(file.path(COMMON, "sim_diag.R"))
source("https://raw.githubusercontent.com/huanweiuq/diagGREML/refs/heads/main/diagGREML.R")

N_CORES = parallel::detectCores() - 6

Svals  = c(0, 1e-4, 1e-3, 1e-2, 3e-2)
H2     = 0.75
CV     = 0.65
N_inds = 5e4
K      = 30
t_max  = 400
LAMBDA = 1000
N_reps = 10
SIGMA_BETA = 1

U = LAMBDA * K / t_max        # new mutations per person per generation

params = CJ(REP = 1:N_reps, S = Svals)

run_once = function(REP, S) {
  message(f("S = {S}, REP = {REP}"))
  set.seed(1e4 * which(Svals == S) + REP)
  ## simulating under the hyperbolic model: SLiM says that is the better description
  ## of selection, so the exponential fit is the misspecified one here
  d = sim_diag(S = S, REP = REP, CV = CV, N_inds = N_inds, K = K, t_max = t_max,
               LAMBDA = LAMBDA, type = "hyperbolic", H2 = H2)
  list(
    S = S, REP = REP, y_var = var(d$Y),
    fit_exp = fit_sel_reml(y = d$Y, C = d$C, t = d$t_mid, curve = "exp"),
    fit_hyp = fit_sel_reml(y = d$Y, C = d$C, t = d$t_mid, curve = "hyp")
  )
}

res = mclapply(seq_len(nrow(params)), function(i) {
  tryCatch(run_once(REP = params$REP[i], S = params$S[i]), error = function(e) {
    message(f("S = {params$S[i]}, REP = {params$REP[i]} failed: {conditionMessage(e)}")); NULL })
}, mc.cores = N_CORES)

res_dt = rbindlist(map(compact(res), function(x) {
  data.table(
    S = x$S, REP = x$REP, y_var = x$y_var,
    h2m_true    = U * SIGMA_BETA^2 / x$y_var,          # truth for this replicate
    s_exp       = x$fit_exp$s,        s_hyp       = x$fit_hyp$s,
    s_exp_se    = x$fit_exp$s_se,     s_hyp_se    = x$fit_hyp$s_se,
    sigma_m_exp = x$fit_exp$sigma2_m, sigma_m_hyp = x$fit_hyp$sigma2_m,
    ll_exp      = x$fit_exp$logLik,   ll_hyp      = x$fit_hyp$logLik,
    bound_exp   = x$fit_exp$at_bound, bound_hyp   = x$fit_hyp$at_bound
  )
}))

## mutational heritability from each fit, and which curve the data preferred
res_dt[, `:=`(h2m_exp = U * sigma_m_exp / y_var,
              h2m_hyp = U * sigma_m_hyp / y_var,
              hyp_wins = ll_hyp > ll_exp)]

## long: one row per replicate per fitted curve
long = rbind(
  res_dt[, .(S, REP, curve = "exponential", s_hat = s_exp, s_se = s_exp_se,
             h2m = h2m_exp, h2m_true, at_bound = bound_exp, hyp_wins)],
  res_dt[, .(S, REP, curve = "hyperbolic",  s_hat = s_hyp, s_se = s_hyp_se,
             h2m = h2m_hyp, h2m_true, at_bound = bound_hyp, hyp_wins)])
long[, curve := factor(curve, levels = c("hyperbolic", "exponential"))]

## bias and RMSE: each estimate against its own truth
summ = long[, {
  zh <- h2m[is.finite(h2m)]   - h2m_true[is.finite(h2m)]
  zs <- s_hat[is.finite(s_hat)] - S
  .(n_fit        = sum(is.finite(s_hat)),
    n_bound      = sum(at_bound, na.rm = TRUE),
    pct_hyp_wins = 100 * mean(hyp_wins, na.rm = TRUE),
    h2m_true     = mean(h2m_true),
    h2m_median   = median(h2m, na.rm = TRUE),
    h2m_bias     = mean(zh),
    h2m_bias_pct = 100 * mean(zh) / mean(h2m_true),
    h2m_rmse     = sqrt(mean(zh^2)),
    s_median     = median(s_hat, na.rm = TRUE),
    s_bias       = mean(zs),
    s_rmse       = sqrt(mean(zs^2)))
}, keyby = .(S, curve)]

print(summ[, .(S, curve, n_fit, n_bound, pct_hyp_wins,
               h2m_true, h2m_median, h2m_bias_pct, h2m_rmse,
               s_median, s_bias, s_rmse)])
