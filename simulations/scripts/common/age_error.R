## ---------------------------------------------------------------------------
## Allele-age error, split into two independently switchable components.
##
##   noise : unbiased random scatter.  E[t_obs | t] = t exactly.
##           Classical error-in-variables -> attenuates s_hat toward zero.
##
##   bias  : deterministic distortion of the age scale, applied BEFORE the
##           scatter.  Two knobs, because ARG-based dating fails in two
##           distinguishable ways:
##             scale : all ages multiplied by a constant
##             gamma : compression toward t_ref on the log scale, i.e.
##                     regression to the mean.  Young variants dated too old,
##                     old variants too young.  gamma < 1 shrinks the range.
##
## Setting cv = 0 gives pure bias; scale = 1, gamma = 1 gives pure noise.
## The two together are a clean 2x2 for attributing the damage.
## ---------------------------------------------------------------------------

## Deterministic map: true age -> MEAN observed age.
## t_ref is the pivot that compression leaves unmoved; use the count-weighted
## mean age so gamma rotates the age scale about the centre of mass of the data.
make_age_bias <- function(scale = 1, gamma = 1, t_ref = 30) {
  force(scale); force(gamma); force(t_ref)
  function(t) t_ref * scale * (t / t_ref)^gamma
}

## Row k of P = distribution of OBSERVED age bin given TRUE age t[k].
## cv is the coefficient of variation of the scatter, mean-preserving about
## bias(t[k]).  cv = 0 collapses to a deterministic relabelling.
make_age_error_matrix <- function(t, cv = 0, bias = make_age_bias()) {
  n  <- length(t)
  mu <- pmax(bias(t), 1e-8)
  P  <- matrix(0, n, n)

  if (cv <= 0) {
    k <- pmin(pmax(1L, as.integer(round(mu))), n)
    P[cbind(seq_len(n), k)] <- 1
    return(P)
  }

  sdlog   <- sqrt(log(1 + cv^2))
  meanlog <- log(mu) - sdlog^2 / 2          # => E[t_obs] = mu
  edges   <- c(0, t + 0.5)
  P <- t(vapply(seq_len(n),
                function(k) diff(plnorm(edges, meanlog[k], sdlog)),
                numeric(n)))
  P / pmax(rowSums(P), .Machine$double.eps)
}

## Diagnostics: what the error actually did, before any model is fitted.
## - realised mean observed age per true age (checks noise is unbiased)
## - E[true age | observed bin], which is the Berkson / regression-calibration
##   target.  Feeding THIS back as the bin age is the attenuation fix.
age_error_report <- function(t, P, p_true) {
  joint   <- p_true * P                      # P(true = k, obs = j)
  p_obs   <- colSums(joint)
  e_obs_g <- as.vector(P %*% t)              # E[t_obs | t_true]
  e_true_g <- ifelse(p_obs > 0,
                     as.vector(t %*% joint) / pmax(p_obs, 1e-300), NA_real_)
  data.table::data.table(
    t_true      = t,
    p_true      = p_true,
    E_obs_given_true  = e_obs_g,
    t_obs       = t,
    p_obs       = p_obs,
    E_true_given_obs  = e_true_g             # <- use as calibrated bin age
  )
}

## ---------------------------------------------------------------------------
## Conditions.  Run in this order: none -> noise -> bias -> both, so each
## column of the result isolates one mechanism.
## ---------------------------------------------------------------------------
ERR_GRID <- list(
  none       = list(cv = 0.0, scale = 1.0, gamma = 1.0),
  noise_lo   = list(cv = 0.3, scale = 1.0, gamma = 1.0),
  noise_hi   = list(cv = 0.7, scale = 1.0, gamma = 1.0),
  bias_scale = list(cv = 0.0, scale = 0.8, gamma = 1.0),
  bias_comp  = list(cv = 0.0, scale = 1.0, gamma = 0.7),
  both       = list(cv = 0.5, scale = 0.8, gamma = 0.7)
)

build_err_matrices <- function(t, t_ref, grid = ERR_GRID) {
  lapply(grid, function(g)
    make_age_error_matrix(t, cv = g$cv,
                          bias = make_age_bias(g$scale, g$gamma, t_ref)))
}
