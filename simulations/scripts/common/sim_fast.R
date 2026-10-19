## ---------------------------------------------------------------------------
## Fast joint simulation of true and misclassified singleton counts.
##
## Replaces  sim_unbinned() + misclassify()  with a single pass.
##
## Because counts are Poisson, the true x observed table thins exactly:
##     M[i,f,k] ~ Pois(lambda_i * p_f * P[f,k]),  independent over (f,k)
## with  C[i,f] = sum_k M[i,f,k]  and  Cobs[i,k] = sum_f M[i,f,k].
## So no conditioning and no sequential binomial chain: one vectorised rpois
## per true bin instead of ncol(C) * |support| sequential rbinom calls.
##
## ~8000 full-length rbinom  ->  ~200 vectorised rpois.
## ---------------------------------------------------------------------------

## P = NULL gives the error-free case. tol trims negligible kernel support,
## which is where most of the remaining cost lives -- 1e-6 keeps >99.9% of each
## row's mass and is far below any bias the kernel itself induces.
sim_counts <- function(REP, pmf, S, P = NULL, tol = 1e-6,
                       n_inds = N_inds, n_sing = N_SING,
                       cv_shared = CV_SHARED, curve = CURVE) {
  set.seed(REP)                                       # common random numbers
  shape <- 1 / cv_shared^2
  lam_i <- rgamma(n_inds, shape, rate = shape / n_sing)

  n_t  <- length(pmf$t)
  fd   <- fdecay(S, pmf$t, curve)
  Cobs <- matrix(0L, n_inds, n_t)
  w    <- numeric(n_inds)

  for (f in seq_len(n_t)) {
    if (pmf$p[f] <= 0) next

    if (is.null(P)) {                                 # no age error
      sup <- f; pr <- 1
    } else {
      pr  <- P[f, ]
      sup <- which(pr > tol)
      if (!length(sup)) next
      pr  <- pr[sup] / sum(pr[sup])
    }

    rate <- outer(lam_i * pmf$p[f], pr)               # n_inds x |sup|
    M    <- matrix(rpois(length(rate), rate), nrow = n_inds)

    w <- w + fd[f] * rowSums(M)                       # true ages drive the truth
    Cobs[, sup] <- Cobs[, sup] + M                    # observed ages drive the fit
  }

  w
  list(C = Cobs, w = w, t = pmf$t)
}

## Phenotype from the TRUE weight, fitted against the OBSERVED counts.
make_pheno <- function(sim, h2 = H2, sigma_beta = SIGMA_BETA) {
  g <- sqrt(sim$w * sigma_beta^2) * rnorm(length(sim$w))
  g + rnorm(length(g), 0, sqrt(var(g) * (1 - h2) / h2))
}

## ---------------------------------------------------------------------------
## Standalone version, for a C you already have and whose counts may not be
## Poisson. The multinomial chain is then unavoidable, but two things help a
## lot: only touch people who still have singletons left to place (the set
## shrinks fast), and drop out once nothing remains.
## ---------------------------------------------------------------------------
misclassify <- function(C, P, tol = 1e-8) {
  n    <- nrow(C)
  Cobs <- matrix(0L, n, ncol(C))

  for (f in seq_len(ncol(C))) {
    sup <- which(P[f, ] > tol)
    if (!length(sup)) next
    ord <- sup[order(P[f, sup], decreasing = TRUE)]
    pr  <- P[f, ord] / sum(P[f, ord])

    idx <- which(C[, f] > 0L)                         # <- the win
    if (!length(idx)) next
    rem    <- C[idx, f]
    p_left <- 1

    for (k in seq_along(ord)) {
      x <- if (k == length(ord)) rem
           else rbinom(length(rem), rem, min(1, pr[k] / p_left))
      Cobs[idx, ord[k]] <- Cobs[idx, ord[k]] + x
      rem    <- rem - x
      p_left <- p_left - pr[k]

      keep <- rem > 0L                                # shrink as we go
      if (!any(keep)) break
      if (!all(keep)) { idx <- idx[keep]; rem <- rem[keep] }
    }
  }
  Cobs
}
