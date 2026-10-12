fit_sel_reml <- function(y, C, t, Z = matrix(1, length(y), 1), s_max = 0.1,
                         curve = c("exp", "hyp")) {
  curve <- match.arg(curve)
  storage.mode(C) <- "double"
  cm <- colMeans(C)
  t0 <- sum(cm * t) / sum(cm)
  vy <- var(y)
  theta <- c(0.05 * vy, 0.95 * vy)              # (a, s2e); warm-started across values of s

  ## per-mutation variance at each bin age, relative to the variance at age t0
  rel <- switch(curve,
                exp = function(s) exp(-s * (t - t0)),
                hyp = function(s) (1 + s * t0) / (1 + s * t))
  ## multiplier taking the variance at age t0 back to age 0
  to_zero <- switch(curve,
                    exp = function(s) exp(s * t0),
                    hyp = function(s) 1 + s * t0)
  ## search range for s: the hyperbolic curve needs 1 + s * t > 0 in every bin
  s_hi <- s_max
  s_lo <- if (curve == "exp") -s_max else max(-s_max, -0.9 / max(t))

  ## REML fit of (a, s2e) for one fixed s, by Fisher scoring with step-halving
  inner <- function(s, th) {
    rk <- rel(s)
    if (any(!is.finite(rk)) || any(rk <= 0))
      return(list(nll = Inf, th = th, b = NA_real_, wbar = NA_real_, converged = FALSE))
    w <- as.vector(C %*% rk); wbar <- mean(w); w <- w / wbar
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

  prof <- function(s) {
    f <- inner(s, theta)
    if (is.finite(f$nll)) theta <<- f$th
    f$nll
  }

  ## coarse search on a grid that is dense near zero, then refine around the best point
  g    <- 10^seq(-3, 0, length.out = 13)
  grid <- c(rev(s_lo * g), 0, s_hi * g)
  pg   <- vapply(grid, prof, 0)
  j    <- which.min(pg)
  s_hat <- optimize(prof, c(grid[max(j - 1, 1)], grid[min(j + 1, length(grid))]),
                    tol = 1e-6)$minimum
  f <- inner(s_hat, theta)
  at_bound <- s_hat > s_hi * (1 - 1e-3) || s_hat < s_lo * (1 - 1e-3)

  ## SE of s from the curvature of the profile; likelihood-ratio test of s = 0
  h    <- min(2e-4, 0.5 * (s_hat - s_lo), 0.5 * (s_hi - s_hat))
  curv <- if (h > 0) (prof(s_hat + h) - 2 * f$nll + prof(s_hat - h)) / h^2 else NA_real_
  s_se <- if (!at_bound && is.finite(curv) && curv > 0) 1 / sqrt(curv) else NA_real_
  lrt  <- max(0, 2 * (pg[grid == 0] - f$nll))

  a_t0 <- f$th[1] / f$wbar                       # per-mutation variance at age t0
  list(curve = curve, s = s_hat, s_se = s_se,
       sigma2_m = a_t0 * to_zero(s_hat), log_sigma2_m_se = NA_real_,
       sigma2_e = f$th[2], b = f$b,
       logLik = -f$nll, convergence = if (f$converged) 0L else 1L,
       lrt = lrt, lrt_p = pchisq(lrt, 1, lower.tail = FALSE),
       at_bound = at_bound, sigma2_t0 = a_t0, t0 = t0, s_lo = s_lo, s_hi = s_hi)
}
