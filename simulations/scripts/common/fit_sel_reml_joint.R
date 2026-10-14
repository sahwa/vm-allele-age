fit_sel_reml <- function(y, C, t, Z = matrix(1, length(y), 1),
                         s_max = 0.1,
                         curve = c("exp", "hyp")) {

  curve <- match.arg(curve)
  storage.mode(C) <- "double"

  if (length(y) != nrow(C)) stop("length(y) must equal nrow(C)")
  if (length(t) != ncol(C)) stop("length(t) must equal ncol(C)")
  if (nrow(Z) != length(y)) stop("nrow(Z) must equal length(y)")

  cm <- colMeans(C)
  if (sum(cm) <= 0) stop("No mutations in C")
  t0 <- sum(cm * t) / sum(cm)
  vy <- var(y)

  ## Starting values: average genetic variance, residual variance
  theta0 <- c(0.05 * vy, 0.95 * vy)

  ## Per-mutation variance at each age relative to age t0
  rel <- switch(
    curve,
    exp = function(s) exp(-s * (t - t0)),
    hyp = function(s) (1 + s * t0) / (1 + s * t)
  )

  ## Convert variance at t0 back to variance at age zero
  to_zero <- switch(
    curve,
    exp = function(s) exp(s * t0),
    hyp = function(s) 1 + s * t0
  )

  ## Hyperbolic model has a pole at 1 + s*t = 0
  s_hi <- s_max
  s_lo <- if (curve == "exp") {
    -s_max
  } else {
    max(-s_max, -0.9 / max(t))
  }

  ## REML fit of (a, sigma2_e) for a fixed value of s
  inner <- function(s, th = theta0) {
    rk <- rel(s)
    if (any(!is.finite(rk)) || any(rk <= 0)) {
      return(list(
        nll = Inf, th = th, b = NA_real_,
        wbar = NA_real_, converged = FALSE
      ))
    }

    ## Individual age-weighted mutation burden
    w <- as.vector(C %*% rk)
    wbar <- mean(w)

    if (!is.finite(wbar) || wbar <= 0) {
      return(list(
        nll = Inf, th = th, b = NA_real_,
        wbar = NA_real_, converged = FALSE
      ))
    }

    ## Scale so that mean(w) = 1
    w <- w / wbar

    obj <- function(th) {

      ## Both variance components must be non-negative
      if (th[1] < 0 || th[2] <= 0)
        return(NULL)

      v <- th[1] * w + th[2]

      if (any(!is.finite(v)) || min(v) <= 0)
        return(NULL)

      ZtViZ <- crossprod(Z, Z / v)
      detZ <- determinant(ZtViZ, logarithm = TRUE)
      if (detZ$sign <= 0)
        return(NULL)
      b <- tryCatch(
        solve(ZtViZ, crossprod(Z, y / v)),
        error = function(e) NULL
      )
      if (is.null(b))
        return(NULL)
      r <- y - as.vector(Z %*% b)
      list(
        v = v,
        r = r,
        b = as.vector(b),
        ZtViZ = ZtViZ,
        nll = 0.5 * (
          sum(log(v)) +
          sum(r^2 / v) +
          as.numeric(detZ$modulus)
        )
      )
    }

    cur <- obj(th)

    if (is.null(cur)) {
      th <- c(0, vy)
      cur <- obj(th)
    }

    if (is.null(cur)) {
      return(list(
        nll = Inf, th = th, b = NA_real_,
        wbar = wbar, converged = FALSE
      ))
    }

    converged <- FALSE

    for (it in 1:50) {

      v <- cur$v
      inv_ZtViZ <- tryCatch(
        solve(cur$ZtViZ),
        error = function(e) NULL
      )
      if (is.null(inv_ZtViZ))
        break
      ## REML leverage correction
      hlev <- rowSums((Z %*% inv_ZtViZ) * Z)
      q <- cur$r^2 / v^2 - (1 - hlev / v) / v

      score <- 0.5 * c(
        sum(q * w),
        sum(q)
      )

      info <- 0.5 * matrix(
        c(
          sum(w^2 / v^2),
          sum(w / v^2),
          sum(w / v^2),
          sum(1 / v^2)
        ),
        2
      )

      step <- tryCatch(
        solve(info, score),
        error = function(e) NULL
      )

      if (is.null(step) || any(!is.finite(step)))
        break

      new <- NULL

      ## Step-halving if proposed step is invalid or makes fit worse
      for (half in 1:25) {

        candidate <- th + step
        new <- obj(candidate)

        if (!is.null(new) &&
            new$nll <= cur$nll + 1e-8)
          break

        new <- NULL
        step <- step / 2
      }

      if (is.null(new))
        break

      th <- th + step
      cur <- new

      if (max(abs(step)) < 1e-8 * vy) {
        converged <- TRUE
        break
      }
    }

    list(
      nll = cur$nll,
      th = th,
      b = cur$b,
      wbar = wbar,
      converged = converged
    )
  }

  ## IMPORTANT: deterministic profile.
  ## Every value of s starts from the same theta0.
  prof <- function(s) {
    inner(s, theta0)$nll
  }

  ## Coarse search, dense near zero
  g <- 10^seq(-3, 0, length.out = 13)

  grid <- sort(unique(c(
    s_lo * g,
    0,
    s_hi * g
  )))

  pg <- vapply(grid, prof, numeric(1))

  j <- which.min(pg)

  lo <- grid[max(j - 1, 1)]
  hi <- grid[min(j + 1, length(grid))]

  opt_s <- optimize(
    prof,
    interval = c(lo, hi),
    tol = 1e-6
  )

  s_hat <- opt_s$minimum
  f <- inner(s_hat, theta0)

  at_bound <-
    abs(s_hat - s_lo) < 1e-3 * max(1, abs(s_lo)) ||
    abs(s_hat - s_hi) < 1e-3 * max(1, abs(s_hi))

  ## Profile curvature SE
  h <- min(
    2e-4,
    0.5 * (s_hat - s_lo),
    0.5 * (s_hi - s_hat)
  )

  curv <- if (h > 0) {
    (prof(s_hat + h) -
       2 * f$nll +
       prof(s_hat - h)) / h^2
  } else {
    NA_real_
  }

  s_se <- if (!at_bound &&
              is.finite(curv) &&
              curv > 0) {
    1 / sqrt(curv)
  } else {
    NA_real_
  }

  ## Explicit null fit at s = 0
  f0 <- inner(0, theta0)

  lrt <- if (is.finite(f0$nll) && is.finite(f$nll)) {
    max(0, 2 * (f0$nll - f$nll))
  } else {
    NA_real_
  }

  ## th[1] is average genetic variance because mean(w) = 1.
  ## Undo the normalization to get per-mutation variance at t0.
  a_t0 <- f$th[1] / f$wbar

  list(
    curve = curve,

    s = s_hat,
    s_se = s_se,

    sigma2_m = a_t0 * to_zero(s_hat),
    sigma2_t0 = a_t0,
    log_sigma2_m_se = NA_real_,

    sigma2_e = f$th[2],
    b = f$b,

    h2_sing = f$th[1] / (f$th[1] + f$th[2]),

    logLik = -f$nll,

    convergence = if (f$converged) 0L else 1L,

    lrt = lrt,
    lrt_p = if (is.finite(lrt))
      pchisq(lrt, 1, lower.tail = FALSE)
    else NA_real_,

    at_bound = at_bound,

    t0 = t0,
    s_lo = s_lo,
    s_hi = s_hi,

    ## Useful for diagnosing odd fits
    profile = data.frame(
      s = grid,
      nll = pg,
      delta2logLik = 2 * (pg - min(pg))
    )
  )
}
