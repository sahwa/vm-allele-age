## ---------------------------------------------------------------------------
## batch_power_fisher.R -- power for lambda across the batching ladder, using
## the REAL count matrices from batch_sweep.py rather than a scaling formula.
##
## Batching is encoded entirely in C and in the representative ages, so the
## Fisher machinery needs no re-parameterisation -- only sigma2_m must be held
## fixed (it is a property of the trait) while h2_sing is allowed to emerge
## from whichever C is fed in. Holding h2_sing fixed would cancel the main
## benefit of batching.
##
## NOT included here: cross-batch variant sharing. Fisher information assumes
## independent observations, so these are NAIVE standard errors -- an upper
## bound on power. The sharing penalty has to come from kdist_nb*.npy, or
## better, from simulating phenotypes off these same C matrices and fitting.
## ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## Fisher SE for lambda, profiling out sigma2_e and sigma2_m.
##   w = f(lambda, t);  dw/dlambda = -t*w (exponential), -t*w^2 (hyperbolic)
## Validated against a numerical Hessian of the expected log-likelihood.
## ---------------------------------------------------------------------------
fisher_se_lambda <- function(C, ages, sigma_m2, sigma_e2, lambda,
                             curve = c("hyperbolic", "exponential")) {
  curve <- match.arg(curve)
  C <- as.matrix(C)
  stopifnot(ncol(C) == length(ages), sigma_m2 > 0, sigma_e2 > 0,
            all(C >= 0), all(is.finite(C)), all(is.finite(ages)))

  w  <- if (curve == "exponential") exp(-lambda * ages) else 1 / (1 + lambda * ages)
  dw <- if (curve == "exponential") ages * w else ages * w^2

  S <- as.vector(C %*% w)
  Q <- as.vector(C %*% dw)
  v <- sigma_e2 + sigma_m2 * S
  if (any(v <= 0)) stop("Non-positive phenotypic variance")

  G <- cbind(sigma_e2 = 1, sigma_m2 = S, lambda = -sigma_m2 * Q)
  I <- crossprod(G / (sqrt(2) * v))
  V <- tryCatch(solve(I), error = function(e) NULL)
  if (is.null(V)) { warning("singular information matrix"); return(NA_real_) }
  sqrt(V["lambda", "lambda"])
}

## ---------------------------------------------------------------------------
## Prepare one level's counts once: drop empty bins, rescale to genome-wide.
## Hoisted out of the lambda loop so a 48MB matrix is not re-copied per lambda.
##
## The tree sequence covers L_SIM, so raw counts are ~1/30 of genome-wide.
## Scaling a realised count preserves CV_total exactly where it should fall
## slightly (0.38 -> 0.34) as Poisson noise shrinks, so this is mildly
## conservative.
## ---------------------------------------------------------------------------
prepare_counts <- function(C, trep, scale = L_GENOME / L_SIM) {
  C  <- as.matrix(C)
  ok <- is.finite(trep) & colSums(C) > 0
  list(C = C[, ok, drop = FALSE] * scale, ages = trep[ok])
}

fisher_se_lambda_trait <- function(prep, vm_vp, lambda, curve = "hyperbolic",
                                   N_BB = 1, N_PER = 5e5) {
  C <- prep$C; ages <- prep$ages
  sigma_m2 <- vm_vp / U                       # V_P standardised to 1
  w  <- if (curve == "exponential") exp(-lambda * ages) else 1 / (1 + lambda * ages)

  ## per-cohort quantities -- independent of N_BB by construction
  h2_sing <- sigma_m2 * mean(as.vector(C %*% w))
  t0      <- sum(colSums(C) * ages) / sum(colSums(C))

  N_total <- N_BB * N_PER

  ## h2_sing >= 1 is not a numerical failure -- it is the design breaking down.
  ## Past this point "batch-singletons" are no longer recent mutations and the
  ## model is incoherent with V_P = 1. Record it rather than stopping, so the
  ## sweep completes and shows where the boundary falls.
  if (!is.finite(h2_sing) || h2_sing >= 1)
    return(list(se = NA_real_, se_sim = NA_real_, h2_sing = h2_sing,
                t0 = t0, N_total = N_total))

  se_sim <- fisher_se_lambda(C, ages, sigma_m2, 1 - h2_sing, lambda, curve)

  list(se      = se_sim * sqrt(nrow(C) / N_total),   # scaled to the real design
       se_sim  = se_sim,                             # raw, for QC
       h2_sing = h2_sing,
       t0      = t0,
       N_total = N_total)
}


power_lambda <- function(lambda, se, z_a = Z_ALPHA) {
  ifelse(is.na(se), NA_real_, pnorm((lambda / se) - z_a))
}

decay <- function(lambda, t, curve = "hyperbolic") {
  switch(curve,
    exponential = exp(-lambda * t),
    hyperbolic  = 1 / (1 + lambda * t),
    stop("curve must be 'exponential' or 'hyperbolic'"))
}


n_biobanks_for_power <- function(prep, vm_vp, lambda, curve = "hyperbolic",
                                 power = 0.8, z_a = Z_ALPHA, N_PER = 5e5) {
  fit <- fisher_se_lambda_trait(prep, vm_vp, lambda, curve,
                                N_BB = 1, N_PER = N_PER)

  if (!is.finite(fit$se_sim) || lambda <= 0)
    return(list(N_BB = NA_real_, N_total = NA_real_, se_need = NA_real_,
                h2_sing = fit$h2_sing, t0 = fit$t0))

  se_need <- lambda / (z_a + qnorm(power))
  N_total <- nrow(prep$C) * (fit$se_sim / se_need)^2

  list(N_BB    = N_total / N_PER,
       N_total = N_total,
       se_need = se_need,
       h2_sing = fit$h2_sing,
       t0      = fit$t0)
}

