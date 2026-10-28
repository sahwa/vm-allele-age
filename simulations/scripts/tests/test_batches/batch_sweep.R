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

library(data.table); library(purrr); library(glue); library(RcppCNPy)

DATA <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/tests/batch_test"

MU       <- 1.25e-8
L_GENOME <- 3e9
L_SIM    <- 1e8                      # sequence length the tree sequence covers
U        <- 2 * MU * L_GENOME        # ~75 mutations per person per generation
C_VM     <- 0.3                      # lambda ~ C_VM * (V_M/V_P) / h2, from SLiM

N_BATCHES   <- c(1, 2, 5, 10, 20, 50, 100, 200)
N_TOTAL     <- 5e5
LAMBDA_VALS <- c(0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2)
VM_VP       <- 2.6e-3
Z_ALPHA     <- 1.96                  # two-sided 0.05

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
prepare_counts <- function(C, trep, kdist, scale = L_GENOME / L_SIM) {
  C  <- as.matrix(C)
  ok <- is.finite(trep) & colSums(C) > 0
  list(C = C[, ok, drop = FALSE] * scale, ages = trep[ok], kdist = kdist)
}

fisher_se_lambda_trait <- function(prep, vm_vp, lambda, curve = "hyperbolic") {
  C <- prep$C; ages <- prep$ages
  sigma_m2 <- vm_vp / U                       # V_P standardised to 1
  w  <- if (curve == "exponential") exp(-lambda * ages) else 1 / (1 + lambda * ages)

  h2_sing <- sigma_m2 * mean(as.vector(C %*% w))
  t0      <- sum(colSums(C) * ages) / sum(colSums(C))

  ## h2_sing >= 1 is not a numerical failure -- it is the design breaking down.
  ## Past this point "batch-singletons" are no longer recent mutations and the
  ## model is incoherent with V_P = 1. Record it rather than stopping, so the
  ## sweep completes and shows where the boundary falls.
  if (!is.finite(h2_sing) || h2_sing >= 1)
    return(list(se = NA_real_, h2_sing = h2_sing, t0 = t0))

  list(se = fisher_se_lambda(C, ages, sigma_m2, 1 - h2_sing, lambda, curve),
       h2_sing = h2_sing, t0 = t0)
}

power_lambda <- function(lambda, se, z_a = Z_ALPHA) {
  ifelse(is.na(se), NA_real_, pnorm((lambda / se) - z_a))
}

decay <- function(lambda, t, curve = "hyerpbolic") {
  switch(curve,
    exponential = exp(-lambda * t),
    hyperbolic  = 1 / (1 + lambda * t),
    stop("curve must be 'exponential' or 'hyperbolic'"))
}

## ---------------------------------------------------------------------------
## sweep
## ---------------------------------------------------------------------------
res <- rbindlist(map(N_BATCHES, function(nb) {
  prep <- prepare_counts(
    npyLoad(file.path(DATA, glue("C_nb{nb}.npy"))),
    npyLoad(file.path(DATA, glue("trep_nb{nb}.npy"))),
    npyLoad(file.path(DATA, glue("kdist_nb{nb}.npy")))
    )
  rbindlist(map(LAMBDA_VALS, function(lam) {

    r <- fisher_se_lambda_trait(prep, vm_vp = VM_VP, lambda = lam)
    data.table(
      n_batches  = nb,
      batch_size = N_TOTAL / nb,
      lambda     = lam,
      ## lambda = C_VM*(V_M/V_P)/h2, so a given (lambda, V_M/V_P) implies a
      ## heritability. Above 1 the combination is not realisable: there is a
      ## floor lambda >= C_VM*(V_M/V_P) = 7.8e-4 at V_M/V_P = 2.6e-3.
      h2_implied = if (lam > 0) C_VM * VM_VP / lam else NA_real_,
      h2_sing    = r$h2_sing,
      t0         = r$t0,
      se         = r$se,
      power      = power_lambda(lam, r$se),
      sigma_m_2 = VM_VP / U,
      kdist = prep$kdist
    )
  }))
}))

res[, realisable := is.na(h2_implied) | h2_implied <= 1]

## lambda = 0 is the false-positive rate (0.025, one tail of a two-sided test),
## not power -- kept as a check.
setorder(res, n_batches, lambda)
print(res[], digits = 3)

## The diagnostics to read BEFORE any power number: if h2_sing climbs linearly
## past 1 then nothing saturated and the T0_CAP question is still open; if
## h2_sing and t0 flatten together, that is the saturation point the design
## hinges on, measured rather than assumed.
cat("\n-- design diagnostics (independent of lambda) --\n")
print(unique(res[, .(n_batches, batch_size, t0,
                     h2_sing = round(h2_sing, 4))]), digits = 3)

fwrite(res, file.path(DATA, "fisher_power_by_batch.csv"))
