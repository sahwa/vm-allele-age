## ---------------------------------------------------------------------------
## Power to detect lambda with N_BB independent biobanks of N_PER people each.
##
## TWO DISTINCT AXES -- do not confuse them:
##   nb    = batches WITHIN one cohort (the C_nb{nb}.npy files). Pays the
##           sqrt(k) variant-sharing penalty AND the fixed-allele-count
##           ascertainment bias that flips the sign of lambda-hat. Keep nb = 1.
##   N_BB  = independent biobanks, each defining singletons within its own
##           cohort. Different people carry different recent mutations, so
##           Fisher information adds cleanly and h2_sing / t0 are unchanged.
##
## Stacking the two (nb > 1 AND N_BB > 1) would multiply the batching
## pathologies by the number of banks. Don't.
##
## Requires power.R with the N_BB-aware fisher_se_lambda_trait() and
## n_biobanks_for_power() (see power_additions.R).
##
## Genome scaling is already handled inside prepare_counts() via
## scale = L_GENOME / L_SIM, so there is nothing to apply here. Note that
## default is lazily evaluated against the globals below -- L_GENOME and
## L_SIM must be defined before any prepare_counts() call.
## ---------------------------------------------------------------------------

library(data.table); library(purrr); library(glue); library(RcppCNPy)

source("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/common/power.R")

DATA <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/tests/batch_test"

MU       <- 1.25e-8
L_GENOME <- 3e9
L_SIM    <- 1e8                      # sequence length the tree sequence covers
U        <- 2 * MU * L_GENOME        # ~75 mutations per person per generation
C_VM     <- 0.3                      # lambda ~ C_VM * (V_M/V_P) / h2, from SLiM

N_PER       <- 5e5                   # people per biobank
N_BBs       <- c(1, 2, 3, 4, 5, 7, 10, 20, 50, 100)
LAMBDA_VALS <- c(0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2)
VM_VP       <- 2.6e-3
Z_ALPHA     <- 1.96

NB <- 1                              # no within-cohort batching

## --- load once, outside the sweep ------------------------------------------
prep <- prepare_counts(
  npyLoad(file.path(DATA, glue("C_nb{NB}.npy"))),
  npyLoad(file.path(DATA, glue("trep_nb{NB}.npy")))
)

cat(glue("genome-wide singletons per person in C: ",
         "{signif(mean(rowSums(prep$C)), 4)}"), "\n")
cat(glue("sigma_m2 = {signif(VM_VP / U, 4)}   ",
         "implied h2_sing ~ {signif(VM_VP / U * mean(rowSums(prep$C)), 4)}"),
    "\n\n")

## --- power over the (lambda, N_BB) grid ------------------------------------
grid <- CJ(lambda = LAMBDA_VALS, N_BB = N_BBs)

res <- grid[, {
  r <- fisher_se_lambda_trait(prep, vm_vp = VM_VP, lambda = lambda,
                              N_BB = N_BB, N_PER = N_PER)
  .(n_per_bb   = N_PER,
    N_total    = r$N_total,
    h2_implied = if (lambda > 0) C_VM * VM_VP / lambda else NA_real_,
    h2_sing    = r$h2_sing,
    t0         = r$t0,
    se         = r$se,
    power      = power_lambda(lambda, r$se))
}, by = .(lambda, N_BB)]

## lambda = C_VM*(V_M/V_P)/h2, so a given (lambda, V_M/V_P) implies an h2.
## Above 1 the combination is not realisable: there is a floor
## lambda >= C_VM*(V_M/V_P), which at V_M/V_P = 2.6e-3 is 7.8e-4.
res[, realisable := is.finite(h2_implied) & h2_implied <= 1]
res[, sigma_m2   := VM_VP / U]

## --- the headline number: biobanks needed for 80% power --------------------
need <- rbindlist(map(LAMBDA_VALS[LAMBDA_VALS > 0], function(lam) {
  q <- n_biobanks_for_power(prep, vm_vp = VM_VP, lambda = lam,
                            power = 0.8, z_a = Z_ALPHA, N_PER = N_PER)
  data.table(lambda     = lam,
             h2_implied = C_VM * VM_VP / lam,
             h2_sing    = q$h2_sing,
             t0         = q$t0,
             se_need    = q$se_need,
             N_BB_80    = q$N_BB,
             N_total_80 = q$N_total)
}))
need[, realisable := is.finite(h2_implied) & h2_implied <= 1]

cat(glue("Realisability floor: lambda >= C_VM * V_M/V_P = ",
         "{signif(C_VM * VM_VP, 3)} for h2 <= 1"), "\n\n")

cat("--- biobanks of", N_PER, "needed for 80% power (Fisher FLOOR) ---\n")
print(need[order(lambda)])

cat("\n--- power grid, realisable combinations only ---\n")
print(res[realisable == TRUE][order(lambda, N_BB)])

## lambda = 0 is a type-I sanity check: power should come back at
## pnorm(-Z_ALPHA) = 0.025. If it doesn't, something upstream is wrong.
cat("\n--- lambda = 0 (should be ~0.025) ---\n")
print(res[lambda == 0][order(N_BB)][, .(N_BB, h2_sing, se, power)])
