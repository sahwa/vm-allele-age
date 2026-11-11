## ---------------------------------------------------------------------------
# power calc for mutliple biobanks
## ---------------------------------------------------------------------------

library(data.table); library(purrr); library(glue); library(RcppCNPy)

source("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/common/power.R")

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
Z_ALPHA     <- 1.96
nb = 1

N_BBs = c(1, 2, 3, 4, 5, 7, 10, 20, 50, 100)

C = npyLoad(file.path(DATA, glue("C_nb{nb}.npy")))
trep = npyLoad(file.path(DATA, glue("trep_nb{nb}.npy")))

paras

rbindlist(map(LAMBDA_VALS, function(lam) {
    prep <- prepare_counts(
        npyLoad(file.path(DATA, glue("C_nb{nb}.npy"))),
        npyLoad(file.path(DATA, glue("trep_nb{nb}.npy")))
    )

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
        sigma_m_2 = VM_VP / U
    )
})

