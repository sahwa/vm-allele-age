COMMON_DIR="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/common"
source(file.path(COMMON_DIR, "fit_sel_reml_joint.R"))
source(file.path(COMMON_DIR, "sim_diag.R"))

N_REPS = 1:1e4
Svals   <- c(0, 1e-4, 1e-3, 1e-2, 3e-2)
H2 = 1
N_inds = 5e5
K_vals = c(3, 8, 20)
t_max = 400
CV = 1
LAMBDA     <- 1000
SIGMA_BETA = 1
params = CJ(REP = N_REPS, S = Svals)

map(1:N_REPS, function(x) {
    d = sim_diag(S = S, REP = REP, t_max = t_max, K = K, CV = CV, N_inds = N_inds, H2 = H2)
    fit_log = fit_sel_reml(y=d$Y, C=d$C, t=d$t_mid, curve="exp")
    fit_hyperbolic = fit_sel_reml(y=d$Y, C=d$C, t=d$t_mid, curve="hyp")
}
