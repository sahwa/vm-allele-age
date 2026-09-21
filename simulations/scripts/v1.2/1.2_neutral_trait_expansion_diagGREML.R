f = glue::glue
library(purrr)

VERSION = "1.2.1"

BASE = "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations"
PROGRAMS = file.path(BASE, "programs")
SCRIPTS = file.path(BASE, "scripts")
SIM_DATA <- file.path(BASE, "data/v1.2")

source(file.path(PROGRAMS, "diagGREML.R"))

params = fread(file.path(SCRIPTS, "v1.2/1.2_neutral_trait_expansion_msprime_n_sweep.params"))
N_final = unique(str_remove_all(params$V1, "_"))
L = unique(params$V2)
Nb = setdiff(1:6, 5)

args = expand.grid(Nsim = N_final, L = L)

summary_DT = rbindlist(purrr::pmap(args, function(Nsim, L) {
    singleton_file = file.path(SIM_DATA, f("n{Nsim}_L{scientific(L)}/bin_summary.csv"))
    if (file.exists(singleton_file)) {
        return(fread(singleton_file))
    }
    }
))

summary_DT[, n_eff_n_dip_ratio := n_eff / n_dip]
summary_DT
