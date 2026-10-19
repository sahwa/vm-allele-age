library(data.table)
library(purrr)
library(ggplot2)
library(patchwork)
f = glue::glue

COMMON="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/common"
DATA="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v4.0"
FIGS="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v4.0"
OUT = file.path(DATA, "collated")
dir.create(OUT, showWarnings = FALSE)
source("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/diagGREML.R")
source(file.path(COMMON, "fit_sel_reml_joint.R"))

S_vals = c(0.002, 0.005, 0.01)
N_REPS = 25
K = 50
n_sample = 2500
H2 = 0.5

get_counts = function(S, REP) {
  set.seed(1e4 * which(S_vals == S) + REP)           # reproducible phenotype noise
  sing <- fread(file.path(DATA, f("sing_K{S}_seed{REP}.tsv")))   # carrier, age, beta
  stopifnot(max(sing$carrier) <= n_sample)

  # ---- Age bins: equal-count, and each age falls in exactly one bin ---------
  breaks <- unique(quantile(sing$age, probs = seq(0, 1, length.out = K + 1), type = 1))
  sing[, bin := cut(age, breaks, include.lowest = TRUE, labels = FALSE)]
  K_eff <- length(breaks) - 1                        # can be < K if many ages are tied

  # ---- Per-bin summaries, including the oracle (true mean beta^2) ----------
  bins <- sing[, .(mid = mean(age), oracle_b2 = mean(beta^2), n_sing = .N), keyby = bin]

  m_orc_hyp = lm(data=bins, (1 / oracle_b2) ~ mid)
  oracle_sigma2_m_hyp = 1 / coef(m_orc_hyp)[[1]]
  oracle_lambda = coef(m_orc_hyp)[[2]] / coef(m_orc_hyp)[[1]]
  data.table(K_SEL = S, lambda = oracle_lambda)
}

params <- CJ(S = S_vals, REP = setdiff(seq_len(N_REPS), 1))

res <- pmap_dfr(params, get_counts)
