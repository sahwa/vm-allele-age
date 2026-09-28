library(future)
library(furrr)

set.seed(123)
source("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/diagGREML.R")

N_inds = 2e4
MU = 1.28e-8
N_REP = 1e3
PI = 0.1
L = 1e9
H2 = 0.5
SIGMA_BETA = 1
K = 8 # number of bins
t_max = 200

cv = 1
SHAPE = 1/cv^2
LAMBDA = 50
RATE = SHAPE / LAMBDA

Svals <- c(0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2)

t_mid <- (seq_len(K) - 0.5) * t_max / K
rgamma
sigma2_e = 0.2

## each individual should have Poisson counts for singletons
sim_diag = function(S, REP, seed) {
    set.seed(REP + 1000 * which(Svals == S) + seed)
    ## Then, for each bin, each individuals score is just the sum of the number of counts of draws from that bin
    C <- sapply(seq_len(K), function(k){
    rpois(N_inds, rgamma(N_inds, shape = SHAPE, rate = SHAPE / LAMBDA))    # independent across bins
    })
    colnames(C) <- paste0("t", round(t_mid))
    M_t = colSums(C)
    s2 = SIGMA_BETA * exp(-S*t_mid)
    G  <- sqrt(sweep(C, 2, s2, `*`)) * matrix(rnorm(N_inds * K), N_inds, K)   # exact: sum of x N(0, s2)
    Y = rowSums(G) + rnorm(N_inds, 0, sd = sqrt(sigma2_e))
    colnames(C) <- paste0("t", round(t_mid))

    list(C=C, Y=Y, t_mid = t_mid, s2_true = s2, M_t = M_t, REP = REP)
}

run_diagGREML = function(S, REP, seed) {
    d = sim_diag(S, REP, seed)
    Y = d$Y; C = d$C; T = d$t_mid; s2_true = d$s2_true; M_t = d$M_t
    m = colMeans(C)
    A = setNames(lapply(seq_len(ncol(C)), function(k) {
        C[, k] = C[,k] / m[k]
    }), colnames(C))
    X = matrix(1, nrow(C), 1)
    fit = fit_diagGREML(y=Y, X=X, A=A, constraint = FALSE, magic0316 = TRUE)
    idx = match(names(A), fit$Vlistnames)
    res <- data.table(comp = names(A), mid = T, m = m, s2_true = s2_true,
                  est = as.numeric(fit$varcmp)[idx],
                  se  = sqrt(diag(fit$Hi))[idx])
    n_neg <- sum(res$est <= 0)

    resnn <- res[est > 0]
    resnn[, rcomp := log(est / m)]
    modelnn = lm(data = res, rcomp~mid)
    cf = coef(model)

    resp = res
    resnn[, rcomp := log(est / m)]
    modelnn = lm(data = res, rcomp~mid)
    cf = coef(model)





    list(intercept = exp(unname(cf[1])), Slope = -unname(cf[2]), s2_true = s2_true, n_pos_comp = nrow(res))
}

N_REPS = 100

params = CJ(S = Svals, REP=1:N_REPS)
params[, seed := .I]

n_workers = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
plan(multicore, workers = n_workers)

res = future_pmap(
  params,
  safely(run_diagGREML),
  .options = furrr_options(
    seed       = TRUE,
    scheduling = Inf,                  # one task per row = dynamic load balancing
    packages   = c("data.table", "stringr", "purrr", "scales", "glue")
  )
)

res_est = rbindlist(map(res, function(x) {
    data.table(V_M_est = x$result$intercept, S_est = x$result$Slope, n_pos_comp = x$result$n_pos_comp)
}))

res_total = cbind(params, res_est)

p1 = res_total %>%
    ggplot(aes(x=factor(S, levels = Svals), y=S_est)) +
    geom_boxplot(aes(group=factor(S, levels = Svals))) +
    ggbeeswarm::geom_quasirandom() +
    geom_point(aes(x=factor(S, levels = Svals), y=S), colour='red', size=4) +
    xlab("True S") +
    ylab("Estimated S")

p2 = res_total %>%
    ggplot(aes(x=factor(S, levels = Svals), y=V_M_est)) +
    geom_boxplot(aes(group=factor(S, levels = Svals))) +
    ggbeeswarm::geom_quasirandom() +
    geom_point(aes(x=factor(S, levels = Svals), y=1), colour='red', size=4) +
    xlab("True S") +
    ylab("Estimated S")


ggsave("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v3.0/S_Est_model1.png", p1)
ggsave("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v3.0/V_M_Est_model1.png", p2)
