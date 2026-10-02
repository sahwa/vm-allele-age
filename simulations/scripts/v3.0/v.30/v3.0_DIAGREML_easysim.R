library(future)
library(furrr)
library(purrr)
library(patchwork)

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
sigma2_e = 0.2

H2_vals = c(0.1, 0.5, 1)
Svals <- c(0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2)
CV_Vals = seq(0.15, 0.5, 0.15)

## each individual should have Poisson counts for singletons
sim_diag = function(S, REP, seed, H2, CV) {
    set.seed(REP + 1000 * which(Svals == S) + seed)
    ## Then, for each bin, each individuals score is just the sum of the number of counts of draws from that bin
    C = sapply(seq_len(K), function(k){
    rpois(N_inds, rgamma(N_inds, shape = SHAPE, rate = SHAPE / LAMBDA))    # independent across bins
    })
    colnames(C) <- paste0("t", round(t_mid))
    M_t = colSums(C)
    s2 = SIGMA_BETA * exp(-S*t_mid)
    G = sqrt(sweep(C, 2, s2, `*`)) * matrix(rnorm(N_inds * K), N_inds, K)   # exact: sum of x N(0, s2)
    g = rowSums(G)
    V_E = var(g) * (1 - H2) / H2
    Y = g + rnorm(N_inds, 0, sd = sqrt(V_E))
    colnames(C) <- paste0("t", round(t_mid))

    list(C=C, Y=Y, t_mid = t_mid, s2_true = s2, M_t = M_t, REP = REP, H2 = H2)
}

run_diagGREML = function(S, REP, seed, H2, CV) {
    d = sim_diag(S, REP, seed, H2, CV)
    Y = d$Y; C = d$C; T = d$t_mid; s2_true = d$s2_true; M_t = d$M_t
    m = colMeans(C)
    A = setNames(lapply(seq_len(ncol(C)), function(k) {
        C[, k] = C[,k] / m[k]
    }), colnames(C))
    X = matrix(1, nrow(C), 1)
    fit = suppressMessages(fit_diagGREML(y=Y, X=X, A=A, constraint = FALSE, magic0316 = TRUE))
    idx = match(names(A), fit$Vlistnames)
    res <- data.table(comp = names(A), mid = T, m = m, s2_true = s2_true,
                  est = as.numeric(fit$varcmp)[idx],
                  se  = sqrt(diag(fit$Hi))[idx])
    n_neg <- sum(res$est <= 0)

    res <- res[est > 0]
    res[, rcomp := log(est / m)]
    model = lm(data = res, rcomp~mid)
    cf = coef(model)
    res[, rcomp := log(est / m)]
    modelnn = lm(data = res, rcomp~mid)
    cf = coef(model)

    list(sigma2_m = exp(unname(cf[1])), s_hat = -unname(cf[2]),
     s2_true = s2_true, n_pos_comp = nrow(res), H2 = H2)
}

N_REPS = 100

params = CJ(S = Svals, REP=1:N_REPS, H2 = H2_vals, CV = CV_Vals)
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
res_est <- rbindlist(map(res, function(x) {
  r <- x$result
  if (is.null(r)) return(data.table(sigma2_m = NA_real_, s_hat = NA_real_, n_pos_comp = NA_integer_, H2 = NA_integer_))
  data.table(sigma2_m = r$sigma2_m, s_hat = r$s_hat, n_pos_comp = r$n_pos_comp)
}))

res_total <- cbind(params, res_est)
res_total[, S_f := factor(S, levels = Svals)]

p1 <- ggplot(res_total, aes(S_f, s_hat)) +
  geom_boxplot(outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom() +
  geom_point(data = unique(res_total[, .(S_f, S)]), aes(y = S),
             colour = "red", size = 3) +
  labs(x = "True s", y = expression(hat(s))) +
  facet_wrap(~H2)

p2 <- ggplot(res_total, aes(S_f, sigma2_m)) +
  geom_boxplot(outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom() +
  geom_hline(yintercept = SIGMA_BETA^2, colour = "red",
             linetype = "dashed", linewidth = 1) +
  labs(x = "True s", y = expression(hat(sigma)[m]^2))

p1 + p2


ggsave("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v3.0/S_Est_model1.png", p1)
ggsave("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v3.0/V_M_Est_model1.png", p2)
