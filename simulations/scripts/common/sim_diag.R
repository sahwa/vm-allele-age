sim_diag <- function(S, REP, H2, CV, N_inds, K, t_max, LAMBDA, type) {
  set.seed(REP + 1e3 * which(Svals == S) + round(CV * 1e4) + K * 1e5 + t_max * 1e6)
  t_mid <- (seq_len(K) - 0.5) * t_max / K
  shape <- 1 / CV^2
  C <- sapply(seq_len(K), function(k) rpois(N_inds, rgamma(N_inds, shape, rate = shape / LAMBDA)))
  colnames(C) <- paste0("t", round(t_mid))
  if (type == "exponential") {
    s2 <- SIGMA_BETA^2 * exp(-S * t_mid)
  } else if (type == "hyperbolic") {
    s2 = SIGMA_BETA^2 / (1 + S * t_mid)
  } else {
    stop("Type should be one of hyperbolic or exponential!")
  }
  g  <- rowSums(sqrt(sweep(C, 2, s2, `*`)) * matrix(rnorm(N_inds * K), N_inds, K))
  Y  <- g + rnorm(N_inds, 0, sqrt(var(g) * (1 - H2) / H2))
  list(C = C, Y = Y, t_mid = t_mid, s2_true = s2)
}
