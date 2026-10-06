FIGS="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v5.0"

p_bias_nneg = dat[, .(mbias = mean(bias), mneg = mean(mean_n_neg)), .(H2)][order(mbias)] %>%
    ggplot(aes(x=mneg, y=mbias, label=H2)) +
    geom_point() +
    geom_label()

ggsave(file.path(FIGS, "p_bias_nneg.png"), p_bias_nneg)


reps = fread("reps_N300000.csv")
  summ <- reps[, {
    sh <- s_hat[is.finite(s_hat)]; z <- sh - S
    sm <- sigma2_m[is.finite(sigma2_m)]
    .(n_total   = .N, n_fit = length(sh), fail_rate = 1 - length(sh) / .N,
      mean_n_neg = mean(n_neg, na.rm = TRUE),
      gm_sigma2_m  = if (length(sm)) exp(mean(log(sm))) else NA_real_,
      gmedian_sigma2_m  = if (length(sm)) exp(median(log(sm))) else NA_real_,
      fold_lo      = if (length(sm)) quantile(sm, 0.025) else NA_real_,
      fold_hi      = if (length(sm)) quantile(sm, 0.975) else NA_real_)
  }, by = .(S, H2, CV, N_inds)]
  setorder(summ, CV, -H2, S)

su
