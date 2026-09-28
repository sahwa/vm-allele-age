library(data.table)
library(stringr)
library(purrr)
library(scales)
library(future)
library(furrr)

f <- glue::glue

BASE_DATA <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.2"
source("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/diagGREML.R")

bin_midpoint <- function(x) {
  out <- rep(NA_real_, length(x))
  ok  <- str_detect(x, "^[0-9.]+-[0-9.]+$")          # "error" and "205+" stay NA
  out[ok] <- map_dbl(str_split(x[ok], "-"), ~ mean(as.numeric(.x)))
  out
}

runs  <- fread("1.2_neutral_trait_expansion_msprime_n_sweep.params",
               col.names = c("n_dip", "L"), colClasses = "character")
Svals <- c(0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2)
NBs   <- c(1, 2, 3, 4, 6)

# keep existing (n, L) pairs; cross only with S and nb
grid <- runs[, CJ(S = Svals, nb = NBs), by = .(n_dip, L)]

fit_one <- function(n_dip, L, S, nb) {
  run_dir <- file.path(BASE_DATA,
                       f("n{str_remove_all(n_dip, '_')}_L{scientific(as.numeric(L))}"))
  c_file  <- file.path(run_dir, f("singleton_counts_nb{nb}.csv"))
  y_file  <- file.path(run_dir, f("phenotypes_S{scientific(S)}.csv"))
  stopifnot(file.exists(c_file), file.exists(y_file))

  C     <- fread(c_file, header = TRUE)
  pheno <- fread(y_file, header = TRUE)
  stopifnot(nrow(C) == nrow(pheno))

  keep <- names(C)[colSums(C) > 0]                    # drop empty bins
  C    <- C[, ..keep]

  mean_c <- sapply(C, mean)
  M_t    <- sapply(C, sum)
  A      <- Map(`/`, as.list(C), mean_c)
  X      <- matrix(1, nrow = nrow(C), ncol = 1)

  fit <- fit_diagGREML(y = pheno$y, A = A, X = X,
                       constraint = FALSE, magic0316 = TRUE)

  d <- data.table(comp = fit$Vlistnames,
                  est  = as.numeric(fit$varcmp),
                  se   = sqrt(diag(fit$Hi)))
  stopifnot(all(setdiff(d$comp, "error") %in% names(C)))

  d <- merge(d, data.table(comp = names(C), M_t = M_t, mean_c = mean_c),
             by = "comp", all.x = TRUE, sort = FALSE)
  d[, `:=`(t_mid    = bin_midpoint(comp),
           sigma2_b = est / mean_c,                    # == est * N / M_t
           n_dip = n_dip, L = L, S = S, nb = nb)]

  reg <- d[!is.na(t_mid) & est > 0]
  out <- list(fit = d, Hi = fit$Hi, n_used = nrow(reg),
              lm = NULL, sigma2_m = NA_real_, s_hat = NA_real_)
  if (nrow(reg) < 3) return(out)

  m <- lm(log(sigma2_b) ~ t_mid, data = reg, weights = (est / se)^2)
  modifyList(out, list(lm = m,
                       sigma2_m = unname(exp(coef(m)[1])),
                       s_hat    = unname(-coef(m)[2])))
}

n_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
plan(multicore, workers = n_workers)   # forks; Linux only, not inside RStudio

res <- future_pmap(
  grid,
  safely(fit_one),
  .options = furrr_options(
    seed       = TRUE,
    scheduling = Inf,                  # one task per row = dynamic load balancing
    packages   = c("data.table", "stringr", "purrr", "scales", "glue")
  )
)


library(data.table)
library(ggplot2)
library(stringr)
library(purrr)

# ---- inputs ---------------------------------------------------------------
# True per-mutation effect variance, matching what goes into D:
#   causal singletons only  -> SIGMA_BETA^2
#   all singletons          -> PI_CAUSAL * SIGMA_BETA^2
TRUE_SIGMA2_M <- 1e-4          # <-- set from your simulation parameters
OUT_DIR <- "figures"
dir.create(OUT_DIR, showWarnings = FALSE)

# ---- tidy ------------------------------------------------------------------
comp_dt <- copy(summary_dt)
comp_dt[, n := as.numeric(str_remove_all(n_dip, "_"))]
comp_dt[, L := as.numeric(L)]
comp_dt[, n_lab := factor(format(n, big.mark = ",", scientific = FALSE),
                          levels = format(sort(unique(n)), big.mark = ",", scientific = FALSE))]
comp_dt[, S_lab := factor(paste0("s = ", S), levels = paste0("s = ", sort(unique(S))))]
comp_dt[, is_old := str_detect(comp, "\\+$")]
comp_dt[is_old == TRUE, t_mid := 250]                  # plotting position only

slopes <- rbindlist(imap(compact(res), function(r, id) {
  if (is.null(r$lm)) return(NULL)
  cf <- summary(r$lm)$coefficients
  r$fit[1, .(n_dip, L, S, nb)][, `:=`(
    s_hat     = r$s_hat,
    s_se      = cf["t_mid", "Std. Error"],
    sigma2_m  = r$sigma2_m,
    int_se    = cf["(Intercept)", "Std. Error"],
    n_used    = r$n_used
  )]
}))
slopes[, n := as.numeric(str_remove_all(n_dip, "_"))]
slopes[, L := as.numeric(L)]
slopes[, n_lab := factor(format(n, big.mark = ",", scientific = FALSE),
                         levels = format(sort(unique(n)), big.mark = ",", scientific = FALSE))]

L_PLOT  <- 1e8   # the complete series
NB_PLOT <- 4

# ---- 1. raw variance components per bin -----------------------------------
p1 <- ggplot(comp_dt[comp != "error" & L == L_PLOT & nb == NB_PLOT],
             aes(t_mid, est, shape = is_old)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_pointrange(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se), size = 0.3) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
                     labels = c("dated bin", "205+ catch-all"), name = NULL) +
  facet_grid(n_lab ~ S_lab, scales = "free_y") +
  labs(x = "Bin midpoint (generations)",
       y = expression("REML component " * sigma[b(t)]^2),
       title = sprintf("Variance components by age bin (nb = %d, L = %g)", NB_PLOT, L_PLOT)) +
  theme_bw() + theme(legend.position = "bottom")

# ---- 2. the grant regression: log per-mutation variance vs t --------------
reg_dt <- comp_dt[comp != "error" & is_old == FALSE & L == L_PLOT & nb == NB_PLOT]
reg_dt[, log_s2 := ifelse(est > 0, log(sigma2_b), NA_real_)]
reg_dt[, log_se := se / est]                            # delta method

truth <- unique(reg_dt[, .(S, S_lab)])[, `:=`(a = log(TRUE_SIGMA2_M), b = -S)]
fits  <- merge(slopes[L == L_PLOT & nb == NB_PLOT],
               unique(comp_dt[, .(S, S_lab)]), by = "S")
fits[, `:=`(a = log(sigma2_m), b = -s_hat)]

p2 <- ggplot(reg_dt[!is.na(log_s2)], aes(t_mid, log_s2)) +
  geom_pointrange(aes(ymin = log_s2 - 1.96 * log_se, ymax = log_s2 + 1.96 * log_se),
                  size = 0.3) +
  geom_abline(data = truth, aes(intercept = a, slope = b),
              colour = "firebrick", linetype = "dashed") +
  geom_abline(data = fits, aes(intercept = a, slope = b), colour = "steelblue") +
  facet_grid(n_lab ~ S_lab) +
  labs(x = "Bin midpoint (generations)",
       y = expression(log(sigma[b(t)]^2 / bar(x)[t])),
       title = "Per-mutation variance by age: fit (blue) vs truth (red dashed)",
       caption = "Bins with negative estimates are omitted from the plot and the fit") +
  theme_bw()

# ---- 3. recovery of s ------------------------------------------------------
p3 <- ggplot(slopes[L == L_PLOT], aes(S, s_hat, colour = n_lab)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dashed") +
  geom_pointrange(aes(ymin = s_hat - 1.96 * s_se, ymax = s_hat + 1.96 * s_se),
                  position = position_dodge(width = 5e-4), size = 0.3) +
  facet_wrap(~ paste0("nb = ", nb)) +
  labs(x = "True s", y = expression(hat(s)), colour = "n",
       title = "Recovery of the selection parameter") +
  theme_bw()

# ---- 4. recovery of sigma2_m ----------------------------------------------
p4 <- ggplot(slopes[L == L_PLOT], aes(S, sigma2_m / TRUE_SIGMA2_M, colour = n_lab)) +
  geom_hline(yintercept = 1, colour = "grey50", linetype = "dashed") +
  geom_pointrange(aes(ymin = exp(log(sigma2_m) - 1.96 * int_se) / TRUE_SIGMA2_M,
                      ymax = exp(log(sigma2_m) + 1.96 * int_se) / TRUE_SIGMA2_M),
                  position = position_dodge(width = 5e-4), size = 0.3) +
  scale_y_log10() +
  facet_wrap(~ paste0("nb = ", nb)) +
  labs(x = "True s", y = expression(hat(sigma)[m]^2 / sigma[m]^2), colour = "n",
       title = "Recovery of per-mutation effect variance (1 = unbiased)") +
  theme_bw()

ggsave(file.path(OUT_DIR, "components_by_bin.pdf"),  p1, width = 11, height = 9)
ggsave(file.path(OUT_DIR, "log_regression.pdf"),     p2, width = 11, height = 9)
ggsave(file.path(OUT_DIR, "s_recovery.pdf"),         p3, width = 9,  height = 5)
ggsave(file.path(OUT_DIR, "sigma2m_recovery.pdf"),   p4, width = 9,  height = 5)
