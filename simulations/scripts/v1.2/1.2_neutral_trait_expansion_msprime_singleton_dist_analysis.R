## ---------------------------------------------------------------------------
## How much does the singleton age PROFILE really vary between people?
##
## Two things the joint estimator needs, and they are not the same thing:
##   MAGNITUDE of a person's burden   -> identifies sigma2_m
##   SHAPE of their age profile       -> identifies s
##
## This script measures both in the real msprime counts, and puts the two
## simulation models (indep = the old bug, shared = the pessimistic default)
## on the same axes so we can see where reality sits between them.
##
## Input: singleton_counts.csv, one row per person, one column per age bin,
##        column names like "0-1", "1-2", ..., "205+".
## ---------------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(patchwork)

DATA_FILE <- "singleton_counts.csv"
FIGS      <- "."                      # set to your figs dir
S_PROBE   <- 1e-2                     # the s whose signal we are measuring
CV_INDEP  <- 0.65                     # per-bin gamma CV used by the old sim
CV_SHARED <- 0.65                     # per-person gamma CV for the shared model
SEED      <- 1

dat <- fread(DATA_FILE)

## --- bin ages ---------------------------------------------------------------
## Parse "lo-hi" labels into edges, then take the geometric midpoint, which is
## the right centre for log-spaced bins. The open top bin gets 1.5x its lower
## edge (it is nearly empty at this sample size, so the choice barely matters).
parse_mid <- function(nm) {
  vapply(nm, function(s) {
    if (grepl("\\+$", s)) return(1.5 * as.numeric(sub("\\+$", "", s)))
    e <- as.numeric(strsplit(s, "-")[[1]])
    if (e[1] <= 0) mean(e) else sqrt(e[1] * e[2])
  }, numeric(1))
}
bins <- names(dat)
t_mid <- parse_mid(bins)
K <- length(bins)

C <- as.matrix(dat)
T_i <- rowSums(C)
keep <- T_i > 0
C <- C[keep, , drop = FALSE]; T_i <- T_i[keep]
N <- nrow(C)
cat(sprintf("N = %d people, %d bins, mean total = %.1f singletons\n", N, K, mean(T_i)))
cat("bin midpoints (generations): ", paste(round(t_mid, 1), collapse = ", "), "\n")

## --- the identifying statistic ----------------------------------------------
## a_i = w_i(s) / w_i(0) = (sum_k C_ik f_k) / T_i
##     = the average decay weight across the singletons that person carries.
## It is a pure summary of SHAPE: it is unchanged if you double everyone's
## counts. If every person had the same profile, a_i would be a constant and
## s would carry no information at all. So CV(a_i) IS the signal for s, and
## CV(T_i) is the signal for sigma2_m.
fdecay <- function(s, t) 1 / (1 + s * t)        # hyperbolic; exp(-s*t) for the other
a_of <- function(Cm, s) as.vector(Cm %*% fdecay(s, t_mid)) / rowSums(Cm)

a_real <- a_of(C, S_PROBE)

## --- matched simulations ----------------------------------------------------
## Same number of people, same mean profile, same total burden. The ONLY thing
## that differs is how count variation is distributed across bins.
set.seed(SEED)
p_bar <- colMeans(C) / sum(colMeans(C))         # population mean profile
lam_k <- mean(T_i) * p_bar                      # mean count per bin

sim_indep <- function() {
  sh <- 1 / CV_INDEP^2
  sapply(seq_len(K), function(k)
    rpois(N, rgamma(N, sh, rate = sh / max(lam_k[k], 1e-8))))
}
sim_shared <- function() {
  sh  <- 1 / CV_SHARED^2
  f_i <- rgamma(N, sh, rate = sh)               # mean-1 per-person factor
  sapply(seq_len(K), function(k) rpois(N, f_i * max(lam_k[k], 1e-8)))
}

C_ind <- sim_indep();  C_shr <- sim_shared()
ok_i  <- rowSums(C_ind) > 0; ok_s <- rowSums(C_shr) > 0
a_ind <- a_of(C_ind[ok_i, , drop = FALSE], S_PROBE)
a_shr <- a_of(C_shr[ok_s, , drop = FALSE], S_PROBE)

cv <- function(x) sd(x) / mean(x)
sig <- data.table(
  source   = c("real (msprime)", "sim: shared", "sim: indep (old)"),
  cv_shape = c(cv(a_real), cv(a_shr), cv(a_ind)),
  cv_total = c(cv(T_i), cv(rowSums(C_shr)[ok_s]), cv(rowSums(C_ind)[ok_i])))
sig[, rel_info_s := (cv_shape / cv_shape[1])^2]   # information for s scales as CV^2
print(sig)

## --- panel A: raw profiles --------------------------------------------------
set.seed(SEED)
idx <- sample(N, 150)
prof <- rbindlist(list(
  data.table(id = rep(idx, each = K), t = rep(t_mid, length(idx)),
             n = as.vector(t(C[idx, , drop = FALSE])), scale = "counts"),
  data.table(id = rep(idx, each = K), t = rep(t_mid, length(idx)),
             n = as.vector(t(C[idx, , drop = FALSE] / T_i[idx])), scale = "proportions")))
prof[, scale := factor(scale, levels = c("counts", "proportions"))]

pA <- ggplot(prof[n > 0], aes(t, n, group = id)) +
  geom_line(alpha = 0.12, linewidth = 0.3) +
  facet_wrap(~scale, scales = "free_y") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "allele age (generations)", y = NULL,
       title = "A. 150 individual age profiles",
       subtitle = "Left: raw counts, spread vertically. Right: same lines as proportions -- they collapse.\nThat gap is the whole argument: people differ in HOW MANY, barely in WHICH AGES.") +
  theme_bw()

## --- panel B: between-bin correlation ---------------------------------------
## The decisive diagnostic. A shared per-person factor makes every pair of bins
## POSITIVELY correlated (both scale with the same lambda_i). Independent
## per-bin draws leave them uncorrelated. Real data should look like the former.
cor_dt <- function(Cm, lab) {
  use <- which(colMeans(Cm) > 0.5)
  r <- cor(Cm[, use, drop = FALSE])
  data.table(source = lab, r = r[upper.tri(r)])
}
cors <- rbindlist(list(cor_dt(C, "real (msprime)"),
                       cor_dt(C_shr, "sim: shared"),
                       cor_dt(C_ind, "sim: indep (old)")))
cors[, source := factor(source, levels = sig$source)]

pB <- ggplot(cors, aes(source, r)) +
  geom_hline(yintercept = 0, colour = "grey60", linetype = "dotted") +
  geom_boxplot(width = 0.5, outlier.size = 0.6) +
  labs(x = NULL, y = "correlation between a pair of age bins",
       title = "B. Are counts in different age bins correlated?",
       subtitle = "Positive = a shared per-person factor scaling the whole profile.\nZero = the independent-per-bin model the old sim used.") +
  theme_bw()

## --- panel C: count vs mean age (the GEL comparison) ------------------------
## GEL report Spearman rho = 0.81 between a person's singleton count and the
## geometric mean age of those singletons. This is our version of that plot.
gm_age <- exp(as.vector(C %*% log(t_mid)) / T_i)
rho <- cor(T_i, gm_age, method = "spearman")
cat(sprintf("Spearman(count, geometric mean age) = %.3f   [GEL report 0.81]\n", rho))

pC <- ggplot(data.table(T_i, gm_age), aes(T_i, gm_age)) +
  geom_hex(bins = 60) + scale_fill_viridis_c(trans = "log10", name = "people") +
  geom_smooth(method = "gam", colour = "red", se = FALSE, linewidth = 0.6) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "singletons carried", y = "geometric mean age of them (gen)",
       title = sprintf("C. Count vs mean age: Spearman rho = %.2f", rho),
       subtitle = "GEL report 0.81 in 47,535 real genomes. Deeper genealogy = more singletons AND older ones,\nso shape and magnitude are entangled -- and both track ancestry.") +
  theme_bw()

## --- panel D: the signal for s ----------------------------------------------
ad <- rbindlist(list(data.table(source = "real (msprime)",   a = a_real),
                     data.table(source = "sim: shared",      a = a_shr),
                     data.table(source = "sim: indep (old)", a = a_ind)))
ad[, source := factor(source, levels = sig$source)]

pD <- ggplot(ad, aes(a, colour = source)) +
  geom_density(linewidth = 0.7) +
  labs(x = sprintf("a_i = average decay weight carried (s = %.3g)", S_PROBE),
       y = "density", colour = NULL,
       title = "D. The information available for s",
       subtitle = paste0("Spread here is the ONLY thing that identifies s. A point mass = s unestimable.\n",
                         sprintf("CV: real %.3f | shared %.3f | indep %.3f  ->  the old sim had %.0fx the information.",
                                 sig$cv_shape[1], sig$cv_shape[2], sig$cv_shape[3],
                                 sig$rel_info_s[3]))) +
  theme_bw() + theme(legend.position = "bottom")

fig <- (pA) / (pB | pC) / (pD) + plot_layout(heights = c(1, 1, 1))
ggsave(file.path(FIGS, "count_shape_diagnostics.png"), fig,
       width = 13, height = 14, dpi = 160)
fwrite(sig, file.path(FIGS, "count_shape_summary.csv"))
cat("written: count_shape_diagnostics.png\n")
