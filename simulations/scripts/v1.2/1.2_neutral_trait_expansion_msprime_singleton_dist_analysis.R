## ---------------------------------------------------------------------------
## Talk figure: how much does a person's singleton AGE PROFILE really vary?
##
## Self-contained -- reads the counts, builds the two simulation models, makes
## one two-row figure. Nothing else from the diagnostic scripts is needed.
##
## The point in one sentence: people differ a lot in HOW MANY singletons they
## carry, and very little in WHICH AGES those are -- and the old simulation had
## that the wrong way round.
## ---------------------------------------------------------------------------

library(data.table); library(ggplot2); library(patchwork)

DATA_FILE <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.2/n500000_L1e+08/singleton_counts.csv"
FIGS      <- "."
CV_SIM    <- 0.65      # the CV both simulation models were run with
N_LINES   <- 100       # individuals to draw as spaghetti
SEED      <- 1

dat <- fread(DATA_FILE)

## bin midpoints from "lo-hi" / "lo+" labels; geometric centre for log bins
t_mid <- vapply(names(dat), function(s) {
  if (grepl("\\+$", s)) return(1.5 * as.numeric(sub("\\+$", "", s)))
  e <- as.numeric(strsplit(s, "-")[[1]])
  if (e[1] <= 0) mean(e) else sqrt(e[1] * e[2])
}, numeric(1))

C <- as.matrix(dat); C <- C[rowSums(C) > 0, , drop = FALSE]
N <- nrow(C); K <- ncol(C)
lam_k <- colMeans(C)                       # mean count per bin, matched by both sims

s

## --- the two simulation models ----------------------------------------------
## Identical per-bin means and identical CV. The ONLY difference is whether the
## gamma factor is drawn once per person or once per person PER BIN.
set.seed(SEED)
sh <- 1 / CV_SIM^2
C_shr <- { f <- rgamma(N, sh, rate = sh)                      # one factor, all bins
           sapply(seq_len(K), function(k) rpois(N, f * lam_k[k])) }
C_ind <- sapply(seq_len(K), function(k)                       # a fresh factor per bin
           rpois(N, rgamma(N, sh, rate = sh / lam_k[k])))

SRC <- c("msprime (real)", "sim: shared factor", "sim: independent (old)")
mats <- list(C, C_shr, C_ind); names(mats) <- SRC

## --- per-person summaries ---------------------------------------------------
## total    = how many singletons you carry        -> the signal for V_M
## mean age = the typical age of those singletons  -> the signal for s
summ_person <- rbindlist(lapply(SRC, function(s) {
  M <- mats[[s]]; Tt <- rowSums(M); M <- M[Tt > 0, , drop = FALSE]; Tt <- Tt[Tt > 0]
  data.table(source = s, total = Tt,
             mean_age = exp(as.vector(M %*% log(t_mid)) / Tt))
}))
summ_person[, source := factor(source, levels = SRC)]

cvtab <- summ_person[, .(`how many (total count)` = sd(total) / mean(total),
                         `which ages (mean age)`  = sd(mean_age) / mean(mean_age)),
                     by = source]
print(cvtab)

## --- panel A: 100 individual profiles, as proportions ------------------------
set.seed(SEED)
spag <- rbindlist(lapply(SRC, function(s) {
  M <- mats[[s]]; Tt <- rowSums(M); M <- M[Tt > 0, , drop = FALSE]; Tt <- Tt[Tt > 0]
  i <- sample(nrow(M), N_LINES)
  data.table(source = s, id = rep(i, each = K), t = rep(t_mid, length(i)),
             p = as.vector(t(M[i, , drop = FALSE] / Tt[i])))
}))
spag[, source := factor(source, levels = SRC)]
mean_prof <- data.table(t = t_mid, p = lam_k / sum(lam_k))

pA <- ggplot(spag[p > 0], aes(t, p, group = id)) +
  geom_line(alpha = 0.10, linewidth = 0.3) +
  geom_line(data = mean_prof, aes(t, p), inherit.aes = FALSE,
            colour = "#D7191C", linewidth = 1) +
  facet_wrap(~source) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "allele age (generations)", y = "fraction of a person's singletons",
       title = "Each line is one person's singleton age profile",
       subtitle = "Red = population mean. Real people hug it; the old simulation scattered them wildly.") +
  theme_bw(base_size = 13)

## --- panel B: the two kinds of between-person variation ----------------------
pB <- ggplot(summ_person, aes(total, colour = source)) +
  geom_density(linewidth = 0.8) + scale_x_log10() +
  labs(x = "singletons carried", y = NULL, colour = NULL,
       title = "HOW MANY  (this is what gives us V_M)") +
  theme_bw(base_size = 13) + theme(legend.position = "none")

pC <- ggplot(summ_person, aes(mean_age, colour = source)) +
  geom_density(linewidth = 0.8) + scale_x_log10() +
  labs(x = "mean age of those singletons (generations)", y = NULL, colour = NULL,
       title = "WHICH AGES  (this is what would give us selection)") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")

lab <- cvtab[, sprintf("%s:  CV %.2f (how many) | CV %.2f (which ages)",
                       source, `how many (total count)`, `which ages (mean age)`)]

fig <- pA / (pB | pC) +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(caption = paste(lab, collapse = "\n"),
                  theme = theme(plot.caption = element_text(hjust = 0, size = 10)))

ggsave(file.path(FIGS, "profile_summary_talk.png"), fig,
       width = 12, height = 8.5, dpi = 160)
fwrite(cvtab, file.path(FIGS, "profile_summary_cv.csv"))
cat("written: profile_summary_talk.png\n")
