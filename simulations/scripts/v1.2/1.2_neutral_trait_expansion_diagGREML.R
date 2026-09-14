source("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/diagGREML.R")
DATA="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.2/replicates/rep0"
FIGS = "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v1.2"


C <- fread(file.path(DATA, "1.2_singleton_counts.csv"))
pheno = fread(file.path(DATA, "1.2_phenotypes.csv"))

bins <- names(C)

A <- c(
  lapply(bins[1:3], function(b) C[[b]]),
  list(old = rowSums(C[, bins[4:8], with = FALSE]))
)

A <- lapply(A, function(x) x / mean(x))

names(A) <- c(bins[1:3], "old")
X <- matrix(1, nrow = nrow(C), ncol = 1)

fit <- fit_diagGREML(y = pheno$y, A = A, X = X,
                     constraint = FALSE, magic0316 = TRUE)
data.table(comp = fit$Vlistnames, est = fit$varcmp,
           se = sqrt(diag(fit$Hi)))

library(ggplot2)
library(data.table)
library(patchwork)

# ---- panel A: variance components with 95% CIs ----
res <- data.table(
    comp = factor(fit$Vlistnames[1:4], levels = fit$Vlistnames[1:4]),
    h2   = fit$h2,
    se   = fit$se_h2
)[, `:=`(lo = h2 - 1.96 * se, hi = h2 + 1.96 * se)]

pA <- ggplot(res, aes(x = comp, y = h2)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.5) +
    geom_point(size = 2.6, colour = "#B2182B") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(
        x = "Allele age bin (generations)",
        y = expression("Proportion of"~italic(V[P])),
        title = "A. Singleton variance components",
        subtitle = sprintf("LRT = %.2f on %d df, p = %.2f (n = %s)",
                           fit$LRT, fit$df, fit$p_LRT,
                           format(fit$n, big.mark = ","))
    ) +
    theme_classic(base_size = 11) +
    theme(
        panel.border = element_rect(colour = "grey75", fill = NA, linewidth = 0.5),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, colour = "grey30")
    )

# ---- panel B: sampling correlation between components ----
sd_v <- sqrt(diag(fit$Hi))
Cor  <- fit$Hi / outer(sd_v, sd_v)
dimnames(Cor) <- list(fit$Vlistnames, fit$Vlistnames)

cor_dt <- as.data.table(as.table(Cor))
setnames(cor_dt, c("row", "col", "r"))
cor_dt[, `:=`(row = factor(row, levels = fit$Vlistnames),
              col = factor(col, levels = rev(fit$Vlistnames)))]

pB <- ggplot(cor_dt, aes(row, col, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", r)), size = 2.9) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), name = "r") +
    labs(x = NULL, y = NULL, title = "B. Sampling correlation of estimates") +
    theme_classic(base_size = 11) +
    theme(
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 11)
    )

p <- pA + pB + plot_layout(widths = c(1, 1.15))

ggsave(file.path(FIGS, "diagGREML_singleton_components.png"),
       p, width = 10, height = 4.2, dpi = 200)