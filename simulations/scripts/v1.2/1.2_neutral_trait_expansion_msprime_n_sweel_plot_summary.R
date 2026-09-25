library(data.table)
library(ggplot2)
library(patchwork)
library(scales)

DATA <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.2"
FIGS <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v1.2"
dir.create(FIGS, recursive = TRUE, showWarnings = FALSE)

N_VALS <- c(20000, 50000, 100000, 200000, 500000)
L_VALS <- c("1e+08", "1e+09")

# ---------------------------------------------------------------
# Load the sweep outputs
# ---------------------------------------------------------------
read_if <- function(f) if (file.exists(f)) fread(f) else NULL

bins <- rbindlist(lapply(N_VALS, function(n) rbindlist(lapply(L_VALS, function(L)
    read_if(file.path(DATA, sprintf("n%d_L%s", n, L), "bin_summary.csv")))
)), fill = TRUE)

supply <- rbindlist(lapply(N_VALS, function(n) rbindlist(lapply(L_VALS, function(L)
    read_if(file.path(DATA, sprintf("n%d_L%s", n, L), "old_variant_supply.csv")))
)), fill = TRUE)

bins[,   L_lab := factor(sprintf("L = %.0e", L))]
supply[, L_lab := factor(sprintf("L = %.0e", L))]

# ---------------------------------------------------------------
# A. Per-person singleton count is flat in n
# ---------------------------------------------------------------
pA <- ggplot(bins[nb == 1 & bin == "0-205"],
             aes(n_dip, mean_c, colour = L_lab)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    scale_x_log10(labels = comma) +
    scale_y_log10() +
    scale_colour_brewer(palette = "Set1", name = NULL) +
    labs(x = "Sample size (diploids)",
         y = "Singletons per person",
         title = "A. Per-person singleton burden barely changes",
         subtitle = "25x more people, ~25% fewer singletons each") +
    theme_classic(base_size = 11) +
    theme(panel.border = element_rect(colour = "grey75", fill = NA, linewidth = 0.5),
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, colour = "grey30"),
          legend.position = "bottom")

# ---------------------------------------------------------------
# B. Old-variant supply collapses
# ---------------------------------------------------------------
pB <- ggplot(supply[threshold %in% c(50, 200, 1000)],
             aes(n_dip, mean_per_person + 1e-4,
                 colour = factor(threshold), linetype = L_lab)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.2) +
    scale_x_log10(labels = comma) +
    scale_y_log10(labels = label_number(accuracy = 0.001)) +
    scale_colour_viridis_d(option = "plasma", end = 0.85,
                           name = "age >") +
    scale_linetype(name = NULL) +
    labs(x = "Sample size (diploids)",
         y = "Old singletons per person",
         title = "B. Old singletons are stripped out",
         subtitle = expression("age > 200 generations declines as "~n^-1.16)) +
    theme_classic(base_size = 11) +
    theme(panel.border = element_rect(colour = "grey75", fill = NA, linewidth = 0.5),
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, colour = "grey30"),
          legend.position = "bottom", legend.box = "vertical")

# ---------------------------------------------------------------
# C. The age window narrows
# ---------------------------------------------------------------
edges <- bins[nb == 6 & is.finite(bin_hi)]
pC <- ggplot(edges, aes(factor(n_dip), ymin = bin_lo, ymax = bin_hi,
                        fill = factor(bin_lo))) +
    geom_linerange(aes(x = factor(n_dip)), linewidth = 6,
                   colour = "grey35", alpha = 0.8,
                   position = position_dodge(width = 0)) +
    facet_wrap(~ L_lab) +
    scale_y_log10() +
    guides(fill = "none") +
    labs(x = "Sample size (diploids)",
         y = "Allele age (generations)",
         title = "C. Six equal-count bins span a narrowing window",
         subtitle = "The axis the selection signal lives on") +
    theme_classic(base_size = 11) +
    theme(panel.border = element_rect(colour = "grey75", fill = NA, linewidth = 0.5),
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, colour = "grey30"))

# ---------------------------------------------------------------
# D. Precision vs lever arm, normalised
# ---------------------------------------------------------------
trade <- merge(
    bins[nb == 1 & bin == "0-205", .(n_dip, L, n_eff)],
    supply[threshold == 200, .(n_dip, L, old = mean_per_person)],
    by = c("n_dip", "L"))
trade[, `:=`(precision = sqrt(n_eff / n_eff[n_dip == min(n_dip)]),
             lever     = old / old[n_dip == min(n_dip)]), by = L]

pD <- ggplot(melt(trade, id.vars = c("n_dip", "L"),
                  measure.vars = c("precision", "lever")),
             aes(n_dip, value, colour = variable,
                 linetype = factor(sprintf("L = %.0e", L)))) +
    geom_hline(yintercept = 1, colour = "grey60", linetype = "dotted") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.2) +
    scale_x_log10(labels = comma) +
    scale_y_log10() +
    scale_colour_manual(values = c(precision = "#2166AC", lever = "#B2182B"),
                        labels = c(precision = "precision (sqrt n_eff)",
                                   lever = "old-variant supply"),
                        name = NULL) +
    scale_linetype(name = NULL) +
    labs(x = "Sample size (diploids)",
         y = "Relative to n = 20,000",
         title = "D. The trade-off",
         subtitle = "Modest gain in precision, large loss of lever arm") +
    theme_classic(base_size = 11) +
    theme(panel.border = element_rect(colour = "grey75", fill = NA, linewidth = 0.5),
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, colour = "grey30"),
          legend.position = "bottom", legend.box = "vertical")

p <- (pA | pB) / (pC | pD)
ggsave(file.path(FIGS, "sample_size_tradeoff.png"), p,
       width = 12, height = 9, dpi = 200)

# ---------------------------------------------------------------
# Numbers for the text
# ---------------------------------------------------------------
print(dcast(bins[nb == 1 & bin == "0-205"], n_dip ~ L_lab, value.var = "mean_c"))
print(dcast(supply[threshold == 200], n_dip ~ L_lab, value.var = "mean_per_person"))
print(trade[, .(n_dip, L, n_eff, old, precision, lever)])
