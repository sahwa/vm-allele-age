library(data.table)
library(ggplot2)
library(patchwork)
library(scales)

FIGS = "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v1.1"

mph <- data.table(
    bin = c("0-100","100-1000","1000-10000","10000-50000","50000-100000",
            "100000-200000","200000-500000","500000-inf","err"),
    var = c(956373, 389544, 85455.1, 27181, 8575.61,
            713.282, 713.282, 713.282, 5706.25),
    pve = c(0.6484, 0.264102, 0.0579367, 0.0184281, 0.00581408,
            0.00048359, 0.00048359, 0.000483589, 0.00386871)
)

BINS <- c(0, 1e2, 1e3, 1e4, 5e4, 1e5, 2e5, 5e5, Inf)
TT   <- 2e5   # 2*Ne for the v1.1 constant-size sims

K <- TT * (exp(-head(BINS, -1) / TT) -
           ifelse(is.infinite(BINS[-1]), 0, exp(-BINS[-1] / TT)))

truth <- data.table(
    bin = mph$bin[1:8],
    true_share = K / sum(K)
)

d <- merge(mph[bin != "err"], truth, by = "bin")
d[, est_share := pve / mph[bin != "err", sum(pve)]]
d[, bin := factor(bin, levels = truth$bin)]
d[, pinned := abs(var - 713.282) < 1e-3]

long <- melt(d, id.vars = c("bin", "pinned"),
             measure.vars = c("true_share", "est_share"),
             variable.name = "source", value.name = "share")
long[, source := factor(source, levels = c("true_share", "est_share"),
                        labels = c("True", "MPH estimate"))]

# ---- A: estimated vs true share by age bin ----
pA <- ggplot(long, aes(bin, share, colour = source, group = source)) +
    geom_line(linewidth = 0.6) +
    geom_point(aes(shape = pinned), size = 2.6) +
    scale_y_log10(labels = percent_format(accuracy = 0.01)) +
    scale_colour_manual(values = c("True" = "grey30",
                                   "MPH estimate" = "#B2182B"), name = NULL) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
                       labels = c("estimated", "pinned"), name = NULL) +
    labs(x = "Allele age bin (generations)",
         y = "Share of genetic variance",
         title = "A. Age profile is inverted",
         subtitle = "MPH assigns 65% to the youngest bin; kernel predicts 0.05%") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(colour = "grey75", fill = NA, linewidth = 0.5),
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, colour = "grey30"),
          legend.position = "bottom")

# ---- B: the pinned components ----
pB <- ggplot(d, aes(bin, var, fill = pinned)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 713.282, linetype = "dashed",
               colour = "#2166AC", linewidth = 0.5) +
    annotate("text", x = 6.5, y = 713.282, label = "713.282",
             vjust = -0.6, size = 3, colour = "#2166AC") +
    scale_y_log10(labels = label_number(scale_cut = cut_short_scale())) +
    scale_fill_manual(values = c(`FALSE` = "grey45", `TRUE` = "#2166AC"),
                      labels = c("estimated", "pinned at floor"), name = NULL) +
    labs(x = "Allele age bin (generations)", y = "Variance component",
         title = "B. Three oldest bins never move",
         subtitle = "Identical to 6 s.f. — a boundary value, not an estimate") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(colour = "grey75", fill = NA, linewidth = 0.5),
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, colour = "grey30"),
          legend.position = "bottom")

p <- pA + pB
ggsave(file.path(FIGS, "mph_failure.png"), p, width = 11, height = 5, dpi = 200)