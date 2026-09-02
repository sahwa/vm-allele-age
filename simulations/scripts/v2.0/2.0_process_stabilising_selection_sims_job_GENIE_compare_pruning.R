library(data.table)
library(ggplot2)
library(stringr)

f = glue::glue
N_REPS = 100
VERSION = "1.1"
SELECTION_TYPE = "neutral"
N_E = 20000
V_S = 5
STEM = f("{VERSION}_{SELECTION_TYPE}")

DATA <- f("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v{VERSION}")
FIGS <- f("/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v{VERSION}")

# ---------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------
do_comparison_stuff = function(comparison, total_h2, V_P) {
    comparison[, `:=`(
        V_estimated = h2_g * V_P,
        SE_raw      = SE * V_P,
        true_share  = V_observed / sum(V_observed),
        est_share   = h2_g / total_h2
    )]

    comparison[, `:=`(
        ratio = V_estimated / V_observed,
        z     = (V_estimated - V_observed) / SE_raw
    )]

    comparison[, bin_label := ifelse(
        is.infinite(bin_hi),
        paste0(scales::comma(bin_lo), "+"),
        paste0(scales::comma(bin_lo), "–", scales::comma(bin_hi))
    )]

    comparison
}

parse_GENIE_LOGFILE = function(PATH, n_bins) {
    log_lines <- readLines(PATH)

    h2_pattern <- "^h2_g\\[(\\d+)\\]\\s*:\\s*(-?[0-9.eE+-]+)\\s*SE\\s*:\\s*([0-9.eE+-]+)"
    h2_lines   <- log_lines[str_detect(log_lines, h2_pattern)]
    h2_matches <- str_match(h2_lines, h2_pattern)

    total_line  <- log_lines[str_detect(log_lines, "^Total h2\\s*:")][1]
    total_match <- str_match(total_line,
        "^Total h2\\s*:\\s*(-?[0-9.eE+-]+)\\s*SE\\s*:\\s*([0-9.eE+-]+)")

    list(
        bins = unique(data.table(
            bin  = as.integer(h2_matches[, 2]),
            h2_g = as.numeric(h2_matches[, 3]),
            SE   = as.numeric(h2_matches[, 4])
        ))[1:n_bins],
        total_h2    = as.numeric(total_match[2]),
        total_h2_SE = as.numeric(total_match[3])
    )
}

# ---------------------------------------------------------------
# Per-replicate extraction: returns one long table with both
# pruned and unpruned results stacked
# ---------------------------------------------------------------
extract_sim_results = function(REP) {

    print(REP)

    TRUTH <- fread(file.path(DATA, "replicates", f("rep{REP}"), f("{VERSION}_bin_truth.csv")))
    TRUTH[, bin := 0:(.N - 1)]
    n_bins <- nrow(TRUTH)

    pheno <- fread(file.path(DATA, "replicates", f("rep{REP}"), f("{VERSION}_phenotypes.csv")))
    V_P   <- var(pheno$y)

    parsed_PRUNED   <- parse_GENIE_LOGFILE(
        file.path(DATA, "replicates", f("rep{REP}"),  f("{STEM}_out_GENIE.pruned")), n_bins)
    parsed_UNPRUNED <- parse_GENIE_LOGFILE(
        file.path(DATA, "replicates", f("rep{REP}"),  f("{STEM}_out_GENIE")), n_bins)

    cmp_PRUNED <- do_comparison_stuff(
        merge(copy(TRUTH), parsed_PRUNED$bins, by = "bin"),
        parsed_PRUNED$total_h2, V_P)
    cmp_PRUNED[, `:=`(REP = REP, pruned = TRUE,
                      total_h2 = parsed_PRUNED$total_h2,
                      total_h2_SE = parsed_PRUNED$total_h2_SE)]

    cmp_UNPRUNED <- do_comparison_stuff(
        merge(copy(TRUTH), parsed_UNPRUNED$bins, by = "bin"),
        parsed_UNPRUNED$total_h2, V_P)
    cmp_UNPRUNED[, `:=`(REP = REP, pruned = FALSE,
                        total_h2 = parsed_UNPRUNED$total_h2,
                        total_h2_SE = parsed_UNPRUNED$total_h2_SE)]

    rbind(cmp_PRUNED, cmp_UNPRUNED)
}

# ---------------------------------------------------------------
# Run across replicates and stack
# ---------------------------------------------------------------
all_sim_bin_h2 <- rbindlist(purrr::map(1:N_REPS, extract_sim_results))

all_sim_bin_h2[, bias := est_share - true_share]

bin_levels <- c("0–100", "100–1,000", "1,000–10,000", "10,000–50,000",
                "50,000–100,000", "100,000–200,000", "200,000–500,000",
                "500,000+")
all_sim_bin_h2[, bin_label := factor(bin_label, levels = bin_levels)]

# ---------------------------------------------------------------
# Paired analysis: within-replicate difference in bias
# ---------------------------------------------------------------
paired <- dcast(all_sim_bin_h2, bin_label + REP ~ pruned, value.var = "bias")
setnames(paired, c("FALSE", "TRUE"), c("bias_unpruned", "bias_pruned"))

paired[, delta := bias_pruned - bias_unpruned]

paired_summary <- paired[, .(
    mean_delta = mean(delta),
    se_delta   = sd(delta) / sqrt(.N),
    n          = .N
), by = bin_label]

paired_summary[, t_stat := mean_delta / se_delta]
paired_summary[, p_val  := 2 * pt(-abs(t_stat), df = n - 1)]

print(paired_summary)

###### plots ########

library(dplyr)

p1 <- ggplot(all_sim_bin_h2, aes(x = bin_label, y = bias, colour = pruned)) +
    geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed") +
    geom_point(
        position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.75),
        alpha = 0.25, size = 1
    ) +
    stat_summary(
        fun = mean, fun.min = ~mean(.x) - sd(.x)/sqrt(length(.x)),
        fun.max = ~mean(.x) + sd(.x)/sqrt(length(.x)),
        geom = "pointrange",
        position = position_dodge(width = 0.75),
        size = 0.5, linewidth = 0.9
    ) +
    scale_colour_manual(
        values = c("FALSE" = "#d95f02", "TRUE" = "#1b9e77"),
        labels = c("FALSE" = "Unpruned", "TRUE" = "Pruned"),
        name = NULL
    ) +
    labs(x = "Allele age bin (generations)",
         y = "Bias (estimated − true share of V_A)") +
    theme_light(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",
          panel.grid.minor = element_blank())

ggsave(file.path(FIGS, "bias_n100_jitter.png"), p1, width = 9, height = 5.5, dpi = 200)


p4 <- ggplot(all_sim_bin_h2, aes(x = true_share, y = est_share, colour = bin_label)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(alpha = 0.2, size = 1) +
    stat_summary(fun = mean, geom = "point", size = 2.5, colour = "black") +
    facet_wrap(~pruned, labeller = as_labeller(
        c("FALSE" = "Unpruned", "TRUE" = "Pruned"))) +
    scale_colour_viridis(discrete=T, option = "plasma", begin = 0.1, end = 0.9,
                            name = "Age bin") +
    labs(x = "True share of V_A", y = "Estimated share of V_A") +
    coord_fixed() +
    theme_light(base_size = 11) +
    theme(panel.grid.minor = element_blank())

ggsave(file.path(FIGS, "obs_vs_exp_faceted_n100.png"), p4, width = 11, height = 5.5, dpi = 200)

p5 <- ggplot(all_sim_bin_h2[est_share > 0],
             aes(x = true_share, y = est_share, colour = bin_label)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(alpha = 0.2, size = 1) +
    facet_wrap(~pruned, labeller = as_labeller(
        c("FALSE" = "Unpruned", "TRUE" = "Pruned"))) +
    scale_x_log10() + scale_y_log10() +
    scale_colour_viridis_d(option = "plasma", begin = 0.1, end = 0.9,
                            name = "Age bin") +
    labs(x = "True share of V_A (log)", y = "Estimated share of V_A (log)") +
    coord_fixed() +
    theme_light(base_size = 11)

ggsave(file.path(FIGS, "obs_vs_exp_faceted_n100.png"), p4, width = 11, height = 5.5, dpi = 200)
