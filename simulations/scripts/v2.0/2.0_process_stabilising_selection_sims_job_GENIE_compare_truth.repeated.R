library(data.table)
library(ggplot2)
library(stringr)

f = glue::glue
N_REPS = 10
VERSION = "2.0"
SELECTION_TYPE = "stabilising_selection"
N_E = 20000
V_S = 5
STEM=f("{VERSION}_{SELECTION_TYPE}_VS_{V_S}_NE_{N_E}")


# ---------------------------------------------------------------
# 1. Load truth and raw GENIE output
# ---------------------------------------------------------------
DATA <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v2.0"
FIGS <- "/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v2.0"

extract_sim_results = function(REP) {

    TRUTH <- fread(file.path(DATA, f("rep{REP}"), f("{STEM}_bin_truth.csv")))
    ESTIMATED <- file.path(DATA, f("rep{REP}"), f("{STEM}_out_GENIE"))
    pheno <- fread(file.path(DATA, f("rep{REP}"), f("{STEM}_phenotypes.csv")))

    log_lines <- readLines(ESTIMATED)

    # ---------------------------------------------------------------
    # 2. Parse per-bin h2_g and SE
    #    Lines look like: "h2_g[3] : 0.0746268 SE : 0.00686682"
    #    We want the FIRST block of these (not the duplicated
    #    "overlapping def" block later in the file), so take the
    #    first n_bins matches only.
    # ---------------------------------------------------------------
    h2_pattern <- "^h2_g\\[(\\d+)\\]\\s*:\\s*(-?[0-9.eE+-]+)\\s*SE\\s*:\\s*([0-9.eE+-]+)"

    h2_lines <- log_lines[str_detect(log_lines, h2_pattern)]
    h2_matches <- str_match(h2_lines, h2_pattern)

    n_bins <- nrow(TRUTH)

    parsed <- unique(data.table(
        bin   = as.integer(h2_matches[, 2]),
        h2_g  = as.numeric(h2_matches[, 3]),
        SE    = as.numeric(h2_matches[, 4])
    ))

    # Keep only the first block (bins 0..n_bins-1 appear twice in the log,
    # once under "Heritabilities" and once under "overlapping def")
    ESTIMATES <- parsed[1:n_bins]

    # ---------------------------------------------------------------
    # 3. Parse total h2 and its SE
    # ---------------------------------------------------------------
    total_line <- log_lines[str_detect(log_lines, "^Total h2\\s*:")][1]
    total_match <- str_match(total_line,
    "^Total h2\\s*:\\s*(-?[0-9.eE+-]+)\\s*SE\\s*:\\s*([0-9.eE+-]+)")

    total_h2 <- as.numeric(total_match[2])
    total_h2_SE <- as.numeric(total_match[3])

    cat(sprintf("Total h2 (GENIE): %.4f (SE %.4f)\n", total_h2, total_h2_SE))

    # ---------------------------------------------------------------
    # 4. Parse phenotypic variance from the sample, needed to convert
    #    GENIE's standardised-scale h2_g back to raw trait units.
    #    (Use the value you have from the pipeline; recompute here if
    #    you'd rather read it directly from the phenotype file.)
    # ---------------------------------------------------------------
    pheno <- fread(file.path(DATA, f("rep{REP}"), f("{STEM}_phenotypes.csv")))
    V_P <- var(pheno$y)
    cat(sprintf("Sample phenotypic variance (V_P): %.2f\n", V_P))

    # ---------------------------------------------------------------
    # 5. Join truth and estimates, bin-for-bin
    #    TRUTH is ordered by bin_lo ascending, which matches GENIE's
    #    bin 0..7 ordering (bins were built from the same BINS vector
    #    in the same order in the Python pipeline).
    # ---------------------------------------------------------------
    TRUTH[, bin := 0:(.N - 1)]

    comparison <- merge(TRUTH, ESTIMATES, by = "bin")

    comparison[, `:=`(
    V_estimated  = h2_g * V_P,
    SE_raw       = SE * V_P,
    true_share   = V_observed / sum(TRUTH$V_observed),
    est_share    = h2_g / total_h2
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

    list(
        binned_res = comparison,
        total_h2 = total_h2,
        total_h2_SE = total_h2_SE
    )
}

all_sim_reps = purrr::map(1:N_REPS, extract_sim_results)

all_sim_h2 = rbindlist(purrr::map(all_sim_reps, function(x) data.table(x$total_h2, x$total_h2_SE)))

all_sim_bin_h2 = rbindlist(purrr::map(all_sim_reps, function(x) {
    bin = x$binned_res
        }
    )
)

all_sim_bin_h2[, bin_label := factor(
    bin_label,
    levels = all_sim_bin_h2 %>% arrange(bin_lo) %>% select(bin_label) %>% unique %>% pull(1)
    )
]

all_sim_bin_h2[, bias := est_share - true_share]

bin_summary = all_sim_bin_h2[, .(
  mean_bias = mean(bias),
  se_bias   = sd(bias) / sqrt(.N),
  mad       = mean(abs(bias)),
  rmse      = sqrt(mean(bias^2)),
  n_reps    = .N
), by = bin_label]


bin_summary[, z := mean_bias / se_bias]


p1 = all_sim_bin_h2 %>%
    mutate(bin_label = factor(bin_label, levels = unique(gtools::mixedsort(all_sim_bin_h2$bin_label)))) %>%
    ggplot(aes(x=true_share, y=est_share)) + 
    geom_point(aes(colour=bin_label, size = n_variants)) + 
    geom_errorbar(aes(ymin=est_share - SE, ymax = est_share + SE)) + 
    geom_abline(colour='red', linetype='dashed') + 
    theme_light() +
    viridis::scale_colour_viridis(option = "plasma", discrete=T)

ggsave(file.path(FIGS, f("{STEM}_genie_vs_truth.repeated.pruned.png")), p1, width = 7, height = 6, dpi = 200)

p2 = all_sim_h2 %>% 
    rename(h2 = V1, h2_se = V2) %>% 
    mutate(rep = 1:n()) %>% 
    ggplot(aes(x=rep, y=h2)) + 
    geom_point() +
    geom_errorbar(aes(ymin=h2 - h2_se, ymax = h2 + h2_se), width=0.3) + 
    geom_hline(yintercept = 0.5, colour = 'red', linetype='dashed') + 
    ylim(0, 1)

ggsave(file.path(DATA, "genie_vs_truth.repeated.h2.pruned.png"), p2, width = 7, height = 6, dpi = 200)

p3 = bin_summary %>%
    mutate(bin_label = factor(bin_label, levels = unique(gtools::mixedsort(all_sim_bin_h2$bin_label)))) %>%
    ggplot(aes(x=bin_label, y=mean_bias)) +
    geom_point() +
    geom_errorbar(aes(ymin=mean_bias - se_bias, ymax = mean_bias + se_bias)) +
    geom_hline(yintercept = 0, colour = 'red', linetype = 'dashed') +
    theme_light() +
    theme(axis.text.x = element_text(size=8, angle = 45))

ggsave(file.path(DATA, "bias.repeated.h2.pruned.png"), p3, width = 7, height = 6, dpi = 200)

# ---------------------------------------------------------------
# 6. Plot: true share vs estimated share, with error bars
# ---------------------------------------------------------------
p <- ggplot(comparison, aes(x = true_share, y = est_share)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50") +
  geom_errorbar(
    aes(ymin = est_share - SE / total_h2,
        ymax = est_share + SE / total_h2),
    width = 0, color = "#ca5422"
  ) +
  geom_point(size = 3, color = "#c1440e") +
  geom_text(aes(label = bin_label), hjust = -0.15, vjust = -0.3,
            size = 3, color = "grey40") +
  labs(
    x = "true share of V_A",
    y = "GENIE-estimated share of h2",
    title = "GENIE recovery of the age-stratified variance profile",
    subtitle = sprintf("Total h2: true = %.3f, GENIE = %.3f (SE %.3f)",
                       0.495, total_h2, total_h2_SE)
  ) +
  theme_minimal(base_size = 12) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  coord_cartesian(clip = "off") +
  theme(
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 10, r = 60, b = 10, l = 10)
    )

ggsave(file.path(DATA, "genie_vs_truth.png"), p, width = 7, height = 6, dpi = 200)
print(p)