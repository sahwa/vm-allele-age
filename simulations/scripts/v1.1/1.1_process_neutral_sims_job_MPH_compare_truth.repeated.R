library(data.table)
library(ggplot2)
library(stringr)

f = glue::glue
N_REPS = 10

# ---------------------------------------------------------------
# 1. Load truth and raw GENIE output
# ---------------------------------------------------------------
DATA <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.1"

extract_sim_results = function(REP) {

    TRUTH <- fread(file.path(DATA, f("rep{REP}"), "1.1_bin_truth.csv"))

    ESTIMATED_STD <- fread(file.path(DATA, f("rep{REP}"), "1.1_MPH_std.mq.vc.csv"))
    ESTIMATED_STD[, vc_name := sub(".*_bin_([^_]+)_MPH.*", "\\1", vc_name)]
    ESTIMATED_STD = ESTIMATED_STD[, .(bin = vc_name, var, pve, seP)]
    ESTIMATED_STD_TOTAL = ESTIMATED_STD[bin != "err", sum(pve)]
    ESTIMATED_STD[, est_share_std := pve / ESTIMATED_STD_TOTAL]

    ESTIMATED_UNSTD <- fread(file.path(DATA, f("rep{REP}"), "1.1_MPH_unstd.mq.vc.csv"))
    ESTIMATED_UNSTD[, vc_name := sub(".*_bin_([^_]+)_MPH.*", "\\1", vc_name)]
    ESTIMATED_UNSTD = ESTIMATED_UNSTD[, .(bin = vc_name, var, pve, seP)]
    ESTIMATED_UNSTD_TOTAL = ESTIMATED_UNSTD[bin != "err", sum(pve)]
    ESTIMATED_UNSTD[, est_share_unstd := pve / ESTIMATED_UNSTD_TOTAL]

    PHENO <- fread(file.path(DATA, f("rep{REP}"), "1.1_phenotypes.csv"))
    V_P <- var(PHENO$y)
    cat(sprintf("Sample phenotypic variance (V_P): %.2f\n", V_P))

    TRUTH[, bin := str_c(bin_lo, "-", bin_hi)]
    bin_order <- c("0-100", "100-1000", "1000-10000", "10000-50000", 
               "50000-100000", "100000-200000", "200000-500000", "500000-inf")

    bin_order <- c("0-100", "100-1000", "1000-10000", "10000-50000", 
                "50000-100000", "100000-200000", "200000-500000", "500000-inf")

    ALL_EST <- merge(
        ESTIMATED_STD[bin %in% bin_order, .(bin, est_share_std)],
        ESTIMATED_UNSTD[bin %in% bin_order, .(bin, est_share_unstd)],
        by = "bin"
    )

    comparison <- merge(TRUTH, ALL_EST, by = "bin")
    comparison[, true_share := V_observed / sum(TRUTH$V_observed)]
    comparison[, REP := REP]
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
    levels = c(
        "0–100",
        "100–1,000",
        "1,000–10,000",
        "10,000–50,000",
        "50,000–100,000",
        "100,000–200,000",
        "200,000–500,000",
        "500,000+"
        )
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
    viridis::scale_colour_viridis(option = "turbo", discrete=T)

ggsave(file.path(DATA, "genie_vs_truth.repeated.png"), p1, width = 7, height = 6, dpi = 200)

p2 = all_sim_h2 %>% 
    rename(h2 = V1, h2_se = V2) %>% 
    mutate(rep = 1:n()) %>% 
    ggplot(aes(x=rep, y=h2)) + 
    geom_point() +
    geom_errorbar(aes(ymin=h2 - h2_se, ymax = h2 + h2_se), width=0.3) + 
    geom_hline(yintercept = 0.5, colour = 'red', linetype='dashed') + 
    ylim(0, 1)

ggsave(file.path(DATA, "genie_vs_truth.repeated.h2.png"), p2, width = 7, height = 6, dpi = 200)

p3 = bin_summary %>%
    mutate(bin_label = factor(bin_label, levels = unique(gtools::mixedsort(all_sim_bin_h2$bin_label)))) %>%
    ggplot(aes(x=bin_label, y=mean_bias)) +
    geom_point() +
    geom_errorbar(aes(ymin=mean_bias - se_bias, ymax = mean_bias + se_bias)) +
    geom_hline(yintercept = 0, colour = 'red', linetype = 'dashed') +
    theme_light() +
    theme(axis.text.x = element_text(size=8, angle = 45))

ggsave(file.path(DATA, "bias.repeated.h2.png"), p3, width = 7, height = 6, dpi = 200)

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