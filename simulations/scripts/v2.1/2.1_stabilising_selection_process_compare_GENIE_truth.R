library(data.table)
library(ggplot2)
library(stringr)
library(purrr)

f <- glue::glue

# ---------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------
VERSION    <- "2.0"                       # directory version
PREFIX_VER <- "2.1"                       # file-name version (differs!)
DATA       <- f("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v{VERSION}")
REPS_DIR   <- file.path(DATA, "replicates")
FIGS       <- f("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v{VERSION}")
dir.create(FIGS, recursive = TRUE, showWarnings = FALSE)

CONDITIONS <- c("VS_5_NE_25000", "VS_20_NE_25000",
                "VS_100_NE_25000", "VS_10000_NE_25000")
N_REPS <- 100
PRUNED <- FALSE                            # TRUE -> read the ".pruned" GENIE log


# ---------------------------------------------------------------
# Extract one replicate
# ---------------------------------------------------------------
extract_sim_results <- function(condition, rep) {

    rep_dir <- file.path(REPS_DIR, condition, rep)
    prefix  <- f("{PREFIX_VER}_stabilising_selection_{condition}")

    truth_file  <- file.path(rep_dir, f("{prefix}_bin_truth.csv"))
    genie_file  <- file.path(rep_dir, f("{prefix}_neutral_out_GENIE"))
    if (PRUNED) genie_file <- paste0(genie_file, ".pruned")
    pheno_file  <- file.path(rep_dir, f("{prefix}_phenotypes.csv"))

    if (!all(file.exists(truth_file, genie_file, pheno_file))) return(NULL)

    TRUTH     <- fread(truth_file)
    log_lines <- readLines(genie_file, warn = FALSE)

    # --- per-bin h2_g and SE ---
    # "h2_g[3] : 0.0746268 SE : 0.00686682"
    # bins appear twice in the log (Heritabilities, then "overlapping def");
    # keep the first block only.
    h2_pattern <- "^h2_g\\[(\\d+)\\]\\s*:\\s*(-?[0-9.eE+-]+)\\s*SE\\s*:\\s*([0-9.eE+-]+)"
    h2_matches <- str_match(log_lines[str_detect(log_lines, h2_pattern)], h2_pattern)
    if (is.null(h2_matches) || nrow(h2_matches) == 0) return(NULL)

    n_bins <- nrow(TRUTH)
    parsed <- unique(data.table(
        bin  = as.integer(h2_matches[, 2]),
        h2_g = as.numeric(h2_matches[, 3]),
        SE   = as.numeric(h2_matches[, 4])
    ))
    if (nrow(parsed) < n_bins) return(NULL)
    ESTIMATES <- parsed[1:n_bins]

    # --- total h2 ---
    total_line  <- log_lines[str_detect(log_lines, "^Total h2\\s*:")][1]
    if (is.na(total_line)) return(NULL)
    total_match <- str_match(total_line,
        "^Total h2\\s*:\\s*(-?[0-9.eE+-]+)\\s*SE\\s*:\\s*([0-9.eE+-]+)")
    total_h2    <- as.numeric(total_match[2])
    total_h2_SE <- as.numeric(total_match[3])

    # --- phenotypic variance, for the raw-units conversion ---
    V_P <- var(fread(pheno_file)$y)

    # --- join, bin for bin (TRUTH is ordered by bin_lo, matching GENIE 0..n-1) ---
    TRUTH[, bin := 0:(.N - 1)]
    comparison <- merge(TRUTH, ESTIMATES, by = "bin")

    comparison[, `:=`(
        V_estimated = h2_g * V_P,
        SE_raw      = SE * V_P,
        true_share  = V_observed / sum(TRUTH$V_observed),
        est_share   = h2_g / total_h2,
        share_SE    = SE / total_h2,
        condition   = condition,
        V_S         = as.numeric(str_match(condition, "VS_([0-9.]+)_")[, 2]),
        REP         = rep,
        total_h2    = total_h2,
        total_h2_SE = total_h2_SE,
        V_P         = V_P
    )]
    comparison[, `:=`(
        ratio = V_estimated / V_observed,
        z     = (V_estimated - V_observed) / SE_raw,
        bin_label = fifelse(
            is.infinite(bin_hi),
            paste0(scales::comma(bin_lo), "+"),
            paste0(scales::comma(bin_lo), "–", scales::comma(bin_hi)))
    )]
    comparison[]
}


# ---------------------------------------------------------------
# Run over all conditions x reps
# ---------------------------------------------------------------
grid <- CJ(condition = CONDITIONS, rep = 0:(N_REPS - 1), sorted = FALSE)

all_sim_bin_h2 <- rbindlist(
    map2(grid$condition, grid$rep, extract_sim_results),
    use.names = TRUE, fill = TRUE
)

cat(sprintf("Loaded %d replicates across %d conditions\n",
            uniqueN(all_sim_bin_h2[, .(condition, REP)]),
            uniqueN(all_sim_bin_h2$condition)))
print(all_sim_bin_h2[, .(n_reps = uniqueN(REP)), by = condition])

# bin labels ordered by age, taken from the data
bin_levels <- all_sim_bin_h2[, .(bin_lo = bin_lo[1]), by = bin_label][order(bin_lo), bin_label]
all_sim_bin_h2[, bin_label := factor(bin_label, levels = bin_levels)]
all_sim_bin_h2[, condition := factor(condition, levels = CONDITIONS)]
all_sim_bin_h2[, bias := est_share - true_share]

fwrite(all_sim_bin_h2, file.path(DATA, "genie_bin_results_all.csv"))

# one row per replicate: total h2
all_sim_h2 <- unique(all_sim_bin_h2[, .(condition, V_S, REP, total_h2, total_h2_SE)])

bin_summary <- all_sim_bin_h2[, .(
    mean_true = mean(true_share),
    mean_est  = mean(est_share),
    mean_bias = mean(bias),
    se_bias   = sd(bias) / sqrt(.N),
    mad       = mean(abs(bias)),
    rmse      = sqrt(mean(bias^2)),
    n_reps    = .N
), by = .(condition, V_S, bin_label)]
bin_summary[, z := mean_bias / se_bias]
fwrite(bin_summary, file.path(DATA, "genie_bin_summary.csv"))


# ---------------------------------------------------------------
# Plots
# ---------------------------------------------------------------
sfx <- if (PRUNED) "pruned" else "unpruned"

p1 <- ggplot(all_sim_bin_h2, aes(true_share, est_share)) +
    geom_abline(colour = "red", linetype = "dashed") +
    geom_errorbar(aes(ymin = est_share - share_SE,
                      ymax = est_share + share_SE),
                  width = 0, alpha = 0.3) +
    geom_point(aes(colour = bin_label, size = n_variants), alpha = 0.7) +
    facet_wrap(~ condition) +
    viridis::scale_colour_viridis(option = "turbo", discrete = TRUE,
                                  name = "Age bin") +
    labs(x = "True share of V_A", y = "GENIE-estimated share of h2") +
    theme_light()
ggsave(file.path(FIGS, f("genie_vs_truth.{sfx}.png")), p1,
       width = 10, height = 8, dpi = 200)

p2 <- ggplot(all_sim_h2, aes(REP, total_h2)) +
    geom_errorbar(aes(ymin = total_h2 - total_h2_SE,
                      ymax = total_h2 + total_h2_SE), width = 0.3) +
    geom_point() +
    geom_hline(yintercept = 0.5, colour = "red", linetype = "dashed") +
    facet_wrap(~ condition) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Replicate", y = expression(Total~h^2)) +
    theme_light()
ggsave(file.path(FIGS, f("genie_total_h2.{sfx}.png")), p2,
       width = 10, height = 8, dpi = 200)

p3 <- ggplot(bin_summary, aes(bin_label, mean_bias)) +
    geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
    geom_errorbar(aes(ymin = mean_bias - 1.96 * se_bias,
                      ymax = mean_bias + 1.96 * se_bias), width = 0.2) +
    geom_point() +
    facet_wrap(~ condition) +
    labs(x = "Allele age bin (generations)",
         y = "Mean bias (estimated − true share)") +
    theme_light() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
ggsave(file.path(FIGS, f("genie_bias.{sfx}.png")), p3,
       width = 10, height = 8, dpi = 200)

# the profile itself, across selection strengths
p4 <- ggplot(bin_summary, aes(bin_label, group = condition, colour = condition)) +
    geom_line(aes(y = mean_true), linetype = "dashed") +
    geom_line(aes(y = mean_est)) +
    scale_y_log10() +
    labs(x = "Allele age bin (generations)",
         y = "Share of genetic variance",
         caption = "dashed = truth, solid = GENIE") +
    theme_light() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
ggsave(file.path(FIGS, f("genie_profile.{sfx}.png")), p4,
       width = 9, height = 6, dpi = 200)
