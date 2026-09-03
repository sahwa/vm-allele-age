library(data.table)
library(ggplot2)
library(stringr)

f = glue::glue
N_REPS = 10
VERSION = "2.0"
SELECTION_TYPE = "stabilising_selection"
N_E = 10000
V_S_vec = c(5, 20, 100)
FIGS="/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs"
GRID = data.table(expand.grid(REP = 1:10, V_S = V_S_vec))

# ---------------------------------------------------------------
# 1. Load truth and raw GENIE output
# ---------------------------------------------------------------
DATA <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v2.0/replicates"
FIGS <- "/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v2.0"

extract_sim_results = function(REP, V_S) {

    STEM=f("{VERSION}_{SELECTION_TYPE}_VS_{V_S}_NE_{N_E}")
    OUT_DIR = file.path(FIGS, STEM)

    REP_VS_PATH = file.path(DATA, f("VS_{V_S}_NE_{N_E}/{REP}"))

    TRUTH <- fread(file.path(REP_VS_PATH, f("{STEM}_bin_truth.csv")))
    ESTIMATED <- file.path(REP_VS_PATH, f("{STEM}_out_GENIE"))
    pheno <- fread(file.path(REP_VS_PATH,  f("{STEM}_phenotypes.csv")))

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
    pheno <- fread(file.path(REP_VS_PATH, f("{STEM}_phenotypes.csv")))
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
        z     = (V_estimated - V_observed) / SE_raw,
        V_S = V_S,
        rep = REP
    )]

    comparison[, bin_label := ifelse(
        is.infinite(bin_hi),
        paste0(scales::comma(bin_lo), "+"),
        paste0(scales::comma(bin_lo), "–", scales::comma(bin_hi))
    )]

    list(
        binned_res = comparison,
        total_h2 = total_h2,
        total_h2_SE = total_h2_SE,
        V_S = V_S,
        N_E = N_E,
        REP = REP
    )
}

all_sim_reps = purrr::pmap(GRID, extract_sim_results)

all_sim_h2 = rbindlist(purrr::map(all_sim_reps, function(x) {
    data.table(x$total_h2, x$total_h2_SE, V_S = x$V_S)
        }
    )
)

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
), by = .(bin_label, V_S)]


bin_summary[, z := mean_bias / se_bias]

p1 = all_sim_bin_h2 %>%  
    ggplot(aes(x=true_share, y=est_share, colour=bin_label, shape=as.factor(V_S))) + 
    geom_point() +
    geom_abline(colour = 'red', linetype = 'dashed') +
    viridis::scale_colour_viridis(discrete=T)
p1_out = file.path(FIGS, f("{VERSION}_{SELECTION_TYPE}_true_v_est_VA_V_S_vary.png"))
ggsave(p1_out, p1)

p2 = all_sim_bin_h2 %>% 
    mutate(midpoint_bin = (bin_hi - bin_lo) / 2) %>% 
    group_by(bin_lo, V_S) %>% summarise(mean_true = mean(true_share)) %>% 
    ggplot(aes(x=bin_lo, y=mean_true, group =  V_S, colour = as.factor(V_S))) + 
    geom_line() +
    geom_point()
p2_out = file.path(FIGS, f("{VERSION}_{SELECTION_TYPE}_true_share_lines_V_S_vary.png"))
ggsave(p2_out, p2)

p3 = bin_summary %>% 
    mutate(V_S = factor(V_S, levels = c(5, 20, 100))) %>%
    ggplot(aes(x=V_S, y=mean_bias, colour=V_S)) + 
    geom_hline(yintercept = 0, colour = 'red', linetype = 'dashed') + 
    geom_jitter() + 
    facet_wrap(~bin_label, nrow=1)
p3_out = file.path(FIGS, f("{VERSION}_{SELECTION_TYPE}_bias_V_S_vary.png"))
ggsave(p3_out, p3)
