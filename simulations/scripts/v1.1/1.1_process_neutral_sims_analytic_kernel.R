library(purrr)

f = glue::glue
N_REPS = 10

# ---------------------------------------------------------------
# 1. Load truth and raw GENIE output
# ---------------------------------------------------------------
DATA <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.1"

extract_sim_results = function(REP, PRUNED) {

    print(REP); print(PRUNED)

    TRUTH <- fread(f("/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.1/replicates/rep{REP}/1.1_bin_truth.csv"))

    ifelse(
        PRUNED == "PRUNED",
        ESTIMATED <- file.path(f("/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.1/replicates/rep{REP}/1.1_neutral_out_GENIE.pruned")), 
        ESTIMATED <- file.path(f("/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.1/replicates/rep{REP}/1.1_neutral_out_GENIE"))
        )

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
    pheno <- fread(f("/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.1/replicates/rep{REP}/1.1_phenotypes.csv"))
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
    est_share    = h2_g / total_h2,
    REP          = REP,
    TYPE         = PRUNED
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

params = expand.grid(REP = 1:100, PRUNED = c("PRUNED", "UNPRUNED"))

sim_res = purrr::pmap(params, extract_sim_results)

all_sim_bin_h2 = rbindlist(purrr::map(sim_res, function(x) purrr::pluck(x, "binned_res")))

fit_persistence <- function(d) {
    setDT(d)
    d <- d[est_share > 0]
    if (nrow(d) < 3) return(NA_real_)

    # integral of exp(-t/T) across each bin, normalised to a share
    pred_share <- function(T) {
        lo <- d$bin_lo; hi <- d$bin_hi
        K <- T * (exp(-lo / T) - ifelse(is.infinite(hi), 0, exp(-hi / T)))
        K / sum(K)
    }

    obj1 <- function(logT) sum((log(pred_share(exp(logT))) - log(d$est_share))^2)
    obj2 <- function(logT) sum(d$est_share * (log(pred_share(exp(logT))) - log(d$est_share))^2)
    obj3 <- function(logT) sum(d$n_eff * (log(pred_share(exp(logT))) - log(d$est_share))^2)
    obj4 <- function(logT) sum((pred_share(exp(logT)) - d$est_share)^2)


    data.table(
        o_1 = exp(optimize(obj1, c(log(100), log(5e6)))$minimum),
        o_2 = exp(optimize(obj2, c(log(100), log(5e6)))$minimum),
        o_3 = exp(optimize(obj3, c(log(100), log(5e6)))$minimum),
        o_4 = exp(optimize(obj4, c(log(100), log(5e6)))$minimum)
        )
}

pers <- all_sim_bin_h2[, {
    est <- fit_persistence(.SD)
    tru <- fit_persistence(copy(.SD)[, est_share := true_share])

    c(
        setNames(est, paste0("T_est_", names(est))),
        setNames(tru, paste0("T_true_", names(tru)))
    )
}, by = .(REP, TYPE)]


pers_long <- melt(
    pers,
    id.vars = c("REP", "TYPE"),
    measure.vars = patterns(
        T_est  = "^T_est_",
        T_true = "^T_true_"
    ),
    variable.name = "obj"
)

pers_summary <- pers_long[, .(
    med_T_true = median(T_true, na.rm = TRUE),
    iqr_T_true = IQR(T_true, na.rm = TRUE),
    med_T_est  = median(T_est, na.rm = TRUE),
    med_rel = median(
        (T_est - T_true) / T_true,
        na.rm = TRUE
    ),
    q10 = quantile(
        (T_est - T_true) / T_true,
        0.1,
        na.rm = TRUE
    ),
    q90 = quantile(
        (T_est - T_true) / T_true,
        0.9,
        na.rm = TRUE
    )
), by = .(TYPE, obj)]
