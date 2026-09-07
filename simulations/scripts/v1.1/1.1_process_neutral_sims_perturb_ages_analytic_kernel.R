## ---------------------------------------------------------------
## Combine dating error + GENIE estimation error, and compare against
## each contribution measured separately:
##   - GENIE alone      (from the existing pruned/unpruned analysis)
##   - dating error alone (perfect-estimator smearing, prior script)
##   - both together     (this script)
##
## If the combined error is close to (GENIE alone + dating alone),
## the two sources are roughly additive/independent. If it's much
## worse than either alone, they interact - most plausibly because
## GENIE's own bias is concentrated in the same bins the dating bias
## distorts, so errors compound rather than average.
## ---------------------------------------------------------------

library(data.table)
library(stringr)

DATA   <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.1/replicates"
N_REPS <- 100
BINS   <- c(0, 1e2, 1e3, 1e4, 5e4, 1e5, 2e5, 5e5, Inf)

bin_labels <- sapply(seq_len(length(BINS) - 1), function(i) {
    hi <- if (is.infinite(BINS[i + 1])) "inf" else format(BINS[i + 1], scientific = FALSE)
    paste0(format(BINS[i], scientific = FALSE), "-", hi)
})

## ---------------------------------------------------------------
## Reuse the GENIE log parser and persistence fit from earlier work
## ---------------------------------------------------------------
parse_GENIE_LOGFILE <- function(PATH, n_bins) {
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
            h2_g = as.numeric(h2_matches[, 3])
        ))[1:n_bins],
        total_h2 = as.numeric(total_match[2])
    )
}

fit_persistence <- function(d) {
    d <- d[share > 0]
    if (nrow(d) < 3) return(NA_real_)
    pred_share <- function(T) {
        K <- T * (exp(-d$bin_lo / T) -
                  ifelse(is.infinite(d$bin_hi), 0, exp(-d$bin_hi / T)))
        K / sum(K)
    }
    obj <- function(logT) sum((pred_share(exp(logT)) - d$share)^2)
    exp(optimize(obj, c(log(100), log(5e6)))$minimum)
}

bin_edges <- data.table(
    bin    = 0:(length(bin_labels) - 1),
    bin_lo = BINS[-length(BINS)],
    bin_hi = BINS[-1]
)

## ---------------------------------------------------------------
## Ground truth persistence time (no dating error, no GENIE error) -
## computed once per replicate from variant_info.csv directly.
## ---------------------------------------------------------------
truth_T <- function(REP) {
    vi <- fread(file.path(DATA, sprintf("rep%d", REP), "1.1_variant_info.csv"))
    vi <- vi[age > 0 & freq > 0 & freq < 1]
    b <- cut(vi$age, breaks = BINS, labels = bin_labels, right = FALSE)
    dt <- data.table(bin = as.integer(b) - 1,
                     contrib = 2 * vi$freq * (1 - vi$freq) * vi$beta^2)
    agg <- dt[, .(V = sum(contrib)), by = bin]
    agg <- merge(bin_edges, agg, by = "bin", all.x = TRUE)
    agg[is.na(V), V := 0]
    agg[, share := V / sum(V)]
    fit_persistence(agg)
}

## ---------------------------------------------------------------
## GENIE-under-dating-error persistence time, per (REP, rmse)
## ---------------------------------------------------------------
genie_perturbed_T <- function(REP, rmse) {
    log_path <- file.path(DATA, sprintf("rep%d", REP), "dating_error",
                          sprintf("1.1_out_GENIE.rmse%.2f", rmse))
    if (!file.exists(log_path)) return(NA_real_)

    parsed <- tryCatch(parse_GENIE_LOGFILE(log_path, nrow(bin_edges)),
                       error = function(e) NULL)
    if (is.null(parsed)) return(NA_real_)

    d <- merge(bin_edges, parsed$bins, by = "bin")
    d[, share := h2_g / sum(h2_g[h2_g > 0])]
    fit_persistence(d)
}

## ---------------------------------------------------------------
## GENIE-no-dating-error persistence time, per REP (unpruned, from
## the existing v1.1 100-rep GENIE runs)
## ---------------------------------------------------------------
genie_baseline_T <- function(REP) {
    log_path <- file.path(DATA, sprintf("rep%d", REP), "1.1_neutral_out_GENIE")
    if (!file.exists(log_path)) return(NA_real_)

    parsed <- tryCatch(parse_GENIE_LOGFILE(log_path, nrow(bin_edges)),
                       error = function(e) NULL)
    if (is.null(parsed)) return(NA_real_)

    d <- merge(bin_edges, parsed$bins, by = "bin")
    d[, share := h2_g / sum(h2_g[h2_g > 0])]
    fit_persistence(d)
}

## ---------------------------------------------------------------
## Assemble
## ---------------------------------------------------------------
rmse_grid <- c(0.27, 0.30, 0.33)

results <- rbindlist(lapply(1:N_REPS, function(REP) {
    T_true <- truth_T(REP)
    T_genie_only <- genie_baseline_T(REP)

    rbindlist(lapply(rmse_grid, function(r) {
        data.table(
            REP = REP, target_rmse = r,
            T_true = T_true,
            T_genie_only = T_genie_only,
            T_genie_plus_dating = genie_perturbed_T(REP, r)
        )
    }))
}))

results[, rel_err_genie_only := (T_genie_only - T_true) / T_true]
results[, rel_err_combined   := (T_genie_plus_dating - T_true) / T_true]

## ---------------------------------------------------------------
## Compare: is combined error close to additive?
## ---------------------------------------------------------------
summary_tab <- results[, .(
    med_rel_genie_only = median(rel_err_genie_only, na.rm = TRUE),
    med_rel_combined   = median(rel_err_combined,   na.rm = TRUE),
    q10_combined        = quantile(rel_err_combined, 0.1, na.rm = TRUE),
    q90_combined         = quantile(rel_err_combined, 0.9, na.rm = TRUE),
    n_missing            = sum(is.na(T_genie_plus_dating))
), by = target_rmse]

print(summary_tab)

## From the earlier "perfect estimator" smearing experiment, the
## dating-error-alone bias at these RMSE values was roughly:
##   0.27 -> -8.7%,  0.30 -> -8.6%,  0.33 -> -8.4%
## and GENIE-alone (unpruned, obj4) was roughly -0.9%.
## If errors were purely additive, combined would be roughly
## dating_alone + genie_alone (they're both negative here, so
## additive predicts something close to -9.5%). Compare that
## prediction against med_rel_combined above.
dating_alone <- c(`0.27` = -0.0873, `0.30` = -0.0861, `0.33` = -0.0838)
genie_alone  <- -0.009

summary_tab[, additive_prediction := dating_alone[sprintf("%.2f", summary_tab$target_rmse)] + genie_alone]
summary_tab[, interaction := med_rel_combined - additive_prediction]

print(summary_tab)