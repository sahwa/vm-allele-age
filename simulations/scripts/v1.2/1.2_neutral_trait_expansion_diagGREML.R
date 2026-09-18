f = glue::glue
library(purrr)

VERSION = "1.2.1"
REP = 0

source("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/diagGREML.R")
DATA=f("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v{VERSION}/replicates/rep{REP}")
FIGS = f("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/v{VERSION}")
pheno = fread(file.path(DATA, f("{VERSION}_phenotypes.csv")), header=T)

sets = setdiff(1:6, 5)

N = 2e4

bin_midpoint <- function(x) {
  ifelse(x == "error" | str_detect(x, "\\+$"), NA_real_,
         map_dbl(str_split(x, "-"), ~ mean(as.numeric(.x))))
}

res <- purrr::map(sets, function(nb) {
    C <- fread(file.path(DATA, f("{VERSION}_singleton_counts_nb{nb}.csv")), header = TRUE)
    p_count_ind = melt(C, measure.vars = names(C), variable.name = "bin", value.name = "count") %>%
        ggplot(aes(count + 1)) +
        geom_histogram(bins = 60, fill = "grey35") +
        facet_wrap(~ bin, ncol = 4) +
        scale_x_log10(labels = scales::comma) +
        labs(x = "Singletons per individual (+1, log scale)", y = "Individuals")

    ggsave(file.path(FIGS, sprintf("singleton_count_per_ind_nb%d.png", nb)),
           p_count_ind,  height = 4.5, dpi = 200)

    mean_counts <- sapply(C, mean)
    M_t <- sapply(C, sum)

    A <- Map(function(x, m) x / m, as.list(C), mean_counts)
    X <- matrix(1, nrow = nrow(C), ncol = 1)

    fit <- fit_diagGREML(y = pheno$y, A = A, X = X,
                         constraint = FALSE, magic0316 = TRUE)

    d <- data.table(
        comp = fit$Vlistnames,
        est  = as.numeric(fit$varcmp),
        se   = sqrt(diag(fit$Hi)),
        n_comp = nb
    )
    d <- merge(d, data.table(comp = names(C), M_t = M_t, mean_c = mean_counts),
               by = "comp", all.x = TRUE, sort = FALSE)
    d[, T := bin_midpoint(comp)]
    d[, sigma2_b := est / mean_c]        # == est * N / M_t

    reg <- d[comp != "error" & !is.na(T) & sigma2_b > 0]
    if (nrow(reg) < 3) {
        return(list(fit = d, lm = NULL, Hi = fit$Hi,
                    comps = fit$Vlistnames, n_used = nrow(reg)))
    }

    m <- lm(log(sigma2_b) ~ T, data = reg)
    list(fit = d, lm = m, Hi = fit$Hi, comps = fit$Vlistnames,
         sigma2_m = exp(coef(m)[1]), s_hat = -coef(m)[2], n_used = nrow(reg))
})

summary_dt <- rbindlist(lapply(res, function(r) {
    data.table(
        n_comp   = r$fit$n_comp[1],
        n_used   = if (is.null(r$lm)) 0L else r$n_used,
        sigma2_m = if (is.null(r$lm)) NA_real_ else r$sigma2_m,
        s_hat    = if (is.null(r$lm)) NA_real_ else r$s_hat,
        s_se     = if (is.null(r$lm)) NA_real_ else summary(r$lm)$coef[2, 2],
        s_p      = if (is.null(r$lm)) NA_real_ else summary(r$lm)$coef[2, 4],
        r2       = if (is.null(r$lm)) NA_real_ else summary(r$lm)$r.squared
    )
}))

fits_dt <- rbindlist(lapply(res, `[[`, "fit"))


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


########

library(patchwork)
H2_SING <- 0.10

purrr::iwalk(res, function(r, i) {
    nb <- r$fit$n_comp[1]
    vy <- var(pheno$y)

    d <- as.data.table(r$fit)[comp != "error"]
    d[, comp := factor(comp, levels = gtools::mixedsort(unique(comp)))]
    d[, `:=`(h2 = est / vy, se_h2 = se / vy)]
    d[, `:=`(lo = h2 - 1.96 * se_h2, hi = h2 + 1.96 * se_h2)]

    n_young <- d[comp != "205+", .N]
    d[, truth := fifelse(comp == "205+", NA_real_, H2_SING / n_young)]

    pA <- ggplot(d, aes(comp, h2)) +
        geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
        geom_point(aes(y = truth), shape = 95, size = 8, colour = "#2166AC") +
        geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.5) +
        geom_point(size = 2.6, colour = "#B2182B") +
        scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
        labs(x = "Allele age bin (generations)",
             y = expression("Proportion of"~italic(V[P])),
             title = sprintf("A. Variance components (%d bins)", nb),
             subtitle = sprintf("blue = truth; sum = %.3f vs %.3f",
                                d[comp != "205+", sum(est)], H2_SING)) +
        theme_classic(base_size = 11) +
        theme(panel.border = element_rect(colour = "grey75", fill = NA, linewidth = 0.5),
              axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(face = "bold", size = 11),
              plot.subtitle = element_text(size = 9, colour = "grey30"))

    p <- pA
    w <- 7

    if (!is.null(r$Hi)) {
        comps <- r$comps
        Cor <- cov2cor(r$Hi)
        dimnames(Cor) <- list(comps, comps)
        cd <- as.data.table(as.table(Cor))
        setnames(cd, c("row", "col", "corr"))
        cd[, `:=`(row = factor(row, levels = comps),
                  col = factor(col, levels = rev(comps)))]

        pB <- ggplot(cd, aes(row, col, fill = corr)) +
            geom_tile(colour = "white", linewidth = 0.5) +
            geom_text(aes(label = sprintf("%.2f", corr)), size = 2.7) +
            scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                 midpoint = 0, limits = c(-1, 1), name = "r") +
            labs(x = NULL, y = NULL, title = "B. Sampling correlation") +
            theme_classic(base_size = 11) +
            theme(axis.line = element_blank(), axis.ticks = element_blank(),
                  axis.text.x = element_text(angle = 45, hjust = 1),
                  plot.title = element_text(face = "bold", size = 11))

        p <- pA + pB + plot_layout(widths = c(1, 1.1))
        w <- 12
    }
    ggsave(file.path(FIGS, sprintf("diagGREML_components_nb%d.png", nb)),
           p, width = w, height = 4.5, dpi = 200)
})
