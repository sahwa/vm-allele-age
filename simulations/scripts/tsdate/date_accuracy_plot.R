## ---------------------------------------------------------------------------
## shrinkage_fig.R -- tsdate shrinks allele ages toward the prior.
##
## Reads bias_by_age_and_mac.csv (written by date_accuracy_dt.R) so it runs
## against the full table rather than a transcription.
##
## Three panels:
##   A  inferred vs true age        -- the collapse onto the prior
##   B  bias fold vs true age       -- the pivot, explicit
##   C  age distribution at MAC=1   -- what it costs the estimator
## ---------------------------------------------------------------------------

library(data.table); library(ggplot2); library(patchwork)

IN  <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/tsdate/date_accuracy/bias_by_age_and_mac.csv"
OUT <- "tsdate_shrinkage.png"
MIN_N <- 1000                      # drop cells too thin to read

## -- palette (dataviz reference instance) ------------------------------------
SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"
GRID <- "#dedcd6"; REFLINE <- "#9a978f"; ORANGE <- "#eb6834"
BLUES <- c("#86b6ef","#6da7ec","#5598e7","#3987e5","#2a78d6",
           "#256abf","#1c5cab","#184f95","#104281")   # ordinal, >= step 250

AGE_EDGES <- c(0,1,2,4,8,16,32,64,128,256,512,1024,2048)
T_REP     <- c(0.5, sqrt(AGE_EDGES[2:12] * AGE_EDGES[3:13]))
MAC_LAB   <- c("1","2-3","4-7","8-15","16-31","32-63","64-127","128-255","256+")

theme_v <- function() theme_minimal(base_size = 11) + theme(
  plot.background  = element_rect(fill = SURFACE, colour = NA),
  panel.background = element_rect(fill = SURFACE, colour = NA),
  panel.grid.major = element_line(colour = GRID, linewidth = .3),
  panel.grid.minor = element_blank(),
  axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK),
  plot.title = element_text(colour = INK, size = 12.5, hjust = 0,
                            margin = margin(b = 8)),
  legend.position = "none")

## -- data --------------------------------------------------------------------
dt <- fread(IN)[n >= MIN_N]
dt[, `:=`(t_true = T_REP[bin_true],
          mac_f  = factor(MAC_LAB[mac_str], levels = MAC_LAB))]
dt[, t_inf := t_true * bias_fold]
setorder(dt, mac_str, bin_true)

## emphasise MAC=1 (current design) and MAC 16-31 (B ~ 30,000)
dt[, lw := fifelse(mac_str %in% c(1L, 5L), 1.1, 0.6)]

## pivot: where the MAC=1 bias curve crosses 1, interpolated in log-log
p1 <- dt[mac_str == 1][order(t_true)]
i  <- which(diff(sign(log(p1$bias_fold))) != 0)[1]
PIV <- exp(with(p1, log(t_true)[i] - log(bias_fold)[i] *
                 diff(log(t_true)[i:(i+1)]) / diff(log(bias_fold)[i:(i+1)])))
cat(sprintf("pivot %.2f generations\n", PIV))

LOGBRK <- c(1, 10, 100, 1000)

## -- A: inferred vs true -----------------------------------------------------
pA <- ggplot(dt, aes(t_true, t_inf, colour = mac_f)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22",
              colour = REFLINE, linewidth = .5) +
  geom_hline(yintercept = PIV, linetype = "12", colour = ORANGE, linewidth = .5) +
  geom_line(aes(linewidth = I(lw))) +
  geom_point(size = 1.3) +
  annotate("text", x = 2300, y = PIV * 1.25, hjust = 1, vjust = 0,
           label = sprintf("prior  %.0f gen", PIV), colour = ORANGE, size = 3.2) +
  annotate("text", x = 2.2, y = 0.75, hjust = 0, size = 3.3, colour = INK,
           label = "MAC = 1:\n10,000-fold range of true age\n→ 1.4-fold range inferred") +
  scale_colour_manual(values = BLUES, drop = FALSE) +
  scale_x_log10(breaks = LOGBRK, limits = c(.3, 2600)) +
  scale_y_log10(breaks = LOGBRK, limits = c(.3, 2600)) +
  labs(  x = "true allele age (generations)", y = "tsdate age (generations)") +
  theme_v() +
  theme(legend.position = c(.02, .98), legend.justification = c(0, 1),
        legend.title = element_text(size = 8.5, colour = INK2),
        legend.text = element_text(size = 8), legend.key.height = unit(9, "pt"),
        legend.background = element_blank()) +
  guides(colour = guide_legend(title = "minor allele count", ncol = 2,
                               override.aes = list(linewidth = 1.1, size = 1.3)))

## -- B: bias fold ------------------------------------------------------------
pB <- ggplot(dt, aes(t_true, bias_fold, colour = mac_f)) +
  geom_hline(yintercept = 1, linetype = "22", colour = REFLINE, linewidth = .5) +
  geom_vline(xintercept = PIV, linetype = "12", colour = ORANGE, linewidth = .5) +
  geom_line(aes(linewidth = I(lw))) +
  geom_point(size = 1.3) +
  annotate("point", x = PIV, y = 1, size = 3, colour = ORANGE) +
  annotate("text", x = PIV * 1.35, y = 13, hjust = 0, vjust = 1, size = 3.3,
           colour = ORANGE, label = sprintf("pivot\n%.0f gen", PIV)) +
  annotate("text", x = .35, y = 8.5, hjust = 0, size = 3.3, colour = INK2,
           fontface = "italic", label = "too old") +
  annotate("text", x = .35, y = .52, hjust = 0, size = 3.3, colour = INK2,
           fontface = "italic", label = "too young") +
  annotate("text", x = 1.1, y = 12, hjust = 0, size = 3.3, colour = BLUES[1],
           fontface = "bold", label = "MAC = 1") +
  annotate("text", x = 300, y = 1.12, size = 3.3, colour = BLUES[9],
           label = "MAC 256+") +
  scale_colour_manual(values = BLUES, drop = FALSE) +
  scale_x_log10(breaks = LOGBRK, limits = c(.3, 2600)) +
  scale_y_log10(breaks = c(.5, 1, 2, 4, 8, 16),
                labels = c("0.5×","1×","2×","4×","8×","16×"),
                limits = c(.45, 16)) +
  labs( x = "true allele age (generations)", y = "tsdate age / true age") +
  theme_v()

## -- C: the age distribution at MAC = 1 --------------------------------------
d1 <- dt[mac_str == 1][, w := n / sum(n)][
  , `:=`(lo = log(pmax(AGE_EDGES[bin_true], .25)),
         hi = log(AGE_EDGES[bin_true + 1]),
         mu = log(t_inf), sdl = log(scatter))]

g <- seq(log(.25), log(2600), length.out = 700)
dens <- rbindlist(list(
  data.table(age = exp(g), src = "true age",
             d = vapply(g, function(x)
               sum(d1$w * (x >= d1$lo & x < d1$hi) / (d1$hi - d1$lo)), 0)),
  data.table(age = exp(g), src = "tsdate age",
             d = vapply(g, function(x)
               sum(d1$w * dnorm(x, d1$mu, d1$sdl)), 0))))
dens[, src := factor(src, levels = c("true age", "tsdate age"))]

mid  <- (d1$lo + d1$hi) / 2
mt   <- sum(d1$w * mid);       vt   <- sum(d1$w * (mid - mt)^2)
mi   <- sum(d1$w * d1$mu);     vi_b <- sum(d1$w * (d1$mu - mi)^2)
vi   <- vi_b + sum(d1$w * d1$sdl^2)
rho  <- sum(d1$w * (mid - mt) * (d1$mu - mi)) / sqrt(vt * vi)

pC <- ggplot(dens, aes(age, d, colour = src, fill = src)) +
  geom_area(position = "identity", alpha = .15, colour = NA) +
  geom_line(linewidth = .9) +
  annotate("label", x = 2400, y = Inf, hjust = 1, vjust = 1, size = 2.9,
           colour = INK2, family = "mono", label.size = 0, fill = SURFACE,
           lineheight = 1.4, label = sprintf(
             "SD of log age\n   true             %.2f\n   tsdate           %.2f\n   tsdate, signal   %.2f\n\nr(log ages)         %.2f",
             sqrt(vt), sqrt(vi), sqrt(vi_b), rho)) +
  scale_colour_manual(values = c("true age" = INK, "tsdate age" = ORANGE)) +
  scale_fill_manual(values  = c("true age" = INK, "tsdate age" = ORANGE)) +
  scale_x_log10(breaks = LOGBRK, limits = c(.3, 2600)) +
  labs(x = "allele age (generations)", y = "density of singletons") +
  theme_v() +
  theme(legend.position = c(.02, .98), legend.justification = c(0, 1),
        legend.title = element_blank(), legend.text = element_text(size = 9),
        legend.background = element_blank())

## -- assemble ----------------------------------------------------------------
fig <- (pA | pB | pC) + plot_annotation(
  subtitle = paste("msprime n=500,000, L=1e8, 40.5M variants. True topology, so this is a",
                   "lower bound on error — tsinfer adds topology error on top."),
  theme = theme(plot.background = element_rect(fill = SURFACE, colour = NA),
                plot.title = element_text(size = 15, face = "bold", colour = INK),
                plot.subtitle = element_text(size = 10, colour = INK2,
                                             margin = margin(b = 10))))

ggsave(OUT, fig, width = 15.5, height = 5.3, dpi = 200, bg = SURFACE)
cat("written:", OUT, "\n")
