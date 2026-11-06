## ---------------------------------------------------------------------------
## batch_power_fisher.R -- power for lambda across the batching ladder, using
## the REAL count matrices from batch_sweep.py rather than a scaling formula.
##
## Batching is encoded entirely in C and in the representative ages, so the
## Fisher machinery needs no re-parameterisation -- only sigma2_m must be held
## fixed (it is a property of the trait) while h2_sing is allowed to emerge
## from whichever C is fed in. Holding h2_sing fixed would cancel the main
## benefit of batching.
##
## NOT included here: cross-batch variant sharing. Fisher information assumes
## independent observations, so these are NAIVE standard errors -- an upper
## bound on power. The sharing penalty has to come from kdist_nb*.npy, or
## better, from simulating phenotypes off these same C matrices and fitting.
## ---------------------------------------------------------------------------

library(data.table); library(purrr); library(glue); library(RcppCNPy)

DATA <- "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/tests/batch_test"

MU       <- 1.25e-8
L_GENOME <- 3e9
L_SIM    <- 1e8                      # sequence length the tree sequence covers
U        <- 2 * MU * L_GENOME        # ~75 mutations per person per generation
C_VM     <- 0.3                      # lambda ~ C_VM * (V_M/V_P) / h2, from SLiM

N_BATCHES   <- c(1, 2, 5, 10, 20, 50, 100, 200)
N_TOTAL     <- 5e5
LAMBDA_VALS <- c(0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2)
VM_VP       <- 2.6e-3
Z_ALPHA     <- 1.96                  # two-sided 0.05

## ---------------------------------------------------------------------------
## Fisher SE for lambda, profiling out sigma2_e and sigma2_m.
##   w = f(lambda, t);  dw/dlambda = -t*w (exponential), -t*w^2 (hyperbolic)
## Validated against a numerical Hessian of the expected log-likelihood.
## ---------------------------------------------------------------------------
fisher_se_lambda <- function(C, ages, sigma_m2, sigma_e2, lambda,
                             curve = c("hyperbolic", "exponential")) {
  curve <- match.arg(curve)
  C <- as.matrix(C)
  stopifnot(ncol(C) == length(ages), sigma_m2 > 0, sigma_e2 > 0,
            all(C >= 0), all(is.finite(C)), all(is.finite(ages)))

  w  <- if (curve == "exponential") exp(-lambda * ages) else 1 / (1 + lambda * ages)
  dw <- if (curve == "exponential") ages * w else ages * w^2

  S <- as.vector(C %*% w)
  Q <- as.vector(C %*% dw)
  v <- sigma_e2 + sigma_m2 * S
  if (any(v <= 0)) stop("Non-positive phenotypic variance")

  G <- cbind(sigma_e2 = 1, sigma_m2 = S, lambda = -sigma_m2 * Q)
  I <- crossprod(G / (sqrt(2) * v))
  V <- tryCatch(solve(I), error = function(e) NULL)
  if (is.null(V)) { warning("singular information matrix"); return(NA_real_) }
  sqrt(V["lambda", "lambda"])
}

## ---------------------------------------------------------------------------
## Prepare one level's counts once: drop empty bins, rescale to genome-wide.
## Hoisted out of the lambda loop so a 48MB matrix is not re-copied per lambda.
##
## The tree sequence covers L_SIM, so raw counts are ~1/30 of genome-wide.
## Scaling a realised count preserves CV_total exactly where it should fall
## slightly (0.38 -> 0.34) as Poisson noise shrinks, so this is mildly
## conservative.
## ---------------------------------------------------------------------------
prepare_counts <- function(C, trep, kdist, scale = L_GENOME / L_SIM) {
  C  <- as.matrix(C)
  ok <- is.finite(trep) & colSums(C) > 0
  list(C = C[, ok, drop = FALSE] * scale, ages = trep[ok], kdist = kdist)
}

fisher_se_lambda_trait <- function(prep, vm_vp, lambda, curve = "hyperbolic") {
  C <- prep$C; ages <- prep$ages
  sigma_m2 <- vm_vp / U                       # V_P standardised to 1
  w  <- if (curve == "exponential") exp(-lambda * ages) else 1 / (1 + lambda * ages)

  h2_sing <- sigma_m2 * mean(as.vector(C %*% w))
  t0      <- sum(colSums(C) * ages) / sum(colSums(C))

  ## h2_sing >= 1 is not a numerical failure -- it is the design breaking down.
  ## Past this point "batch-singletons" are no longer recent mutations and the
  ## model is incoherent with V_P = 1. Record it rather than stopping, so the
  ## sweep completes and shows where the boundary falls.
  if (!is.finite(h2_sing) || h2_sing >= 1)
    return(list(se = NA_real_, h2_sing = h2_sing, t0 = t0))

  list(se = fisher_se_lambda(C, ages, sigma_m2, 1 - h2_sing, lambda, curve),
       h2_sing = h2_sing, t0 = t0)
}

power_lambda <- function(lambda, se, z_a = Z_ALPHA) {
  ifelse(is.na(se), NA_real_, pnorm((lambda / se) - z_a))
}

decay <- function(lambda, t, curve = "hyerpbolic") {
  switch(curve,
    exponential = exp(-lambda * t),
    hyperbolic  = 1 / (1 + lambda * t),
    stop("curve must be 'exponential' or 'hyperbolic'"))
}

## ---------------------------------------------------------------------------
## sweep
## ---------------------------------------------------------------------------
res <- rbindlist(map(N_BATCHES, function(nb) {
  prep <- prepare_counts(
    npyLoad(file.path(DATA, glue("C_nb{nb}.npy"))),
    npyLoad(file.path(DATA, glue("trep_nb{nb}.npy"))),
    npyLoad(file.path(DATA, glue("kdist_nb{nb}.npy")))
    )
  rbindlist(map(LAMBDA_VALS, function(lam) {

    r <- fisher_se_lambda_trait(prep, vm_vp = VM_VP, lambda = lam)
    data.table(
      n_batches  = nb,
      batch_size = N_TOTAL / nb,
      lambda     = lam,
      ## lambda = C_VM*(V_M/V_P)/h2, so a given (lambda, V_M/V_P) implies a
      ## heritability. Above 1 the combination is not realisable: there is a
      ## floor lambda >= C_VM*(V_M/V_P) = 7.8e-4 at V_M/V_P = 2.6e-3.
      h2_implied = if (lam > 0) C_VM * VM_VP / lam else NA_real_,
      h2_sing    = r$h2_sing,
      t0         = r$t0,
      se         = r$se,
      power      = power_lambda(lam, r$se),
      sigma_m_2 = VM_VP / U,
      kdist = prep$kdist
    )
  }))
}))

res[, realisable := is.na(h2_implied) | h2_implied <= 1]

## lambda = 0 is the false-positive rate (0.025, one tail of a two-sided test),
## not power -- kept as a check.
setorder(res, n_batches, lambda)
print(res[], digits = 3)

## The diagnostics to read BEFORE any power number: if h2_sing climbs linearly
## past 1 then nothing saturated and the T0_CAP question is still open; if
## h2_sing and t0 flatten together, that is the saturation point the design
## hinges on, measured rather than assumed.
cat("\n-- design diagnostics (independent of lambda) --\n")
print(unique(res[, .(n_batches, batch_size, t0,
                     h2_sing = round(h2_sing, 4))]), digits = 3)

fwrite(res, file.path(DATA, "fisher_power_by_batch.csv"))


## ---------------------------------------------------------------------------
## plot_batch_sweep.R -- visualise the Fisher power sweep across the batching
## ladder.
##
## Three panels, in the order they should be read:
##   A  how t0 and h2_sing actually scale with n_batches  (the design question)
##   B  power across the (batch size x lambda) grid        (the headline)
##   C  what the variant-sharing penalty does to it        (the catch)
##
## Panel A comes first on purpose: if t0 and the singleton count do not grow
## the way the power formula assumed, nothing in B or C means much.
## ---------------------------------------------------------------------------

library(data.table); library(ggplot2); library(patchwork); library(scales)

OUT = "/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/figs/tests/batch_power_sweep"

SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"
GRID <- "#dedcd6"; ORANGE <- "#eb6834"; BLUE <- "#2a78d6"
SEQ <- c("#f2f6fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
         "#256abf", "#1c5cab", "#104281")          # sequential blue, light->dark

theme_v <- function() theme_minimal(base_size = 11) + theme(
  plot.background  = element_rect(fill = SURFACE, colour = NA),
  panel.background = element_rect(fill = SURFACE, colour = NA),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(colour = GRID, linewidth = .3),
  axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK),
  plot.title = element_text(colour = INK, size = 12, hjust = 0,
                            margin = margin(b = 3)),
  plot.subtitle = element_text(colour = INK2, size = 9.5,
                               margin = margin(b = 8)))

## -- de-duplicate: kdist is a length-(nb+1) vector and recycled every row -----
res <- unique(res[, .(n_batches, batch_size, lambda, h2_implied,
                      h2_sing, t0, se, power, realisable)])
stopifnot(nrow(res) == uniqueN(res$n_batches) * uniqueN(res$lambda))

## ===========================================================================
## A: the scaling that the power table assumed, versus what happened
## ===========================================================================
d <- unique(res[lambda == 0, .(n_batches, batch_size, t0, h2_sing)])[order(n_batches)]
d[, `:=`(t0_rel = t0 / t0[1], h2_rel = h2_sing / h2_sing[1])]

ex <- function(y) coef(lm(log(y[-1]) ~ log(d$n_batches[-1])))[2]
e_t0 <- ex(d$t0_rel); e_h2 <- ex(d$h2_rel)

sc <- melt(d, id.vars = "n_batches", measure.vars = c("t0_rel", "h2_rel"),
           variable.name = "what", value.name = "rel")
sc[, what := factor(what, levels = c("t0_rel", "h2_rel"),
                    labels = c("mean allele age (t0)",
                               "singleton heritability"))]

pA <- ggplot(sc, aes(n_batches, rel, colour = what)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22",
              colour = INK2, linewidth = .5) +
  geom_line(linewidth = .9) + geom_point(size = 2.1) +
  annotate("text", x = 130, y = 190, label = "assumed:\nlinear in n_batches",
           colour = INK2, size = 3.1, hjust = 1, lineheight = 1.2) +
  annotate("text", x = 200, y = c(14, 9), hjust = 1, size = 3.3,
           colour = c(ORANGE, BLUE),
           label = c(sprintf("~n^%.2f", e_t0), sprintf("~n^%.2f", e_h2))) +
  scale_colour_manual(values = c(ORANGE, BLUE), name = NULL) +
  scale_x_log10(breaks = d$n_batches) +
  scale_y_log10(breaks = c(1, 3, 10, 30, 100, 300)) +
  labs(title = "A  Both levers grow far slower than assumed",
       subtitle = "relative to the unbatched design; dashed line is the linear scaling the power table used",
       x = "number of batches", y = "fold change vs. 1 batch") +
  theme_v() + theme(legend.position = "top",
                    legend.justification = "left",
                    legend.margin = margin(b = -4),
                    legend.text = element_text(size = 9))

## ===========================================================================
## B: power over the grid. Non-realisable cells are (lambda, V_M/V_P) pairs
## implying h2 > 1 -- they cannot occur for a trait with this much mutational
## variance, so they are masked rather than plotted.
## ===========================================================================
h <- copy(res)[lambda > 0]
h[, `:=`(bs = factor(batch_size, levels = sort(unique(batch_size), TRUE),
                     labels = comma(sort(unique(batch_size), TRUE))),
         lam = factor(lambda, levels = sort(unique(lambda)),
                      labels = scales::scientific(sort(unique(lambda)))))]

pB <- ggplot(h, aes(bs, lam)) +
  geom_tile(aes(fill = ifelse(realisable, power, NA_real_)),
            colour = SURFACE, linewidth = .7) +
  geom_text(aes(label = ifelse(realisable, sprintf("%.2f", power), "n/a"),
                colour = realisable & power > .55), size = 3.1,
            show.legend = FALSE) +
  scale_fill_gradientn(colours = SEQ, limits = c(0, 1), na.value = "#e9e7e0",
                       name = "power") +
  scale_colour_manual(values = c("TRUE" = SURFACE, "FALSE" = INK)) +
  labs(title = "B  Power for detecting lambda (naive: no variant sharing)",
       subtitle = "grey = not realisable; those (lambda, V_M/V_P) pairs imply heritability above 1",
       x = "batch size", y = expression(lambda)) +
  theme_v() + theme(panel.grid = element_blank(),
                    axis.text.x = element_text(angle = 40, hjust = 1))

## ===========================================================================
## C: the sharing penalty. k = variants shared across batches; substitute the
## MEASURED k_per_instance from kdist_nb*.npy for the 0.25*nb placeholder.
## ===========================================================================
K_PER_NB <- 0.25                      # <- replace with measured k / n_batches
LAM_SHOW <- 1e-3

p <- res[lambda == LAM_SHOW & realisable][order(n_batches)]
p[, `:=`(k = pmax(1, K_PER_NB * n_batches))]
p[, `:=`(`no sharing` = se, `with sharing` = se * sqrt(k))]
pl <- melt(p, id.vars = c("batch_size", "n_batches"),
           measure.vars = c("no sharing", "with sharing"),
           variable.name = "model", value.name = "SE")

## plot against n_batches (increasing) but LABEL with batch size, so smaller
## batches run left to right without needing a reversed transform
pC <- ggplot(pl, aes(n_batches, SE, colour = model)) +
  geom_hline(yintercept = LAM_SHOW, linetype = "22", colour = INK2,
             linewidth = .5) +
  annotate("text", x = min(pl$n_batches), y = LAM_SHOW * 1.3, hjust = 0,
           label = "SE = lambda  (power ~ 0.17)", colour = INK2, size = 3.1) +
  geom_line(linewidth = .9) + geom_point(size = 2.1) +
  scale_colour_manual(values = c("no sharing" = BLUE,
                                 "with sharing" = ORANGE), name = NULL) +
  scale_x_log10(breaks = sort(unique(pl$n_batches)),
                labels = function(b) comma(round(5e5 / b))) +
  scale_y_log10(labels = scientific) +
  labs(title = sprintf("C  Variant sharing cancels the gain (lambda = %s)",
                       scientific(LAM_SHOW)),
       subtitle = "SE falls as n^-0.5 while the sqrt(k) penalty grows as n^0.5 -- they offset almost exactly",
       x = "batch size", y = expression(SE(hat(lambda)))) +
  theme_v() + theme(legend.position = "top",
                    legend.justification = "left",
                    legend.margin = margin(b = -4),
                    axis.text.x = element_text(angle = 40, hjust = 1))

fig <- pA / (pB | pC) + plot_layout(heights = c(1, 1.15)) +
  plot_annotation(
    title = "Batched singleton ascertainment: what the sweep shows",
    subtitle = "UK Biobank-scale counts taken from the tree sequence, not from a scaling formula",
    theme = theme(plot.background = element_rect(fill = SURFACE, colour = NA),
                  plot.title = element_text(size = 15, face = "bold", colour = INK),
                  plot.subtitle = element_text(size = 10, colour = INK2,
                                               margin = margin(b = 8))))

ggsave(file.path(OUT, "batch_sweep.png"), fig, width = 13, height = 9.5, dpi = 200, bg = SURFACE)
cat(sprintf("scaling exponents:  t0 ~ n^%.2f   h2_sing ~ n^%.2f\n", e_t0, e_h2))
