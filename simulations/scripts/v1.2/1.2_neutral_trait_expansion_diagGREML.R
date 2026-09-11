source("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/diagGREML.R")

C <- fread("1.2_singleton_counts.csv")
pheno = fread("1.2_phenotypes.csv")

bins <- names(C)

A <- c(
  lapply(bins[1:3], function(b) C[[b]]),
  list(old = rowSums(C[, bins[4:8], with = FALSE]))
)

A <- lapply(A, function(x) x / mean(x))

names(A) <- c(bins[1:3], "old")
X <- matrix(1, nrow = nrow(C), ncol = 1)

fit <- fit_diagGREML(y = pheno$y, A = A, X = X,
                     constraint = FALSE, magic0316 = TRUE)
data.table(comp = fit$Vlistnames, est = fit$varcmp,
           se = sqrt(diag(fit$Hi)))