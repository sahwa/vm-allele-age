# Mutational variance ($V_{M}$) from age-stratified GRMs

Estimating $V_{M}$ in humans by partitioning heritability by allele age.

## Pipeline
1. `slim/` - SLiM scripts (neutral trait, constant N)
2. `python/` - tree sequence processing: recapitate, overlay mutations, extract ages
3. `R/` - GRM construction, REML, kernel fitting

## Status
- [x] Neutral simulation validated: $V_{A}$ empirical/analytic = 1.004
- [x] Per-bin K̄ matches analytic kernel
- [x] GRM + REML recovery of $V_{bin}$
- [x] GENIE estimate of variance components
- [ ] Selection models
