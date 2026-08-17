# Mutational variance (V_M) from age-stratified GRMs

Estimating V_M in humans by partitioning heritability by allele age.

## Pipeline
1. `slim/` - SLiM scripts (neutral trait, constant N)
2. `python/` - tree sequence processing: recapitate, overlay mutations, extract ages
3. `R/` - GRM construction, REML, kernel fitting

## Status
- [x] Neutral simulation validated: V_A empirical/analytic = 1.004
- [x] Per-bin K̄ matches analytic kernel
- [ ] GRM + REML recovery of V_bin
- [ ] Selection models
