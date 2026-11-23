import tskit, numpy as np

ts = tskit.load("/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/daiki_sim_data/tszip/chr20p_yri_10000_ceu_500000_chb_10000_jpt_10000_neutralmu_1.7e-08_neutralseed_12.tree")

# Optional: subsample to UKB-like f (50k / 7M ≈ 0.007)
rng = np.random.default_rng(1)
samp = ts.samples()
ind_of = ts.tables.nodes.individual[samp]
keep = rng.choice(np.unique(ind_of), 50_000, replace=False)
ts = ts.simplify(samples=samp[np.isin(ind_of, keep)])

nmut = np.bincount(ts.tables.mutations.site, minlength=ts.num_sites)
k_, s_, t_ = [], [], []
for tree in ts.trees():
    for m in tree.mutations():
        if nmut[m.site] != 1:
            continue                                   # skip stacked sites
        ml = m.metadata["mutation_list"]
        if len(ml) != 1:
            continue
        k = tree.num_samples(m.node)                   # derived allele count
        if k == 0:
            continue
        k_.append(k); s_.append(abs(ml[0]["selection_coeff"])); t_.append(m.time)
k_, s_, t_ = map(np.array, (k_, s_, t_))

den = s_[t_ <= 2].mean()                               # youngest = unfiltered DFE
print("E[s | new] =", den, " from", (t_ <= 2).sum(), "mutations")
for k in range(1, 6):
    sel = k_ == k
    print(f"MAC {k}: kappa = {s_[sel].mean()/den:.4f}   t0 = {t_[sel].mean():.1f}   n = {sel.sum()}")
