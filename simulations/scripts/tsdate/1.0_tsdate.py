import tskit
import tsdate
import pathlib
import daiquiri
import numpy as np
import pandas as pd


daiquiri.setup(level="INFO")

N=500000
L="1e+08"

MU = 1.25e-8                  # per-bp per-generation mutation rate
REC = 1e-8                    # per-bp per-generation recombination rate


DATA = f"/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.2/n{N}_L{L}"
DATES_FILE = pathlib.Path(DATA, "sim.trees.dates")

ts = tskit.load(pathlib.Path(DATA, "sim.trees"))

if DATES_FILE.exists():
    dates = tskit.load(DATES_FILE)
else:
    dates = tsdate.date(ts, mutation_rate = MU)
    dates.dump(DATES_FILE)
    mut_node = tabs.mutations.node


def mutation_ages(ts):
    """Age of each mutation: its own time when recorded, else its node's time."""
    tabs = ts.tables
    mt = np.asarray(tabs.mutations.time, dtype=np.float64)
    node_t = np.asarray(tabs.nodes.time, dtype=np.float64)
    return np.where(np.isnan(mt), node_t[tabs.mutations.node], mt)


def align_mutations(ts, dated, verbose=True):
    """Row index in `dated` for each mutation row in `ts`.

    tsdate preserves topology, sites and mutations, but tskit requires the
    mutations at a site to be ordered parent-before-child. Assigning new times
    changes that order, so TableCollection.sort() transposes mutations within
    multi-mutation sites. Row i of the two mutation tables is therefore NOT the
    same mutation; (site, node) is the stable key.

    Returns `to_dated`, so that:
        ages_dated[to_dated[m]]  is the inferred age of ts mutation m
    """
    a, b = ts.tables.mutations, dated.tables.mutations
    if len(a) != len(b):
        raise ValueError(f"mutation tables differ in length: {len(a)} vs {len(b)}")

    ka = np.lexsort((a.node, a.site))
    kb = np.lexsort((b.node, b.site))
    if not (np.array_equal(a.site[ka], b.site[kb])
            and np.array_equal(a.node[ka], b.node[kb])):
        raise ValueError("(site, node) multisets differ -- 'dated' is not "
                         "tsdate output for this tree sequence")

    to_dated = np.empty(len(a), dtype=np.int64)
    to_dated[ka] = kb

    if verbose:
        moved = int((to_dated != np.arange(len(a))).sum())
        key = a.site[ka].astype(np.int64) * (ts.num_nodes + 1) + a.node[ka]
        dup = int((np.diff(key) == 0).sum())
        print(f"aligned {len(a):,} mutations; {moved:,} ({moved/len(a):.2%}) "
              f"had moved row")
        if dup:
            print(f"  note: {dup:,} share a (site, node) with another mutation "
                  f"and are matched arbitrarily within that pair -- same branch, "
                  f"so the age difference is bounded by the branch length")
    return to_dated


def derived_counts(ts):
    """Derived allele count for every mutation, as an array indexed by row.

    The only part that has to be a Python loop: tree.num_samples(u) is O(1)
    but there is no batch API for it. Writing into a preallocated int32 array
    keeps the loop body to a single store, so the walk costs ~1.5us per
    mutation and everything downstream stays vectorised. Accumulating Python
    lists here instead costs ~7GB and most of the runtime at 40M mutations.
    """
    mut_node = ts.tables.mutations.node
    site_pos = ts.tables.sites.position[ts.tables.mutations.site]
    q = np.zeros(ts.num_mutations, dtype=np.int32)
    for tree in ts.trees():
        lo, hi = np.searchsorted(site_pos, tree.interval)
        for m in range(lo, hi):
            q[m] = tree.num_samples(mut_node[m])
    return q


def variant_age_table(ts, dated, mac_max=0, mac_keep_all=0, thin=0.005, seed=1):
    """One row per sampled variant: MAC, true age, inferred age, weight.

    MAC comes from the TRUE tree sequence; tsdate changes node times, not
    topology or genotypes, so allele counts are identical in both. The mutation
    ROWS are not aligned though -- see align_mutations().

    thin=1 keeps everything (the default 0.005 is for quick looks). Set it as
    target_rows / ts.num_mutations.
    """
    t_true = mutation_ages(ts)
    t_inf = mutation_ages(dated)[align_mutations(ts, dated)]

    n = ts.num_samples
    if mac_max <= 0:
        mac_max = n // 2
    site_pos = ts.tables.sites.position[ts.tables.mutations.site]

    q = derived_counts(ts)
    mac = np.minimum(q, n - q).astype(np.int32)

    keep = (q > 0) & (q < n) & (mac <= mac_max)
    w = np.ones(ts.num_mutations)
    if thin < 1.0:
        rng = np.random.default_rng(seed)
        thinned = mac > mac_keep_all
        keep &= ~thinned | (rng.random(ts.num_mutations) <= thin)
        w[thinned] = 1.0 / thin

    df = pd.DataFrame(dict(
        position=site_pos[keep], mac=mac[keep],
        age_true=t_true[keep], age_inferred=t_inf[keep], weight=w[keep]))
    df["maf"] = df["mac"] / n
    print(f"scanned {ts.num_mutations:,} mutations, kept {len(df):,}")
    return df

df = variant_age_table(ts, dates, thin=1)

df.to_csv("/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/tsdate/variant_age_table.csv")
