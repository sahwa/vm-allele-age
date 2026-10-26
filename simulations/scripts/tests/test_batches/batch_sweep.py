#!/usr/bin/env python3
"""
batch_sweep.py -- measure the three ceilings of the batched-singleton design
directly from the tree sequence instead of assuming them.

ONE streaming pass over the trees fills, for every batch level n_b at once:

  C[n_b]      (n_ind x K) int32   batch-singletons per person per age bin
  kdist[n_b]  histogram           how many batches each variant is a singleton in
  diagnostics                     t0, T per person, cv_total, V_struct, V_samp

What each output settles:
  ceiling 2 (saturation)  -- t0(n_b) and T(n_b) flatten on their own; no T0_CAP
  ceiling 1 (sharing)     -- kdist is the real k, replacing the 0.37*n_b guess
  ceiling 3 (validity)    -- feed C into the estimator and check sigma2_m is
                             flat across n_b; where it drifts, batch-singletons
                             have stopped being recent mutations

Age bins are ONE fixed log2 grid shared by every level, so the C matrices and
every derived quantity are directly comparable across n_b.

Notes on cost: the bottleneck is the per-mutation Python loop, so start at the
reduced-L tree sequence. Everything here is either per-bp invariant (t0,
CV_age, kdist) or linear in L (T per person, hence h2_sing), so the shape of
all three ceilings is visible at L=1e7 and only the absolute counts need the
full run.
"""


import argparse
import json
import pathlib

import numpy as np
import tskit

# ---------------------------------------------------------------------------
# config
# ---------------------------------------------------------------------------
MU = 1.25e-8                       # per-bp per-generation mutation rate
N_BATCHES = [1, 2, 5, 10, 20, 50, 100, 200]   # NUMBER of batches, not size
MAC_SLACK = 5                      # keep MAC <= MAC_SLACK * max(n_b); see below

# One shared log2 age grid. Fixed across levels on purpose.
AGE_EDGES = np.array([0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, np.inf])


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--trees", required=True, type=pathlib.Path)
    p.add_argument("--out", required=True, type=pathlib.Path)
    p.add_argument("--n-batches", type=int, nargs="+", default=N_BATCHES)
    p.add_argument("--mac-slack", type=int, default=MAC_SLACK)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--progress-every", type=int, default=200_000)
    return p.parse_args()


# ---------------------------------------------------------------------------
# sample -> individual map
# ---------------------------------------------------------------------------
def build_node_to_ind(ts):
    """Map each sample NODE to an individual index in [0, n_ind).

    Uses the individual table when msprime wrote one; otherwise pairs
    consecutive sample nodes as the two haplotypes of one person.
    """
    samples = ts.samples()
    node_ind = ts.tables.nodes.individual[samples]

    if np.all(node_ind >= 0) and ts.num_individuals > 0:
        uniq = np.unique(node_ind)
        remap = np.full(uniq.max() + 1, -1, dtype=np.int32)
        remap[uniq] = np.arange(uniq.size, dtype=np.int32)
        ind_of_sample = remap[node_ind]
        n_ind = uniq.size
    else:
        # haploid sample nodes, paired in order
        assert samples.size % 2 == 0, "odd number of sample nodes; set ploidy explicitly"
        ind_of_sample = np.repeat(np.arange(samples.size // 2, dtype=np.int32), 2)
        n_ind = samples.size // 2

    # index by NODE id, not by position in samples()
    node_to_ind = np.full(ts.num_nodes, -1, dtype=np.int32)
    node_to_ind[samples] = ind_of_sample
    return node_to_ind, n_ind


# ---------------------------------------------------------------------------
# nested batch assignment
# ---------------------------------------------------------------------------
def build_batch_maps(n_ind, levels, seed):
    """batch_of[level][i] = which batch person i is in, at that level.

    One random permutation drives every level: person of permutation-rank r
    goes in batch r*n_b//n_ind. Levels that divide each other are then exactly
    nested, so varying n_b varies granularity rather than reshuffling people.
    """
    rank = np.empty(n_ind, dtype=np.int64)
    rank[np.random.default_rng(seed).permutation(n_ind)] = np.arange(n_ind)
    return {nb: (rank * nb // n_ind).astype(np.int32) for nb in levels}


# ---------------------------------------------------------------------------
# age binning
# ---------------------------------------------------------------------------
def bin_ages(ages, edges=AGE_EDGES):
    """Index into edges for each age; ages are clamped to >= 0."""
    return np.clip(np.searchsorted(edges, ages, side="right") - 1,
                   0, len(edges) - 2).astype(np.int32)


# ---------------------------------------------------------------------------
# the pass
# ---------------------------------------------------------------------------
def sweep(ts, levels, mac_slack, seed, progress_every):
    node_to_ind, n_ind = build_node_to_ind(ts)
    batch_of = build_batch_maps(n_ind, levels, seed)
    K = len(AGE_EDGES) - 1
    q_max = mac_slack * max(levels)

    tabs = ts.tables
    mut_node = tabs.mutations.node
    mut_time = tabs.mutations.time
    mut_site = tabs.mutations.site
    # mutation time is nan when msprime did not record it; fall back to the
    # node's time, which is the age of the branch the mutation sits on
    ages = np.where(np.isnan(mut_time), tabs.nodes.time[mut_node], mut_time)
    mut_bin = bin_ages(ages)
    site_pos = tabs.sites.position[mut_site]

    # accumulators, one per level
    C = {nb: np.zeros((n_ind, K), dtype=np.int32) for nb in levels}
    age_sum = {nb: np.zeros(K) for nb in levels}       # for count-weighted t_rep
    age_n = {nb: np.zeros(K, dtype=np.int64) for nb in levels}
    kdist = {nb: np.zeros(nb + 1, dtype=np.int64) for nb in levels}

    n_seen = n_kept = 0
    for tree in ts.trees():
        lo, hi = np.searchsorted(site_pos, tree.interval)
        for m in range(lo, hi):
            n_seen += 1
            if progress_every and n_seen % progress_every == 0:
                print(f"  {n_seen:,} mutations scanned, {n_kept:,} kept", flush=True)

            u = mut_node[m]
            q = tree.num_samples(u)          # O(1), maintained incrementally
            if q == 0 or q > q_max:
                continue
            n_kept += 1

            carrier_nodes = np.fromiter(tree.samples(u), dtype=np.int64, count=q)
            carrier_ind = node_to_ind[carrier_nodes]
            b = mut_bin[m]
            a = ages[m]

            for nb in levels:
                if nb == 1:
                    # one batch = the whole cohort; singleton iff MAC == 1
                    if q != 1:
                        kdist[1][0] += 1
                        continue
                    hit = carrier_ind
                    k = 1
                else:
                    bid = batch_of[nb][carrier_ind]
                    vals, cnt = np.unique(bid, return_counts=True)
                    solo = vals[cnt == 1]
                    k = solo.size
                    if k == 0:
                        kdist[nb][0] += 1
                        continue
                    hit = carrier_ind[np.isin(bid, solo)]

                kdist[nb][k] += 1
                np.add.at(C[nb], (hit, b), 1)
                age_sum[nb][b] += a * hit.size
                age_n[nb][b] += hit.size

    print(f"scanned {n_seen:,} mutations, kept {n_kept:,} (MAC <= {q_max})")
    return C, age_sum, age_n, kdist, n_ind


# ---------------------------------------------------------------------------
# diagnostics -- the python mirror of measure_knobs() in sim_common.R
# ---------------------------------------------------------------------------
def diagnose(Cm, t_rep):
    ok_bin = np.isfinite(t_rep) & (t_rep > 0)
    C = Cm[:, ok_bin].astype(np.float64)
    lt = np.log(t_rep[ok_bin])

    Tt = C.sum(1)
    ok = Tt >= 2
    if ok.sum() < 100:
        return dict(mean_total=float(Tt.mean()), cv_total=np.nan, t0=np.nan,
                    var_logage_total=np.nan, var_logage_sampling=np.nan,
                    cv_age=np.nan, cv_age_eff=np.nan, frac_usable=float(ok.mean()))

    Lbar = (C[ok] @ lt) / Tt[ok]
    s2 = (C[ok] @ lt**2) / Tt[ok] - Lbar**2
    v_tot = float(np.var(Lbar, ddof=1))
    v_samp = float(np.mean(s2 / (Tt[ok] - 1)))
    v_struct = max(v_tot - v_samp, 0.0)

    colsum = C.sum(0)
    return dict(
        mean_total=float(Tt.mean()),
        cv_total=float(Tt.std(ddof=1) / Tt.mean()),
        t0=float((colsum * t_rep[ok_bin]).sum() / colsum.sum()),
        var_logage_total=v_tot,
        var_logage_sampling=v_samp,
        cv_age=float(np.sqrt(v_struct)),
        # effective CV_age after regression dilution -- what the power
        # formula's denominator actually wants
        cv_age_eff=float(v_struct / np.sqrt(v_tot)) if v_tot > 0 else np.nan,
        frac_usable=float(ok.mean()),
    )


def main():
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    levels = sorted(set(args.n_batches))

    ts = tskit.load(args.trees)
    print(f"{args.trees.name}: {ts.num_samples:,} samples, "
          f"{ts.num_mutations:,} mutations, {ts.sequence_length:,.0f} bp")

    C, age_sum, age_n, kdist, n_ind = sweep(
        ts, levels, args.mac_slack, args.seed, args.progress_every)

    rows = []
    for nb in levels:
        t_rep = np.where(age_n[nb] > 0, age_sum[nb] / np.maximum(age_n[nb], 1), np.nan)
        d = diagnose(C[nb], t_rep)

        kd = kdist[nb]
        shared = kd[1:]                       # variants that were a singleton somewhere
        kk = np.arange(1, nb + 1)
        n_var = max(shared.sum(), 1)
        ek = (kk * shared).sum() / n_var              # E[k | k >= 1], per VARIANT
        ek2 = (kk**2 * shared).sum() / n_var
        d.update(
            n_batches=nb,
            batch_size=n_ind // nb,
            k_per_variant=float(ek),
            # instance-weighted redundancy E[k^2]/E[k]: the number of OTHER rows
            # a typical batch-singleton row shares its variant with. This is the
            # one the SE penalty needs, and the one 0.37*n_b was standing in for.
            k_per_instance=float(ek2 / ek) if ek > 0 else np.nan,
            k_pred_0p37=0.37 * nb,            # the old assumed value
            k_pred_sfs=0.25 * nb,             # SFS-weighted prediction for E[k^2]/E[k]
            n_variants_used=int(shared.sum()),
            n_instances=int((kk * shared).sum()),
        )
        rows.append(d)

        np.save(args.out / f"C_nb{nb}.npy", C[nb])
        np.save(args.out / f"trep_nb{nb}.npy", t_rep)
        np.save(args.out / f"kdist_nb{nb}.npy", kd)
        print(f"nb={nb:>4}  B={d['batch_size']:>7,}  T={d['mean_total']:>8.1f}  "
              f"t0={d['t0']:>7.1f}  cv_tot={d['cv_total']:.3f}  "
              f"cv_age={d['cv_age']:.3f}  cv_eff={d['cv_age_eff']:.3f}  "
              f"k/var={d['k_per_variant']:.2f}  k/inst={d['k_per_instance']:.2f}  "
              f"(0.37nb={0.37*nb:.2f})")

    np.save(args.out / "age_edges.npy", AGE_EDGES)
    with open(args.out / "diagnostics.json", "w") as fh:
        json.dump(rows, fh, indent=2)
    print(f"\nwritten to {args.out}")


if __name__ == "__main__":
    main()
