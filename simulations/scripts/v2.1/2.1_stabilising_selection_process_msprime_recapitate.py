#!/usr/bin/env python
"""Recapitate a SLiM tree sequence and strip multiallelic sites."""

import argparse

import numpy as np
import pyslim
import tskit

p = argparse.ArgumentParser()
p.add_argument("--tree", required=True, help="SLiM .trees file")
p.add_argument("--ancestral-ne", type=float, default=25000,
               help="Ne for the recapitated (pre-SLiM) history")
p.add_argument("--rec-rate", type=float, default=1e-8)
p.add_argument("--seed", type=int, default=1)
p.add_argument("--out", default=None)
args = p.parse_args()

OUT = args.out or args.tree.replace(".trees", ".biallelic.recapitated.trees")

# --- load ---
ts = tskit.load(args.tree)
n_multi = sum(1 for t in ts.trees() if t.num_roots > 1)
print(f"Loaded: {ts.num_samples} samples, {ts.num_trees} trees, "
      f"{ts.num_sites} sites, {n_multi} ({100*n_multi/ts.num_trees:.1f}%) "
      f"uncoalesced", flush=True)

# --- recapitate ---
rts = pyslim.recapitate(ts,
                        ancestral_Ne=args.ancestral_ne,
                        recombination_rate=args.rec_rate,
                        random_seed=args.seed)
assert max(t.num_roots for t in rts.trees()) == 1, "still uncoalesced"
print("Recapitated; all trees coalesced", flush=True)

# --- Ne from branch-mode diversity (exact, no sampling) ---
# mean pairwise branch length = 2*E[TMRCA] = 4*Ne
div = rts.diversity(mode="branch")
print(f"Mean pairwise TMRCA: {div/2:.0f}   Ne estimate: {div/4:.0f} "
      f"(census {args.ancestral_ne:.0f})", flush=True)

# --- strip multiallelic sites ---
multi = np.array([s.id for s in rts.sites() if len(s.mutations) > 1],
                 dtype=np.int32)
print(f"Removing {len(multi)} multiallelic of {rts.num_sites} sites", flush=True)

if len(multi):
    tables = rts.dump_tables()
    tables.delete_sites(multi)
    tables.sort()
    tables.build_index()
    mts = tables.tree_sequence()
else:
    mts = rts

assert mts.num_mutations == mts.num_sites, "sites still carry >1 mutation"
print(f"{mts.num_sites} biallelic sites retained", flush=True)

mts.dump(OUT)
print(f"Wrote {OUT}", flush=True)
