
#!/usr/bin/env python3
"""Audit SLiM singleton extraction and permutation nulls; no plotting dependency."""
import argparse
import glob
import os
import re
from collections import Counter, defaultdict

import numpy as np
import tskit

ap = argparse.ArgumentParser()
ap.add_argument('directory')
ap.add_argument('--perms', type=int, default=2000)
ap.add_argument('--seed', type=int, default=42)
ap.add_argument('--nominal', type=float, default=0.01)
a = ap.parse_args()

paths = [p for p in sorted(glob.glob(os.path.join(a.directory, 'null_rep*.trees')))
         if re.search(r'null_rep\d+\.trees$', os.path.basename(p))]
if not paths:
    ap.error('No null_rep*.trees files found')

# Counts and sums of beta squared for each category.
stats = defaultdict(lambda: [0, 0.0])
audit = Counter()
runs = []


def record(label, values):
    stats[label][0] += len(values)
    stats[label][1] += float(np.sum(values))


def age_bin(t):
    if not np.isfinite(t):
        return 'unknown'
    if t <= 1:
        return '0-1'
    if t <= 5:
        return '1-5'
    if t <= 20:
        return '5-20'
    if t <= 100:
        return '20-100'
    return '>100'


for path in paths:
    ts = tskit.load(path)
    b2, counts, ages = [], [], []
    local = Counter()
    for v in ts.variants():
        local['sites'] += 1
        k = int(np.count_nonzero(v.genotypes))
        if k == 0 or k == ts.num_samples:
            local['absent_or_fixed'] += 1
            continue
        label = 'singleton' if k == 1 else 'other'
        local['candidate_' + label] += 1
        if len(v.site.mutations) != 1:
            local['multiple_site_mutations_' + label] += 1
            continue
        mut = v.site.mutations[0]
        try:
            ml = mut.metadata['mutation_list']
            if len(ml) != 1:
                local['stacked_metadata_' + label] += 1
                continue
            beta = float(ml[0]['selection_coeff'])
            if not np.isfinite(beta):
                raise ValueError('nonfinite beta')
        except (KeyError, TypeError, ValueError, IndexError):
            local['invalid_metadata_' + label] += 1
            continue
        local['retained_' + label] += 1
        b2.append(beta * beta)
        counts.append(k)
        ages.append(float(mut.time))

    audit.update(local)
    b2, counts, ages = map(np.asarray, (b2, counts, ages))
    if not len(b2):
        continue
    sing = counts == 1
    groups = np.array([age_bin(t) for t in ages])
    runs.append((b2, sing, groups))
    record('all', b2)
    record('singleton', b2[sing])
    record('non-singleton', b2[~sing])
    for name, mask in [('AC=1', sing), ('AC=2', counts == 2),
                       ('AC=3-5', (counts >= 3) & (counts <= 5)),
                       ('AC=6-20', (counts >= 6) & (counts <= 20)),
                       ('AC=21-100', (counts >= 21) & (counts <= 100)),
                       ('AC>100', counts > 100)]:
        record(name, b2[mask])
    for name in ['0-1', '1-5', '5-20', '20-100', '>100', 'unknown']:
        mask = groups == name
        record('age=' + name, b2[mask])
        record('age=' + name + ', singleton', b2[mask & sing])
        record('age=' + name + ', other', b2[mask & ~sing])
    print(f'{os.path.basename(path)}: {sing.sum():,} singletons; '
          f'sing/nom={b2[sing].mean()/a.nominal:.4f}')

print('\nEXTRACTION AUDIT (counts across all replicates)')
for key in sorted(audit):
    print(f'  {key:38s} {audit[key]:,}')
for label in ('singleton', 'other'):
    total = audit['candidate_' + label]
    kept = audit['retained_' + label]
    print(f'  {label} retention: {kept:,}/{total:,} = {kept/total:.3%}' if total else '')

print('\nEFFECT-SIZE DIAGNOSTICS (mean beta^2 / nominal)')
print(f'{"group":32s} {"n":>10s} {"ratio":>10s}')
for name, (n, total) in stats.items():
    if n:
        print(f'{name:32s} {n:10,d} {total/n/a.nominal:10.5f}')

# Both nulls shuffle effects only within the same simulation replicate.
# Age-adjusted additionally keeps each age category's singleton count fixed.
all_n, all_sum = stats['all']
sing_n, sing_sum = stats['singleton']
observed = (sing_sum / sing_n) / (all_sum / all_n)
rng = np.random.default_rng(a.seed)
results = np.empty((a.perms, 2))
for p in range(a.perms):
    ordinary = adjusted = 0.0
    for b2, sing, groups in runs:
        ordinary += b2[rng.choice(len(b2), size=int(sing.sum()), replace=False)].sum()
        for group in np.unique(groups):
            ids = np.flatnonzero(groups == group)
            ns = int(np.count_nonzero(sing[ids]))
            if ns:
                adjusted += b2[rng.choice(ids, size=ns, replace=False)].sum()
    results[p] = [ordinary, adjusted]
results /= sing_n * (all_sum / all_n)

print(f'\nObserved singleton/all = {observed:.6f}')
for j, label in enumerate(('within replicate', 'within replicate + age bin')):
    null = results[:, j]
    pval = (1 + np.count_nonzero(null >= observed)) / (len(null) + 1)
    print(f'{label}: null mean={null.mean():.6f}, '
          f'SD={null.std(ddof=1):.6f}, one-sided p={pval:.6g}')
