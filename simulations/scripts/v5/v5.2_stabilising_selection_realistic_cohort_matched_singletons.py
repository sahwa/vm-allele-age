
#!/usr/bin/env python3
"""
Compare neutral singletons with all mutations born in the
same generations.

Requires v5.2:
  null_rep*.trees
  null_rep*.birth_recent.tsv

Cohorts (generations before END_TICK):
  0
  1-4
  5-19
  20-99

Birth denominators are obtained by subtracting the
cumulative birth windows recorded by SLiM.

Run:
python3 cohort_matched_singletons.py /path/to/v5.2/armA_oligo
"""

import argparse
import csv
import glob
import os
import re

import numpy as np
import tskit

ap = argparse.ArgumentParser()
ap.add_argument("directory")
ap.add_argument("--end-tick", type=int, default=75000)
ap.add_argument("--boot", type=int, default=2000)
ap.add_argument("--seed", type=int, default=42)
ap.add_argument("--out", default="cohort_matched.tsv")
args = ap.parse_args()

# Non-overlapping cohorts.
COHORTS = [
    ("0", 0, 1),
    ("1-4", 1, 5),
    ("5-19", 5, 20),
    ("20-99", 20, 100),
]

# Per-replicate rows:
# birth_n, birth_sum_beta2, sing_n, sing_sum_beta2
data = {name: [] for name, _, _ in COHORTS}


def get_births(path):
    """Read cumulative birth windows and subtract them."""
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    cumulative = {
        int(float(r["window_generations"])): (
            float(r["n_birth"]),
            float(r["sum_beta2"])
        )
        for r in rows
    }

    result = {}
    previous_n = 0.0
    previous_s2 = 0.0

    for name, _, upper in COHORTS:
        n, s2 = cumulative[upper]
        result[name] = (
            n - previous_n,
            s2 - previous_s2
        )
        previous_n, previous_s2 = n, s2

    return result


def get_singletons(path):
    """Extract singleton effects by SLiM mutation birth tick."""
    ts = tskit.load(path)

    result = {
        name: [0, 0.0]
        for name, _, _ in COHORTS
    }

    for v in ts.variants():

        if len(v.site.mutations) != 1:
            continue

        mut = v.site.mutations[0]

        try:
            ml = mut.metadata["mutation_list"]
            if len(ml) != 1:
                continue

            beta = float(ml[0]["selection_coeff"])
            birth_tick = int(ml[0]["slim_time"])

        except (KeyError, TypeError, ValueError, IndexError):
            continue

        if not np.isfinite(beta):
            continue

        # Same singleton definition as original extraction.
        if np.count_nonzero(v.genotypes) != 1:
            continue

        age = args.end_tick - birth_tick

        for name, lower, upper in COHORTS:
            if lower <= age < upper:
                result[name][0] += 1
                result[name][1] += beta**2
                break

    return result


# Read all matching replicates.
paths = sorted(
    p for p in glob.glob(
        os.path.join(args.directory, "null_rep*.trees")
    )
    if re.fullmatch(
        r"null_rep\d+\.trees",
        os.path.basename(p)
    )
)

if not paths:
    ap.error("No null replicate tree files found")

for path in paths:

    prefix = path[:-6]
    birth_path = prefix + ".birth_recent.tsv"

    if not os.path.exists(birth_path):
        print("Missing birth file:", birth_path)
        continue

    births = get_births(birth_path)
    singletons = get_singletons(path)

    for name, _, _ in COHORTS:
        bn, bs2 = births[name]
        sn, ss2 = singletons[name]

        data[name].append([bn, bs2, sn, ss2])

    print(
        os.path.basename(path),
        "matched successfully"
    )


# Bootstrap whole replicates, retaining pairing.
rng = np.random.default_rng(args.seed)
output = []

print("\nCOHORT-MATCHED SINGLETON RECOVERY")
print("Expected ratio under neutrality = 1\n")

for name, _, _ in COHORTS:

    if not data[name]:
        continue

    R = np.asarray(data[name], dtype=float)

    bn, bs2, sn, ss2 = R.sum(axis=0)

    if bn <= 0 or sn <= 0:
        print(name, "insufficient mutations")
        continue

    mean_birth = bs2 / bn
    mean_sing = ss2 / sn

    ratio = mean_sing / mean_birth

    # Replicate bootstrap.
    estimates = []

    if len(R) > 1:
        for _ in range(args.boot):

            idx = rng.integers(
                0, len(R), len(R)
            )

            b_n, b_s2, s_n, s_s2 = R[idx].sum(axis=0)

            if b_n > 0 and s_n > 0:
                estimates.append(
                    (s_s2 / s_n) / (b_s2 / b_n)
                )

    se = (
        np.std(estimates, ddof=1)
        if len(estimates) > 1 else np.nan
    )

    print(
        f"Age {name:6s} | "
        f"births={int(bn):9,d} | "
        f"singletons={int(sn):7,d} | "
        f"birth/nom={mean_birth/0.01:.5f} | "
        f"sing/nom={mean_sing/0.01:.5f} | "
        f"sing/birth={ratio:.5f} +/- {se:.5f}"
    )

    output.append({
        "cohort": name,
        "reps": len(R),
        "birth_n": int(bn),
        "singleton_n": int(sn),
        "birth_over_nom": mean_birth / 0.01,
        "singleton_over_nom": mean_sing / 0.01,
        "singleton_over_birth": ratio,
        "se": se
    })


if output:
    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(output[0].keys()),
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(output)

    print("\nWrote", args.out)

print("\nNote: cohorts use SLiM origin ticks,")
print("not tskit mutation.time.")



# Pooled cohort-adjusted ratio with replicate bootstrap.
# Each row is one replicate; columns are:
# birth_n, birth_sum_beta2, singleton_n, singleton_sum_beta2

R = np.array([data[name] for name, _, _ in COHORTS])
assert len({len(data[name]) for name, _, _ in COHORTS}) == 1

def adjusted_ratio(x):
    totals = x.sum(axis=1)
    bn, bs2, sn, ss2 = totals.T
    return ss2.sum() / np.sum(sn * bs2 / bn)

observed = adjusted_ratio(R)

rng = np.random.default_rng(args.seed)
boots = []
for _ in range(args.boot):
    idx = rng.integers(0, R.shape[1], R.shape[1])
    boots.append(adjusted_ratio(R[:, idx, :]))

se = np.std(boots, ddof=1)
ci = np.percentile(boots, [2.5, 97.5])

print("\nPOOLED COHORT-ADJUSTED RECOVERY")
print(f"Ratio:   {observed:.6f}")
print(f"SE:      {se:.6f}")
print(f"95% CI:  [{ci[0]:.6f}, {ci[1]:.6f}]")
print(f"z vs 1:  {(observed - 1) / se:.2f}")
