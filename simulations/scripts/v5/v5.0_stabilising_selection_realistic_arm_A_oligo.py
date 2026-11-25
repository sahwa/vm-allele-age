#!/usr/bin/env python3
"""Compare singleton squared effects with mutational input in SLiM tree sequences.

Usage:
    python singleton_recovery_birth_v3.py /path/to/trees \
        --birth-dir /path/to/trees --out diagnostics.tsv --boot 2000

For each arm, report pooled (mutation-count-weighted) mean(beta**2) ratios:
  sing/nom, young/nom, sing/all, sing/young (original diagnostics)
  birth/nom, sing/birth, young/birth (new birth diagnostics)

Birth comparisons use ONLY replicates with matching .birth.tsv files.
Bootstrap uncertainty resamples whole simulation replicates, not mutations.
A single replicate cannot give a replicate-bootstrap SE.

CAUTION: age <= 1 does not guarantee a mutation has avoided selection.
CAUTION: these are tree-sequence sample singletons, which may be population
         singletons if the entire simulated population was recorded.
"""
import argparse
import csv
import glob
import os
import re
import sys

import numpy as np
import tskit

# arm: (V_S, MEAN_EFFECT, SD_EFFECT, NEUTRAL)
ARMS = {
    "null":       (5.0, 0.00, 0.1, True),
    "realistic":  (5.0, 0.00, 0.1, False),
    "loud":       (1.0, 0.00, 0.1, False),
    "veryloud":   (0.2, 0.00, 0.1, False),
    "skewmild":   (5.0, 0.05, 0.1, False),
    "skewstrong": (5.0, 0.10, 0.1, False),
}

FIELDS = ("n_s", "S_s", "n_a", "S_a", "n_y", "S_y", "n_b", "S_b")
METRICS = (
    "sing_over_nom", "young_over_nom", "sing_over_all", "sing_over_young",
    "birth_over_nom", "sing_over_birth", "young_over_birth",
)


def arm_of(path):
    match = re.fullmatch(r"([a-z]+)_rep\d+\.trees", os.path.basename(path))
    return match.group(1) if match and match.group(1) in ARMS else None


def extract(path):
    """Extract effects/counts/ages from single-mutation segregating sites."""
    ts = tskit.load(path)
    n_samples = ts.num_samples
    tick = (ts.metadata.get("SLiM", {}).get("tick")
            if isinstance(ts.metadata, dict) else None)
    counts, effects, ages = [], [], []
    for variant in ts.variants():
        if len(variant.site.mutations) != 1:
            continue
        mutation = variant.site.mutations[0]
        try:
            mutation_list = mutation.metadata["mutation_list"]
            if len(mutation_list) != 1:
                continue
            beta = float(mutation_list[0]["selection_coeff"])
        except (KeyError, IndexError, TypeError, ValueError):
            continue
        count = int(np.count_nonzero(variant.genotypes))
        if count == 0 or count == n_samples:
            continue
        age = float(mutation.time)
        if not np.isfinite(age) and tick is not None:
            age = float(tick - mutation_list[0]["slim_time"])
        counts.append(count)
        effects.append(beta)
        ages.append(age)
    return (np.asarray(counts, dtype=np.int64),
            np.asarray(effects, dtype=float),
            np.asarray(ages, dtype=float))


def read_birth(path):
    if not os.path.isfile(path):
        return None
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"expected exactly one data row in {path}")
    row = rows[0]
    n = float(row["n_birth"])
    total = float(row["sum_beta2"])
    if not np.isfinite(n) or not np.isfinite(total) or n <= 0 or total < 0:
        raise ValueError(f"invalid birth data in {path}")
    return n, total


def safe_ratio(numerator, denominator):
    return numerator / denominator if denominator > 0 else float("nan")


def summarize(runs, nominal):
    """Return ratios; paired birth comparisons use birth-logged replicates only."""
    totals = {f: sum(r[f] for r in runs) for f in FIELDS}
    sing = safe_ratio(totals["S_s"], totals["n_s"])
    young = safe_ratio(totals["S_y"], totals["n_y"])
    all_sites = safe_ratio(totals["S_a"], totals["n_a"])
    out = {
        "sing_over_nom": sing / nominal,
        "young_over_nom": young / nominal,
        "sing_over_all": safe_ratio(sing, all_sites),
        "sing_over_young": safe_ratio(sing, young),
        "birth_over_nom": float("nan"),
        "sing_over_birth": float("nan"),
        "young_over_birth": float("nan"),
    }
    paired = [r for r in runs if r["n_b"] > 0]
    if paired:
        sums = {f: sum(r[f] for r in paired) for f in FIELDS}
        birth = safe_ratio(sums["S_b"], sums["n_b"])
        paired_sing = safe_ratio(sums["S_s"], sums["n_s"])
        paired_young = safe_ratio(sums["S_y"], sums["n_y"])
        out.update(birth_over_nom=birth / nominal,
                   sing_over_birth=safe_ratio(paired_sing, birth),
                   young_over_birth=safe_ratio(paired_young, birth))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("directory", nargs="?", default=".")
    ap.add_argument("--birth-dir", default=None,
                    help="Directory containing matching *_repN.birth.tsv files (default: tree directory)")
    ap.add_argument("--out", default=None, help="Write arm-level TSV summary")
    ap.add_argument("--young-age", type=float, default=1.0)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    birth_dir = args.birth_dir or args.directory

    paths = sorted(glob.glob(os.path.join(args.directory, "*_rep*.trees")))
    if not paths:
        sys.exit(f"No *_rep*.trees files in {args.directory}")

    runs = {arm: [] for arm in ARMS}
    for path in paths:
        arm = arm_of(path)
        if arm is None:
            continue
        try:
            counts, betas, ages = extract(path)
            birth_path = os.path.join(
                birth_dir, os.path.splitext(os.path.basename(path))[0] + ".birth.tsv")
            birth = read_birth(birth_path)
        except Exception as exc:
            print(f"SKIP {path}: {exc}", file=sys.stderr)
            continue
        if len(counts) == 0:
            print(f"SKIP {path}: no usable segregating sites", file=sys.stderr)
            continue
        b2 = betas * betas
        sing = counts == 1
        young = np.isfinite(ages) & (ages >= 0) & (ages <= args.young_age)
        r = dict(n_s=int(sing.sum()), S_s=float(b2[sing].sum()),
                 n_a=len(b2), S_a=float(b2.sum()),
                 n_y=int(young.sum()), S_y=float(b2[young].sum()),
                 n_b=birth[0] if birth else 0.0,
                 S_b=birth[1] if birth else 0.0)
        runs[arm].append(r)
        nominal = ARMS[arm][1] ** 2 + ARMS[arm][2] ** 2
        print(f"{os.path.basename(path):27s} sing={r['n_s']:5d} "
              f"young={r['n_y']:5d} birth={'yes' if birth else 'NO '} "
              f"sing/nom={safe_ratio(r['S_s'],r['n_s'])/nominal:.5f}",
              file=sys.stderr)

    rng = np.random.default_rng(args.seed)
    summary_rows = []
    print("\nRatios use pooled sums across replicates; birth ratios use only paired replicates.")
    print("SEs are replicate bootstrap SEs (NA when only one replicate).\n")
    for arm, arm_runs in runs.items():
        if not arm_runs:
            continue
        nominal = ARMS[arm][1] ** 2 + ARMS[arm][2] ** 2
        point = summarize(arm_runs, nominal)
        se = {m: float("nan") for m in METRICS}
        if len(arm_runs) >= 2 and args.boot >= 2:
            estimates = {m: [] for m in METRICS}
            for _ in range(args.boot):
                indices = rng.integers(0, len(arm_runs), len(arm_runs))
                b = summarize([arm_runs[i] for i in indices], nominal)
                for m in METRICS:
                    if np.isfinite(b[m]):
                        estimates[m].append(b[m])
            for m in METRICS:
                if len(estimates[m]) >= 2:
                    se[m] = float(np.std(estimates[m], ddof=1))

        n_paired = sum(r["n_b"] > 0 for r in arm_runs)
        print(f"{arm}: reps={len(arm_runs)}, paired birth logs={n_paired}, "
              f"singletons={sum(r['n_s'] for r in arm_runs):,}")
        for m in METRICS:
            val, err = point[m], se[m]
            print(f"  {m:19s} {val:9.5f} +/- "
                  f"{err:.5f}" if np.isfinite(val) else f"  {m:19s} NA")
        summary_rows.append({"arm": arm, "reps": len(arm_runs),
                             "birth_reps": n_paired,
                             "singletons": sum(r["n_s"] for r in arm_runs),
                             "young": sum(r["n_y"] for r in arm_runs),
                             "birth_mutations": int(sum(r["n_b"] for r in arm_runs)),
                             "lambda_pred": (0.0 if ARMS[arm][3]
                                             else nominal / ARMS[arm][0]),
                             **point, **{m + "_se": se[m] for m in METRICS}})
    if args.out and summary_rows:
        with open(args.out, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(summary_rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(summary_rows)
        print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
