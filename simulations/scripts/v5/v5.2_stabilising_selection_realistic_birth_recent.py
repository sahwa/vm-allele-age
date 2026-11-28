
#!/usr/bin/env python3
"""Summarise recent mutation-at-birth diagnostics across SLiM replicates."""

import argparse
import csv
import glob
import os
import re

import numpy as np

ARMS = {
    "null":       (0.00, 0.1),
    "realistic":  (0.00, 0.1),
    "loud":       (0.00, 0.1),
    "veryloud":   (0.00, 0.1),
    "skewmild":   (0.05, 0.1),
    "skewstrong": (0.10, 0.1),
}

WINDOWS = (1, 5, 20, 100)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("directory")
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default="recent_birth_summary.tsv")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    data = {arm: {w: [] for w in WINDOWS} for arm in ARMS}

    files = glob.glob(
        os.path.join(args.directory, "*_rep*.birth_recent.tsv")
    )

    if not files:
        ap.error("No .birth_recent.tsv files found")

    for path in files:
        name = os.path.basename(path)
        match = re.fullmatch(
            r"([a-z]+)_rep\d+\.birth_recent\.tsv", name
        )

        if not match or match.group(1) not in ARMS:
            continue

        arm = match.group(1)

        with open(path, newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))

        for row in rows:
            w = int(float(row["window_generations"]))

            if w not in WINDOWS:
                continue

            n = float(row["n_birth"])
            s2 = float(row["sum_beta2"])

            if n > 0:
                data[arm][w].append((n, s2))

    output = []

    print("\nRECENT MUTATION-AT-BIRTH DIAGNOSTICS")
    print("Ratios relative to nominal E[beta^2]\n")

    for arm, windows in data.items():
        mean, sd = ARMS[arm]
        nominal = mean**2 + sd**2

        for w, values in windows.items():
            if not values:
                continue

            R = np.asarray(values, dtype=float)

            n_total = R[:, 0].sum()
            s2_total = R[:, 1].sum()

            ratio = (s2_total / n_total) / nominal

            if len(R) > 1 and args.boot > 1:
                estimates = []

                for _ in range(args.boot):
                    idx = rng.integers(
                        0, len(R), size=len(R)
                    )
                    sample = R[idx]

                    est = (
                        sample[:, 1].sum() /
                        sample[:, 0].sum() /
                        nominal
                    )
                    estimates.append(est)

                se = np.std(estimates, ddof=1)
            else:
                se = np.nan

            print(
                f"{arm:12s} last {w:3d} gen  "
                f"reps={len(R):3d}  "
                f"n={int(n_total):10,d}  "
                f"ratio={ratio:.6f}  "
                f"SE={se:.6f}"
            )

            output.append({
                "arm": arm,
                "window": w,
                "reps": len(R),
                "n_birth": int(n_total),
                "ratio_to_nominal": ratio,
                "se": se,
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

        print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
