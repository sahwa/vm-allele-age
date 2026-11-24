#!/usr/bin/env python3
"""Recovery of sigma2_m from singletons, pooled per simulation arm,
against THREE denominators.

    python3 singleton_recovery_v2.py [DIR] [--out results.tsv]
                                     [--young-age 1] [--boot 2000]

Numerator, every run:  E[beta^2 | singleton]

Denominators:
  nominal  sigma2_m typed into the config (MEAN_EFFECT^2 + SD_EFFECT^2).
           What v1 used.
  young    E[beta^2] over mutations no older than --young-age generations.
           Selection has not acted on them yet, so this is an unbiased
           estimate of the DFE as it was actually drawn in that run.
           ** sing / young is the quantity you want: kappa. **
  all      E[beta^2] over every segregating mutation (fixed sites excluded).
           Under neutrality this must equal the singleton value, so for the
           NULL arm sing / all is a clean extraction test.  For SELECTED arms
           it is NOT kappa: older segregating variants have been filtered
           harder than singletons, so this denominator is too small and the
           ratio is pushed above kappa.

Reading the output:
  young / nominal ~ 1    the effect sizes were drawn and stored correctly.
  null  sing / young ~ 1 the extraction is unbiased.  If nominal says 1.017
                         but this says 1.000, the anomaly was the denominator.
  null  sing / young ~ 1.017 and young / nominal ~ 1
                         singletons genuinely differ from new mutations under
                         neutrality, which cannot happen biologically: an
                         extraction bug to hunt for.

All three ratios are formed from POOLED sums (sum of beta^2 / count across
runs), not by averaging per-run ratios: the young denominator has only a few
hundred mutations per run, and averaging 1/denominator over runs would add a
bias of roughly CV^2 ~ 0.3%, the size of the effects being measured.  SEs are
from a bootstrap over replicate runs, so they include between-run variance.
"""
import sys, os, glob, re, argparse
import numpy as np
import tskit

#  arm          V_S   MEAN_EFFECT  SD_EFFECT  NEUTRAL      (from armA_jobs.tsv)
ARMS = {
    "null":       (5.0, 0.00, 0.1, True),
    "realistic":  (5.0, 0.00, 0.1, False),
    "loud":       (1.0, 0.00, 0.1, False),
    "veryloud":   (0.2, 0.00, 0.1, False),
    "skewmild":   (5.0, 0.05, 0.1, False),
    "skewstrong": (5.0, 0.10, 0.1, False),
}


def sigma2_m(arm):
    _, mean, sd, _ = ARMS[arm]
    return mean ** 2 + sd ** 2


def lambda_pred(arm):
    v_s, _, _, neutral = ARMS[arm]
    return 0.0 if neutral else sigma2_m(arm) / v_s


def extract(path):
    """(count, beta, age) for every segregating single-mutation site.

    count  derived-allele count in the samples
    beta   trait effect (SLiM selection_coeff; locus fitness is neutralised
           by mutationEffect(m2) { return 1.0; })
    age    generations before the end of the run
    Fixed sites (count == number of samples) are excluded.
    """
    ts = tskit.load(path)
    n_samp = ts.num_samples
    tick = ts.metadata.get("SLiM", {}).get("tick") if isinstance(ts.metadata, dict) else None
    ks, bs, ages = [], [], []
    for v in ts.variants():
        if len(v.site.mutations) != 1:                 # skip stacked sites
            continue
        mut = v.site.mutations[0]
        ml = mut.metadata["mutation_list"]
        if len(ml) != 1:
            continue
        k = int(np.count_nonzero(v.genotypes))
        if k == 0 or k == n_samp:                       # absent or fixed
            continue
        age = mut.time
        if not np.isfinite(age) and tick is not None:   # fall back to SLiM origin tick
            age = tick - ml[0]["slim_time"]
        ks.append(k)
        bs.append(ml[0]["selection_coeff"])
        ages.append(age)
    return np.array(ks), np.array(bs, dtype=float), np.array(ages, dtype=float)


def arm_of(path):
    m = re.match(r"([a-z]+)_rep\d+", os.path.basename(path))
    return m.group(1) if m and m.group(1) in ARMS else None


def ratios(R):
    """Pooled ratios from per-run sums R (runs x 7 columns)."""
    n_s, S_s, n_a, S_a, n_y, S_y, _ = R.sum(axis=0)
    e_s = S_s / n_s
    e_a = S_a / n_a if n_a > 0 else np.nan
    e_y = S_y / n_y if n_y > 0 else np.nan
    return e_s, e_a, e_y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("directory", nargs="?", default=".")
    ap.add_argument("--out", default=None)
    ap.add_argument("--young-age", type=float, default=1.0,
                    help="max age (generations) counted as 'new'; default 1")
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    files = sorted(glob.glob(os.path.join(a.directory, "*_rep*.trees")))
    if not files:
        sys.exit(f"no *_rep*.trees found in {a.directory}")

    # per-run sums: n_sing, S_sing, n_all, S_all, n_young, S_young, Sbeta_young
    runs = {arm: [] for arm in ARMS}

    for path in files:
        arm = arm_of(path)
        if arm is None:
            continue
        try:
            k, b, age = extract(path)
        except Exception as e:
            print(f"  skipped {os.path.basename(path)}: {e}", file=sys.stderr)
            continue
        if k.size == 0:
            continue
        b2 = b ** 2
        sing = k == 1
        young = age <= a.young_age
        row = [sing.sum(), b2[sing].sum(), k.size, b2.sum(),
               young.sum(), b2[young].sum(), b[young].sum()]
        runs[arm].append(row)
        s2m = sigma2_m(arm)
        es = b2[sing].mean() if sing.any() else np.nan
        ey = b2[young].mean() if young.any() else np.nan
        print(f"  {os.path.basename(path):26s} sing {sing.sum():5d}  young {young.sum():5d}  "
              f"sing/nom {es/s2m:.4f}  sing/all {es/b2.mean():.4f}  "
              f"sing/young {es/ey:.4f}", file=sys.stderr)

    rng = np.random.default_rng(a.seed)
    rows = []
    print(f"\nyoung = mutations aged <= {a.young_age:g} generations\n")
    hdr = (f"{'arm':<11}{'reps':>5}{'sing':>9}{'young':>8}{'lambda':>8}"
           f"{'young/nom':>11}{'sing/nom':>18}{'sing/all':>18}{'sing/young':>18}")
    print(hdr)
    print("-" * len(hdr))

    def fmt(x, se):
        z = (x - 1) / se if se > 0 else np.nan
        return f"{x:7.4f}±{se:.4f}({z:+5.1f})"

    for arm in ARMS:
        if not runs[arm]:
            continue
        R = np.array(runs[arm], dtype=float)
        s2m = sigma2_m(arm)
        e_s, e_a, e_y = ratios(R)
        est = np.array([e_s / s2m, e_s / e_a, e_s / e_y, e_y / s2m])

        boots = []
        for _ in range(a.boot):
            idx = rng.integers(0, len(R), len(R))
            bs_, ba_, by_ = ratios(R[idx])
            boots.append([bs_ / s2m, bs_ / ba_, bs_ / by_, by_ / s2m])
        se = np.nanstd(np.array(boots), axis=0, ddof=1)

        n_s, n_y = int(R[:, 0].sum()), int(R[:, 4].sum())
        mean_b_young = R[:, 6].sum() / R[:, 4].sum() if n_y else np.nan
        print(f"{arm:<11}{len(R):5d}{n_s:9,d}{n_y:8,d}{lambda_pred(arm):8.4f}"
              f"{est[3]:6.3f}±{se[3]:.3f} "
              f"{fmt(est[0], se[0]):>18}{fmt(est[1], se[1]):>18}{fmt(est[2], se[2]):>18}")
        rows.append(dict(arm=arm, reps=len(R), singletons=n_s, young=n_y,
                         lam=lambda_pred(arm), sigma2_m=s2m,
                         young_over_nom=est[3], young_over_nom_se=se[3],
                         sing_over_nom=est[0], sing_over_nom_se=se[0],
                         sing_over_all=est[1], sing_over_all_se=se[1],
                         sing_over_young=est[2], sing_over_young_se=se[2],
                         mean_b_young=mean_b_young))

    print("\nratio±SE(z vs 1).  sing/young = kappa (the number you want).")
    print("young/nom ~ 1 confirms the effect sizes were drawn as configured.")
    print("For the NULL: sing/young and sing/all must both be ~1.")
    print("For SELECTED arms: sing/all is biased upward (older variants are")
    print("filtered harder), so read sing/young, not sing/all.")

    if a.out:
        keys = list(rows[0].keys())
        with open(a.out, "w") as fh:
            fh.write("\t".join(keys) + "\n")
            for r in rows:
                fh.write("\t".join(str(r[k_]) for k_ in keys) + "\n")
        print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
