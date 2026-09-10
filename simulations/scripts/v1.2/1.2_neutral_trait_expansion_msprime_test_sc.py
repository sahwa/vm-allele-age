import msprime
import stdpopsim
import numpy as np

MU = 1.25e-8
L = 1e7
REC = 1e-8

species = stdpopsim.get_species("HomSap")

for name, pop in [("OutOfAfrica_3G09", "CEU"), ("OutOfAfrica_2T12", "EUR")]:
    demog = species.get_demographic_model(name).model

    for n_dip in [1_000, 5_000, 20_000]:
        ts = msprime.sim_ancestry(
            samples = {pop: n_dip},
            demography = demog,
            sequence_length = L,
            recombination_rate = REC,
            model = [msprime.DiscreteTimeWrightFisher(duration=100),
                     msprime.StandardCoalescent()],
            random_seed = 1
        )

        ts = msprime.sim_mutations(ts, rate=MU, random_seed=1)
        afs = ts.allele_frequency_spectrum(polarised = True, span_normalise = False)
        n_sing, n_sites = int(afs[1]), ts.num_sites

        mut_time = ts.tables.mutations.time
        site_of_mut = ts.tables.mutations.site
        counts = np.bincount(site_of_mut, minlength=ts.num_sites)

        ok = counts == 1
        ac = np.array([v.genotypes.sum() for v in ts.variants()])
        sing_sites = ok & (ac == 1)
        ages = mut_time[np.isin(site_of_mut, np.flatnonzero(sing_sites))]

        print(f"{name:18s} n={n_dip:6d} sites={n_sites:7d} "
              f"sing={n_sing:7d} frac={n_sing/n_sites:.3f} "
              f"age_q={np.percentile(ages, [50, 90, 95, 99, 99.9, 99.99]).round(0)}")
