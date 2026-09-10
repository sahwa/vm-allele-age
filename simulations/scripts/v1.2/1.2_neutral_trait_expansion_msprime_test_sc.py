# Here we are testing which demographic model produces the fraction of singletons that match the observed data


import msprime
import stdpopsim
import numpy as np

import copy


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



species = stdpopsim.get_species("HomSap")
base = species.get_demographic_model("OutOfAfrica_2T12").model

def scale_eur_final_size(base_demog, target_final):
    demog = copy.deepcopy(base_demog)
    demog.add_population_parameters_change(
        time=0, population="EUR",
        initial_size=target_final, growth_rate=0.0195,
    )
    demog.sort_events()
    return demog

for target in [500_000, 2_000_000, 5_000_000]:
    demog = scale_eur_final_size(base, target)
    dbg = demog.debug()
    print(f"\n=== target = {target:,} ===")
    print(dbg)

    

for target_final in [500_000, 2_000_000, 5_000_000]:
    demog = scale_eur_final_size(base, target_final)
    ts = msprime.sim_ancestry(
        samples={"EUR": 20_000},       # n_dip fixed
        demography=demog,
        sequence_length=1e7,
        recombination_rate=1e-8,
        model=[msprime.DiscreteTimeWrightFisher(duration=100),
               msprime.StandardCoalescent()],
        random_seed=1,
    )
    ts = msprime.sim_mutations(ts, rate=1.25e-8, random_seed=1)
    afs = ts.allele_frequency_spectrum(polarised = True, span_normalise = False)
    n_sing, n_sites = int(afs[1]), ts.num_sites

    mut_time = ts.tables.mutations.time
    site_of_mut = ts.tables.mutations.site
    counts = np.bincount(site_of_mut, minlength=ts.num_sites)

    ok = counts == 1
    ac = np.array([v.genotypes.sum() for v in ts.variants()])
    sing_sites = ok & (ac == 1)
    ages = mut_time[np.isin(site_of_mut, np.flatnonzero(sing_sites))]

    print(
        f"target={target_final:>9,} sites={n_sites:7d} "
        f"sing={n_sing:7d} frac={n_sing/n_sites:.3f} "
        f"age_q={np.percentile(ages, [50,90,95,99,99.9,99.99]).round(0)}"
    )


for target_final in [350_000, 450_000, 550_000]:
    demog = scale_eur_final_size(base, target_final)
    ts = msprime.sim_ancestry(
        samples={"EUR": 20_000},       # n_dip fixed
        demography=demog,
        sequence_length=1e7,
        recombination_rate=1e-8,
        model=[msprime.DiscreteTimeWrightFisher(duration=100),
               msprime.StandardCoalescent()],
        random_seed=1,
    )
    ts = msprime.sim_mutations(ts, rate=1.25e-8, random_seed=1)
    afs = ts.allele_frequency_spectrum(polarised = True, span_normalise = False)
    n_sing, n_sites = int(afs[1]), ts.num_sites

    mut_time = ts.tables.mutations.time
    site_of_mut = ts.tables.mutations.site
    counts = np.bincount(site_of_mut, minlength=ts.num_sites)

    ok = counts == 1
    ac = np.array([v.genotypes.sum() for v in ts.variants()])
    sing_sites = ok & (ac == 1)
    ages = mut_time[np.isin(site_of_mut, np.flatnonzero(sing_sites))]

    print(
        f"target={target_final:>9,} sites={n_sites:7d} "
        f"sing={n_sing:7d} frac={n_sing/n_sites:.3f} "
        f"age_q={np.percentile(ages, [50,90,95,99,99.9,99.99]).round(0)}"
    )


