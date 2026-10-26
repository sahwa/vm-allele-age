import tskit
import tsdate
import pathlib
import daiquiri

daiquiri.setup(level="DEBUG")

N=500000
L="1e+08"

MU = 1.25e-8                  # per-bp per-generation mutation rate
REC = 1e-8                    # per-bp per-generation recombination rate


DATA = f"/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.2/n{N}_L{L}"

ts = tskit.load(pathlib.Path(DATA, "sim.trees"))
afs = ts.allele_frequency_spectrum(polarised=True, span_normalise=False)
n_sing, n_sites = int(afs[1]), ts.num_sites

# singleton ages straight from the tables (no variant traversal)
tabs = ts.tables
mut_node = tabs.mutations.node
mt = tabs.mutations.time
ages_all = np.where(np.isnan(mt), tabs.nodes.time[mut_node], mt)
is_sample = np.zeros(ts.num_nodes, dtype=bool)
is_sample[ts.samples()] = True
sing_ages = ages_all[is_sample[mut_node]]


dates = tsdate.date(ts, mutation_rate = MU)
dates.dump(pathlib.Path(DATA, "sim.trees.dates"))

afs = ts.allele_frequency_spectrum(polarised=True, span_normalise=False)
n_sing, n_sites = int(afs[1]), ts.num_sites

# singleton ages straight from the tables (no variant traversal)
tabs = ts.tables
mut_node = tabs.mutations.node
mt = tabs.mutations.time
ages_all = np.where(np.isnan(mt), tabs.nodes.time[mut_node], mt)
is_sample = np.zeros(ts.num_nodes, dtype=bool)
is_sample[ts.samples()] = True
sing_ages = ages_all[is_sample[mut_node]]
