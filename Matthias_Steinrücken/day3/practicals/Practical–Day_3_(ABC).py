import msprime
import demesdraw
from matplotlib import pyplot as plt
from tqdm import tqdm
import numpy as np
import tskit

# add demography of 2 populations
def repeat_simulations(mut, sample_sizes, demography, length, reco, num_simulations, seed=None):
    results = []
    for i in tqdm(range(num_simulations), desc="Running simulations"):
        if seed is not None:
            np.random.seed(seed + i)
        # Simulate 10 diploid samples under the coalescent with recombination on a 10kb region.
        demography = msprime.Demography()
        demography.add_population(name="pop1", initial_size=60000)
        demography.add_population(name="pop2", initial_size=2000)
        # instantaneous reduction of size
        demography.add_population(name="NA", initial_size=7000000)
        while True:
            Time = np.random.normal(8000, 3000, 1)[0]  # Generate a new Time value
            if Time > 1000:
                break
        demography.add_population_split(time=Time, derived=["pop1", "pop2"], ancestral="NA")
        ts = msprime.sim_ancestry(demography=demography,
                                  samples=sample_sizes,
                                  recombination_rate=reco,
                                  sequence_length=length,
                                  random_seed=np.random.randint(99999999))
        # we can add mutations
        mutated_ts = msprime.sim_mutations(ts, rate=mut, random_seed=np.random.randint(99999999))
        pop_id = {p.metadata["name"]: p.id for p in mutated_ts.populations()}
        sample_sets=[mutated_ts.samples(pop_id["pop1"]), mutated_ts.samples(pop_id["pop2"])]
        FST=mutated_ts.Fst(sample_sets)
        dxy=mutated_ts.divergence(sample_sets)
        seg=mutated_ts.segregating_sites(sample_sets)
        div=mutated_ts.diversity(sample_sets)
        TajD=mutated_ts.Tajimas_D(sample_sets)
        Time=Time
        results.append((FST, dxy, seg, div, TajD, Time))
    return results

mut= 3.5e-9
sample_sizes = {"pop1" : 50; "pop2" : 50}
length = 1000
seed = 4711
reco = 8.4e-9
num_simulations = 100

results = repeat_simulations(mut,sample_sizes, demography, length, reco, num_simulations, seed=seed)

print(results)
