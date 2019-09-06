import slime
import unittest
import msprime, tskit

class TestAncestryTables(unittest.TestCase):

    # A simple example to test on.
    pop_configs_ex = [
    msprime.PopulationConfiguration(sample_size=3, initial_size = 500, growth_rate = 0),
    msprime.PopulationConfiguration(sample_size=3, initial_size = 500, growth_rate = 0)]

    migration_rates_ex = np.array([
    [0, 0.05],
    [0.02, 0]])

    demographic_events_ex = [
        msprime.MigrationRateChange(time = 100, rate = 0.01, matrix_index=(0, 1))]

    # Simulate!
    ts_ex = msprime.simulate(
        population_configurations = pop_configs_ex, 
        migration_matrix = migration_rates_ex, length = 1000, 
        demographic_events = demographic_events_ex, recombination_rate = 1e-6,
        random_seed= 23)

    '''
    Here's the ancestry table (unsorted) 
    id  left        right       parent  child
0   898.33110615    1000.00000000   1   3
1   898.33110615    1000.00000000   1   5
2   0.00000000  1000.00000000   0   2
3   0.00000000  898.33110615    0   5
5   0.00000000  1000.00000000   0   0
6   0.00000000  1000.00000000   0   4
9   0.00000000  898.33110615    0  1
10  0.00000000  898.33110615    0  3
11  898.33110615    1000.00000000   0  1
    '''


