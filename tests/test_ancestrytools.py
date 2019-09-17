import slime
import unittest
import msprime, tskit
import collections
import numpy as np

class TestAncestryTableClass(unittest.TestCase):

    def test_initialize_works(self):
        tab = slime.AncestryTable()
        self.assertEqual(tab.num_rows, 0)

    def test_add_row(self):
        tab = slime.AncestryTable()
        tab.add_row(0, 1, 0, 1, 5)
        self.assertEqual(tab.num_rows, 1)

    def test_add_row_no_ancestors(self):
        tab = slime.AncestryTable()
        tab.add_row(0, 1, 0, 1)
        self.assertEqual(tab.num_rows, 1)
        self.assertEqual(tab.ancestor, [-1])

    def test_set_columns(self):
        tab = slime.AncestryTable()
        tab.set_columns(
            left = [0, 0, 0, 0],
            right = [1, 1, 1, 1],
            population = [0, 1, 0, 1],
            child = [0, 1, 2, 3],
            ancestor = [4, 5, 6, 7])
        self.assertEqual(tab.num_rows, 4)

    def test_set_columns_no_ancestors(self):
        tab = slime.AncestryTable()
        tab.set_columns(
            left = [0, 0, 0, 0],
            right = [1, 1, 1, 1],
            population = [0, 1, 0, 1],
            child = [0, 1, 2, 3])
        self.assertEqual(tab.num_rows, 4)  
        self.assertEqual(list(tab.ancestor), [-1, -1, -1, -1])


class TestGetAncestryTables(unittest.TestCase):

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

    samples_ex = ts_ex.samples()
    populations_ex = ts_ex.populations()

    # Here's the ancestry table (unsorted) 
    #     id  left        right       parent  child
    # 0   898.33110615    1000.00000000   1   3
    # 1   898.33110615    1000.00000000   1   5
    # 2   0.00000000  1000.00000000   0   2
    # 3   0.00000000  898.33110615    0   5
    # 5   0.00000000  1000.00000000   0   0
    # 6   0.00000000  1000.00000000   0   4
    # 9   0.00000000  898.33110615    0  1
    # 10  0.00000000  898.33110615    0  3
    # 11  898.33110615    1000.00000000   0  1

    def test_simple_case(self):
        slime.get_ancestry_table(self.ts_ex, self.samples_ex, self.populations_ex)


