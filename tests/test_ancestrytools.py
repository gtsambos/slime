import slime
import unittest
import msprime, tskit
import collections
import numpy as np
import io

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

    # Tree sequence  | Population
    #     4          | 2
    #    / \         |
    #   /   3        | 1
    #  /   / \       |
    # 0   1   2      | 0
    nodes0 = io.StringIO("""\
    id      is_sample   population      time
    0       1       0               0.00000000000000
    1       1       0               0.00000000000000
    2       1       0               0.00000000000000
    3       0       1               1.00000000000000
    4       0       2               2.00000000000000
    """)
    edges0 = io.StringIO("""\
    id      left            right           parent  child
    0       0.00000000      1.00000000      3       1,2
    1       0.00000000      1.00000000      4       0,3
    """)

    # Tree sequence                                       | Population
    #          9                        10                | 3
    #         / \                      / \                |
    #        /   \                    /   8               | 2
    #       /     \                  /   / \              |
    #      7       \                /   /   \             | 1
    #     / \       6              /   /     6            | 0
    #    /   5     / \            /   5     / \           | 0
    #   /   / \   /   \          /   / \   /   \          |
    #  4   0   1 2     3        4   0   1 2     3         | 0
    #
    # 0 ------------------ 0.5 ------------------ 1.0

    nodes1 = io.StringIO("""\
    id      is_sample   population      time
    0       1       0               0.00000000000000
    1       1       0               0.00000000000000
    2       1       0               0.00000000000000
    3       1       0               0.00000000000000
    4       1       0               0.00000000000000
    5       0       0               0.14567111023387
    6       0       0               0.21385545626353
    7       0       1               0.43508024345063
    8       0       2               0.60156352971203
    9       0       3               0.90000000000000
    10      0       3               1.20000000000000
    """)
    edges1 = io.StringIO("""\
    id      left            right           parent  child
    0       0.00000000      1.00000000      5       0,1
    1       0.00000000      1.00000000      6       2,3
    2       0.00000000      0.50000000      7       4,5
    3       0.50000000      1.00000000      8       5,6
    4       0.00000000      0.50000000      9       6,7
    5       0.50000000      1.00000000      10      4,8
    """)

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
    populations_ex = [0, 1]

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
        tab = slime.get_ancestry_table(self.ts_ex, list(self.populations_ex))
        # print(tab)


    def test_simple_example(self):
        ts = tskit.load_text(nodes=self.nodes0, edges=self.edges0, strict=False)
        tab = slime.get_ancestry_table(ts, [1, 2])
        self.assertEqual(list(tab.left), [0, 0, 0])
        self.assertEqual(list(tab.right), [1, 1, 1])
        self.assertEqual(list(tab.population), [2, 1, 1])
        self.assertEqual(list(tab.child), [0, 1, 2])

    def test_one_population(self):
        ts = tskit.load_text(nodes=self.nodes0, edges=self.edges0, strict=False)
        tab = slime.get_ancestry_table(ts, [1, 3])
        self.assertEqual(list(tab.left), [0, 0])
        self.assertEqual(list(tab.right), [1, 1])
        self.assertEqual(list(tab.population), [1, 1])
        self.assertEqual(list(tab.child), [1, 2])

    def test_no_rows(self):
        ts = tskit.load_text(nodes=self.nodes0, edges=self.edges0, strict=False)
        tab = slime.get_ancestry_table(ts, samples=[0], populations=[1])
        self.assertEqual(tab.num_rows, 0)

    def test_no_ancestors(self):
        ts = tskit.load_text(nodes=self.nodes0, edges=self.edges0, strict=False)
        with self.assertRaises(ValueError):
            slime.get_ancestry_table(ts, samples=[0], populations=[3])

    def test_some_samples_have_relevant_population_labels(self):
        # Tree sequence  | Population
        #     4          | 2
        #    / \         |
        #   /   3        | 1
        #  /   / \       |
        # 0   1   2      | 0, 1, 2
        nodes = io.StringIO("""\
        id      is_sample   population      time
        0       1       0               0.00000000000000
        1       1       1               0.00000000000000
        2       1       2               0.00000000000000
        3       0       1               1.00000000000000
        4       0       2               2.00000000000000
        """)
        edges = io.StringIO("""\
        id      left            right           parent  child
        0       0.00000000      1.00000000      3       1,2
        1       0.00000000      1.00000000      4       0,3
        """)
        ts = tskit.load_text(nodes=nodes, edges=edges, strict=False)
        tab = slime.get_ancestry_table(ts, populations=[1, 2])
        self.assertEqual(list(tab.left), [0, 0, 0])
        self.assertEqual(list(tab.right), [1, 1, 1])
        self.assertEqual(list(tab.population), [2, 1, 2])
        self.assertEqual(list(tab.child), [0, 1, 2])

    def test_multiple_trees(self):
        ts = tskit.load_text(nodes=self.nodes1, edges=self.edges1, strict=False)
        tab = slime.get_ancestry_table(ts, [1, 2, 3])
        self.assertEqual(list(tab.left), [0, .5, 0, .5, 0, .5, 0, .5, 0, .5])
        self.assertEqual(list(tab.right), [.5, 1, .5, 1, .5, 1, .5, 1, .5, 1])
        self.assertEqual(list(tab.population), [1, 2, 1, 2, 3, 2, 3, 2, 1, 3])
        self.assertEqual(list(tab.child), [0, 0, 1, 1, 2, 2, 3, 3, 4, 4])

    def test_multiple_trees_edges_are_squashed(self):
        ts = tskit.load_text(nodes=self.nodes1, edges=self.edges1, strict=False)
        tab = slime.get_ancestry_table(ts, [3])
        self.assertEqual(list(tab.left), [0, 0, 0, 0, 0])
        self.assertEqual(list(tab.right), [1, 1, 1, 1, 1])
        self.assertEqual(list(tab.population), [3, 3, 3, 3, 3])
        self.assertEqual(list(tab.child), [0, 1, 2, 3, 4])

    def test_multiple_trees_some_missing_segments(self):
        ts = tskit.load_text(nodes=self.nodes1, edges=self.edges1, strict=False)
        tab = slime.get_ancestry_table(ts, [1, 2])
        self.assertEqual(list(tab.left), [0, .5, 0, .5, .5, .5, 0])
        self.assertEqual(list(tab.right), [.5, 1, .5, 1, 1, 1, .5])
        self.assertEqual(list(tab.population), [1, 2, 1, 2, 2, 2, 1])
        self.assertEqual(list(tab.child), [0, 0, 1, 1, 2, 3, 4])

    def verify(self, ts, samples, populations):
        tab = slime.get_ancestry_table(ts, populations=populations, samples=samples)
        # Loop through the rows of the ancestral branch table.
        for row in range(0, tab.num_rows):
            current_sample = tab.child[row]
            current_left = tab.left[row]
            current_right = tab.right[row]
            current_pop = ts.tables.nodes.population[current_sample]
            if current_pop in populations:
                self.assertEqual(tab.population[row], current_pop)
            else:
                for tree in ts.trees():
                    if tree.interval[0] >= current_right:
                        break
                    while tree.interval[1] <= current_left:
                        tree.next()
                    # Check that that most recent node from a relevant population has the
                    # same population ID as listed in the ancestor table.
                    par = tree.get_parent(current_sample)
                    ancestor_pop = ts.tables.nodes.population[par]
                    while ancestor_pop not in populations:
                        par = tree.get_parent(par)
                        ancestor_pop = ts.tables.nodes.population[par]
                    self.assertEqual(tab.population[row], ancestor_pop)


    def test_simple_case(self):
        self.verify(self.ts_ex, list(self.ts_ex.samples()), list(self.populations_ex))
        self.verify(self.ts_ex, list(self.ts_ex.samples()), [0])

    def test_admixture_two_populations(self):
        pop0 = msprime.PopulationConfiguration(
            sample_size=20, initial_size = 500, growth_rate = 0.00)
        pop1 = msprime.PopulationConfiguration(
            sample_size=0, initial_size = 500, growth_rate = 0.00)
        admixture_event  = msprime.MassMigration(
            time=50, source=0, dest=1, proportion=0.5)
        migration_rate_change = msprime.MigrationRateChange(
            time=500, rate=0.01, matrix_index=(0, 1))
        ts = msprime.simulate(
                population_configurations = [pop0, pop1], 
                length = 1000, 
                demographic_events = [admixture_event, migration_rate_change], 
                random_seed = 23, 
                recombination_rate = 5e-6)
        self.verify(ts, list(ts.samples()), [1])

    def test_admixture_five_populations(self):
        pop0 = msprime.PopulationConfiguration(
            sample_size=20, initial_size = 500, growth_rate = 0.00)
        pop1 = msprime.PopulationConfiguration(
            sample_size=0, initial_size = 500, growth_rate = 0.00)
        pop2 = msprime.PopulationConfiguration(
            sample_size=0, initial_size = 500, growth_rate = 0.00)
        pop3 = msprime.PopulationConfiguration(
            sample_size=0, initial_size = 500, growth_rate = 0.00)
        pop4 = msprime.PopulationConfiguration(
            sample_size=0, initial_size = 500, growth_rate = 0.00)
        migration_rate_change = msprime.MigrationRateChange(
            time=500, rate=0.01)
        ts = msprime.simulate(
                population_configurations = [pop0, pop1, pop2, pop3, pop4], 
                length = 1000, 
                demographic_events = [migration_rate_change], 
                random_seed = 3536, 
                recombination_rate = 1e-6)
        self.verify(ts, list(ts.samples()), [1, 2, 3, 4])
        self.verify(ts, list(ts.samples()), [2, 4])

