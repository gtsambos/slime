
import slime
import msprime, tskit
import unittest
import re
import random
import numpy as np

class TestRecentHistory(unittest.TestCase):

    config = msprime.PopulationConfiguration(0, 100, growth_rate = 0)
    scr = slime.RecentHistory(final_gen = 5, chrom_length = 10,
        reference_configs = [config, config], adm_configs = config,
        prop = [0.5, 0.5])
    scr.add_event((1, 'early'), 'start of simulation')
    scr.add_event((3, 'late'), 'second event')
    scr.add_event((3, 'late'), 'first event', insert_at_start=True)
    scr.add_event((3, 'late'), 'third event')
    scr.add_event((2, 'late'), 'the middle of the script')
    scr.add_event((5, 'late'), 'the end of the script')
    scr.add_event((4, 'early'), 'more middle of the script')

    # """
    # Tests the methods in the RecentHistory class.
    # """

    def test_initialize(self):
        scr = self.scr
        scr.initialize('end of initialize')
        scr.initialize('start of initialize', start = True) 
        # scr.print_script()

    def test_time_is_continuous(self):
        scr = self.scr
        self.assertEqual(scr.time_is_continuous((1, 'early')), False)
        self.assertEqual(scr.time_is_continuous((2, 'late')), False)
        self.assertEqual(scr.time_is_continuous((1, 3)), True)

    def test_add_continuous_process(self):
        config = msprime.PopulationConfiguration(0, 100, growth_rate = 0)
        scr = slime.RecentHistory(final_gen = 10, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.add_continuous_process((3, 8), 'continuous event')
        # scr.print_script()

    def test_generation_order(self):
        scr = self.scr.dump_script()
        test_regex = re.compile("^\d+", flags=re.MULTILINE)
        generations = test_regex.findall(scr)
        generations = [int(gen) for gen in generations]
        self.assertEqual(generations, sorted(generations))

    def test_ordering_inside_events(self):
        scr = self.scr.dump_script()
        event = '3 late(){'
        event_loc = scr.find(event)
        end_of_scr = scr[event_loc:]
        first_loc = end_of_scr.find('first event')
        second_loc = end_of_scr.find('second event')
        third_loc = end_of_scr.find('third event')
        self.assertLess(first_loc, second_loc)
        self.assertLess(second_loc, third_loc)
        # self.scr.print_script()

    def test_add_event_over_time_range(self):
        config = msprime.PopulationConfiguration(0, 10, growth_rate = 0)
        scr = slime.RecentHistory(final_gen = 10, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.add_event_over_time_range(start = 3, end = 7, event = 'continuous event')  
        # scr.print_script()      


    def test_all_generations_and_times(self):
        config = msprime.PopulationConfiguration(0, 10, growth_rate = 0)
        scr = slime.RecentHistory(final_gen = 5, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.add_event_over_time_range(start = 3, end = 4, event = 'continuous event') 
        # scr.print_script()
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '1 late', '2 late', '3:4', '5 late'])

    def test_all_early_events(self):
        scr = self.scr
        self.assertEqual(scr.all_early_events(), 
            ['1 early', '4 early'])

    def test_all_late_events(self):
        scr = self.scr
        self.assertEqual(scr.all_late_events(), 
            ['1 late', '2 late', '3 late', '5 late'])

    def test_all_early_or_late_events(self):
        scr = self.scr
        self.assertEqual(scr.all_early_or_late_events(),
            ['1 early', '1 late', '2 late', '3 late', '4 early', '5 late'])

    def test_all_continuous_processes(self):
        # scr = slime.RecentHistory(final_gen = 10, chrom_length = 10)
        config = msprime.PopulationConfiguration(0, 10, growth_rate = 0)
        scr = slime.RecentHistory(final_gen = 10, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.add_event((5, 'late'), 'event in gen 5')
        scr.add_continuous_process((3, 10), event = 'event')
        self.assertEqual(scr.all_continuous_processes(), [(3,4), (6,9)])

    def test_all_events_at_a_given_time(self):
        scr = self.scr
        scr.all_events_at_a_given_time('3 late')
        self.assertEqual(scr.all_events_at_a_given_time('3 late'),
            ['first event', 'second event', 'third event'])

    def test_break_up_range(self):
        config = msprime.PopulationConfiguration(0, 10, growth_rate = 0)
        scr = slime.RecentHistory(final_gen = 5, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.add_continuous_process((2, 4), event = 'growth')
        scr.add_event((2, 'early'), event = 'event')
        events = scr.all_generations_and_times()
        # scr.print_script()
        self.assertEqual(events, ['1 early', '1 late', '2 early', '2 late', '3:4', '5 late'])

        scr = slime.RecentHistory(final_gen = 5, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.add_continuous_process((2, 3), event = 'growth')
        scr.add_event((2, 'early'), event = 'event')
        # scr.print_script()
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '1 late', '2 early', '2 late', '3 early', '5 late'])

        scr = slime.RecentHistory(final_gen = 5, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.add_continuous_process((2, 4), event = 'growth')
        scr.add_event((4, 'early'), event = 'event')
        # scr.print_script()
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '1 late', '2 early', '2 late', '3 early', '4 early', '5 late'])

        scr = slime.RecentHistory(final_gen = 5, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.add_continuous_process((3, 4), event = 'growth')
        scr.add_event((4, 'early'), event = 'event')
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '1 late', '2 late', '3 early', '4 early', '5 late'])

        scr = slime.RecentHistory(final_gen = 10, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.add_continuous_process((2, 8), event = 'continuous')
        scr.add_event((5, 'early'), event = 'middle')
        # scr.print_script()
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '1 late', '2 early', '2 late', '3:4', '5 early', '6:8', '10 late'])

        scr = slime.RecentHistory(final_gen = 10, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.add_continuous_process((2, 8), event = 'continuous')
        scr.add_event((3, 'early'), event = 'middle')
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '1 late', '2 early', '2 late', '3 early', '4:8', '10 late'])

        scr = slime.RecentHistory(final_gen = 10, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.add_continuous_process((2, 8), event = 'continuous')
        scr.add_event((8, 'early'), event = 'middle')
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '1 late', '2 early', '2 late', '3:7', '8 early', '10 late'])


class TestDemographyConfig(unittest.TestCase):

    def test_basic_input(self):
        config = msprime.PopulationConfiguration(0, 10, growth_rate = 0)
        scr = slime.RecentHistory(final_gen = 5, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.run_slim(verbose=False)
        # scr.print_script()
    def test_nonconstant_recombination(self):
        config = msprime.PopulationConfiguration(0, 10, growth_rate = 0)
        scr = slime.RecentHistory(final_gen = 5, chrom_length = 9,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5], recombination = "examples/recomb_map_chromlength10.txt")
        # scr.print_script()
        scr.run_slim(verbose=False)

    def test_add_demographic_events(self):
        config = msprime.PopulationConfiguration(0, 10, growth_rate = 0)
        scr = slime.RecentHistory(final_gen = 10, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        change1 = msprime.PopulationParametersChange(3, growth_rate =  .5, population_id = 0)
        change2 = msprime.PopulationParametersChange(6, growth_rate =  .6, population_id = 0)
        change3 = msprime.PopulationParametersChange(2, growth_rate =  .3, population_id = 1)
        change4 = msprime.PopulationParametersChange(8, growth_rate =  0, population_id = 1, initial_size = 20)
        scr.add_demographic_events([change3, change1, change2, change4])
        # scr.print_script()
        scr.run_slim(verbose=False)

    def test_delete_event(self):
        config = msprime.PopulationConfiguration(0, 10, growth_rate = 0)
        scr = slime.RecentHistory(final_gen = 10, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.delete_event("simulationFinished")
        self.assertEqual(len(scr.all_events_at_a_given_time('10 late')), 1)
        change1 = msprime.PopulationParametersChange(3, growth_rate =  .5, population_id = 0)
        scr.add_demographic_events([change1])
        scr.delete_event("p0.setSubpopulationSize", time = (3, 9))
        self.assertEqual(len(scr.all_events_at_a_given_time("3:9")), 0)
        scr.run_slim(verbose=False)


class TestExamplesInDocs(unittest.TestCase):
    """
    Holds examples of RecentHistory shown in the documentation.
    """

    def test_quickstart(self):
        ref0_config = msprime.PopulationConfiguration(sample_size=0, initial_size=10, growth_rate=0)
        ref1_config = msprime.PopulationConfiguration(sample_size=0, initial_size=15, growth_rate=0)
        adm_config = msprime.PopulationConfiguration(sample_size=5, initial_size=10, growth_rate = 0)
        adm_props = [0.3, 0.7]
        rho = 0.1
        length = 10
        gens = 15
        script = slime.RecentHistory(final_gen=gens, chrom_length=length,
            reference_configs=[ref0_config, ref1_config], adm_configs=adm_config,
            prop=adm_props)
        script.run_slim(verbose=False)
        # Check output.
        ts = tskit.load("recent-history.trees")
        self.assertEqual(len(ts.samples(0)), 10*2)
        self.assertEqual(len(ts.samples(1)), 15*2)
        self.assertEqual(len(ts.samples(2)), 10*2)
        root0 = ts.first().roots[0]
        self.assertEqual(ts.first().time(root0), 15)


class TestDemography(unittest.TestCase):
    """
    Runs SLiM on randomly-generated scripts under various models of demography.
    """

    def test_many_populations(self):
        config = msprime.PopulationConfiguration(sample_size=0, initial_size=10, growth_rate=0)
        num_pops = random.randint(5, 15)
        ref_configs = [config for i in range(0, num_pops)]
        props = np.zeros(num_pops)
         # Ensure rounding error doesn't cause an error.
         # This is giving me quite a bit of grief!!
        while sum(props) != 1 or any(props) == 0:
            props = np.random.dirichlet([1 for i in range(0, num_pops)])
            props = np.round(props, decimals=2)
            props[-1] = np.round(1 - sum(props[:-1]), decimals=2)
        # props = np.round(props, decimals=2)
        assert sum(props) == 1
        script = slime.RecentHistory(final_gen=20, chrom_length=10,
            reference_configs=ref_configs, adm_configs=config,
            prop=props)
        # script.print_script()
        script.run_slim(verbose=False)
        # Check output
        ts = tskit.load("recent-history.trees")
        for pop in range(0, num_pops + 1):
            self.assertEqual(len(ts.samples(pop)), 10*2)
        root0 = ts.first().roots[0]
        self.assertEqual(ts.first().time(root0), 20)

    # def test_constantly_growing_populations(self):
    #     ref_configs = []
    #     for pop in range(0,3):
    #         ref_configs.append(msprime.PopulationConfiguration(0, initial_size=10, 
    #             growth_rate=random.randint(1,10)/100))
    #     adm_config = msprime.PopulationConfiguration(0, initial_size=10,
    #         growth_rate = random.randint(1,10)/100)
    #     script = slime.RecentHistory(final_gen=20, chrom_length=10,
    #         reference_configs=ref_configs, adm_configs=adm_config,
    #         prop=[0.2,0.3,0.5])
    #     script.print_script()
    #     script.run_slim(verbose=False)


