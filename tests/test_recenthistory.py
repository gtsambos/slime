
import slime
import msprime, tskit
import unittest
import re
import random
import numpy as np
import utils

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

    def test_generation_order(self):
        scr = self.scr.dump_script()
        test_regex = re.compile("^\d+", flags=re.MULTILINE)
        generations = test_regex.findall(scr)
        generations = [int(gen) for gen in generations]
        self.assertEqual(generations, sorted(generations))

    def test_ordering_inside_events(self):
        config = msprime.PopulationConfiguration(0, 100, growth_rate = 0)
        scr = slime.RecentHistory(final_gen = 5, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.add_event((3, 'late'), 'second event')
        scr.add_event((3, 'late'), 'first event', insert_at_start=True)
        scr.add_event((3, 'late'), 'third event')
        script = scr.dump_script()
        event = '3 late(){'
        event_loc = script.find(event)
        end_of_scr = script[event_loc:]
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

    def test_delete_time(self):
        scr = self.scr
        scr.delete_time('3 late')
        scr.delete_time('4 early')
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '1 late', '2 late', '5 late'])
        scr.add_continuous_process((3,4), event='to delete')
        scr.delete_time('3:4 {')
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '1 late', '2 late', '5 late'])

    def test_all_generations_and_times(self):
        config = msprime.PopulationConfiguration(0, 10, growth_rate = 0)
        scr = slime.RecentHistory(final_gen = 5, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.add_event_over_time_range(start = 3, end = 4, event = 'continuous event') 
        # scr.print_script()
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '1 late', '2 late', '3:4', '5 late'])
        self.assertEqual(scr.all_generations_and_times(return_as_strings=False),
            [(1, 'early'), (1, 'late'), (2, 'late'), (3,4), (5, 'late')])

        scr = slime.RecentHistory(final_gen = 23, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.add_event((5, 'late'), 'event in gen 5')
        scr.add_continuous_process((3, 20), event = 'event')
        scr.add_continuous_process((4, 10), event = 'another event')
        # print(scr.all_generations_and_times())
        # scr.print_script()


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

    def test_add_continuous_processes(self):
        # scr = slime.RecentHistory(final_gen = 10, chrom_length = 10)
        config = msprime.PopulationConfiguration(0, 10, growth_rate = 0)
        scr = slime.RecentHistory(final_gen = 23, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.add_event((5, 'late'), 'event in gen 5')
        scr.add_continuous_process((3, 20), event = 'event')
        scr.add_continuous_process((4, 10), event = 'another event')
        # scr.print_script()
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '1 late', '2 late', '3 early', '4 early', '5 early', '5 late',
            '6:10', '11:20', '23 late'])
        self.assertEqual(scr.all_events_at_a_given_time(time='3 early'),
            ['event'])
        self.assertEqual(scr.all_events_at_a_given_time(time='4 early'),
            ['event', 'another event'])
        self.assertEqual(scr.all_events_at_a_given_time(time='5 early'),
            ['event', 'another event'])
        self.assertEqual(scr.all_events_at_a_given_time(time='5 late'),
            ['event in gen 5'])
        self.assertEqual(scr.all_events_at_a_given_time(time='6:10'),
            ['event', 'another event'])
        self.assertEqual(scr.all_events_at_a_given_time(time='11:20'),
            ['event'])

        scr.add_continuous_process((10,19), event='yet another event')
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '1 late', '2 late', '3 early', '4 early', '5 early', '5 late',
            '6:9', '10 early', '11:19', '20 early', '23 late'])
        self.assertEqual(scr.all_events_at_a_given_time(time='10 early'),
            ['event', 'another event', 'yet another event'])
        self.assertEqual(scr.all_events_at_a_given_time(time='11:19'),
            ['event', 'yet another event'])
        self.assertEqual(scr.all_events_at_a_given_time(time='20 early'),
            ['event'])

    def test_all_events_at_a_given_time(self):
        scr = self.scr
        scr.all_events_at_a_given_time('3 late')
        self.assertEqual(scr.all_events_at_a_given_time('3 late'),
            ['first event', 'second event', 'third event'])

    def test_break_up_range(self):
        config = msprime.PopulationConfiguration(0, 10, growth_rate = 0)
        config1 = msprime.PopulationConfiguration(0, 10, growth_rate = 0.01)
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

        scr = slime.RecentHistory(final_gen = 10, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.add_continuous_process((2, 8), event = 'continuous')
        scr.add_event((7, 'early'), event = 'middle')
        # scr.print_script()
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '1 late', '2 early', '2 late', '3:6', '7 early',
            '8 early', '10 late'])

        scr = slime.RecentHistory(final_gen = 11, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.add_continuous_process((3, 11), event = 'continuous')
        scr.add_continuous_process((6, 11), event = 'continuous2')
        scr.add_continuous_process((5, 11), event = 'continuous3')
        # scr.print_script()
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '1 late', '2 late', '3:4', '5 early', '6:10',
            '11 early', '11 late'])
        self.assertEqual(scr.all_events_at_a_given_time('3:4'), ['continuous'])
        self.assertEqual(scr.all_events_at_a_given_time('5 early'),
            ['continuous', 'continuous3'])
        self.assertEqual(scr.all_events_at_a_given_time('6:10'),
            ['continuous', 'continuous2', 'continuous3'])
        self.assertEqual(scr.all_events_at_a_given_time('11 early'),
            ['continuous', 'continuous2', 'continuous3'])

    def test_break_up_initial_growth(self):
        config0 = msprime.PopulationConfiguration(0, 10, growth_rate = 0)
        config1 = msprime.PopulationConfiguration(0, 10, growth_rate = 0.1)
        scr = slime.RecentHistory(final_gen = 10, chrom_length = 10,
            reference_configs = [config1, config1], adm_configs = config0,
            prop = [0.4, 0.6])
        # scr.print_script()

    def test_time_already_in_script(self):
        config = msprime.PopulationConfiguration(0, 10, growth_rate = 0)
        scr = slime.RecentHistory(final_gen = 21, chrom_length = 10,
            reference_configs = [config, config, config], adm_configs = config,
            prop = [0.5, 0.4, 0.1])
        # scr.print_script()
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '1 late', '2 late', '21 late'])

        config1 = msprime.PopulationConfiguration(0, 10, growth_rate = 0.1)
        scr = slime.RecentHistory(final_gen = 11, chrom_length = 10,
            reference_configs = [config1, config], adm_configs = config1,
            prop = [0.5, 0.5])
        # scr.print_script()


class TestDemographyConfig(unittest.TestCase):

    def test_basic_input(self):
        config = msprime.PopulationConfiguration(0, 10, growth_rate = 0)
        scr = slime.RecentHistory(final_gen = 5, chrom_length = 10,
            reference_configs = [config, config], adm_configs = config,
            prop = [0.5, 0.5])
        scr.run_slim(verbose=False)
        # scr.print_script()
    def test_nonconstant_recombination(self):
        scr = utils.basic_two_way_admixture(length=9, 
            rho = "examples/recomb_map_chromlength10.txt")
        scr.run_slim(verbose=False)

    def test_add_demographic_events(self):
        scr = utils.basic_two_way_admixture(prop=[.5,.5])
        change1 = msprime.PopulationParametersChange(3, growth_rate =  .5, population_id = 0)
        change2 = msprime.PopulationParametersChange(6, growth_rate =  .6, population_id = 0)
        change3 = msprime.PopulationParametersChange(2, growth_rate =  .3, population_id = 1)
        change4 = msprime.PopulationParametersChange(8, growth_rate =  0, population_id = 1, initial_size = 20)
        scr.add_demographic_events([change3, change1, change2, change4])
        # scr.print_script()
        scr.run_slim(verbose=False)

    def test_delete_event(self):
        scr = utils.basic_two_way_admixture(prop=[.5,.5], final_gen=10)
        scr.delete_event("simulationFinished")
        self.assertEqual(len(scr.all_events_at_a_given_time('10 late')), 1)
        change1 = msprime.PopulationParametersChange(3, growth_rate =  .5, population_id = 0)
        scr.add_demographic_events([change1])
        scr.delete_event("p0.setSubpopulationSize", time = (3, 9))
        self.assertEqual(len(scr.all_events_at_a_given_time("3:9")), 0)
        scr.run_slim(verbose=False)

    def test_size_changes(self):
        config = msprime.PopulationConfiguration(sample_size=0, initial_size=10)
        config1 = msprime.PopulationConfiguration(sample_size=0, initial_size=10, growth_rate=0.1)
        script = slime.RecentHistory(final_gen=15, chrom_length=10,
            reference_configs=[config, config1], adm_configs=config1,
            prop=[0.3,0.7])
        script.add_size_change(msprime.PopulationParametersChange(time=10,
                initial_size=15, population=2))
        script.add_size_change(msprime.PopulationParametersChange(time=10,
                initial_size=15, population=0))
        # script.print_script()
        e_gen10 = script.all_events_at_a_given_time('10 early')
        self.assertEqual(len(e_gen10), 3)
        assert 'p0.setSubpopulationSize(15)' in e_gen10
        assert 'p2.setSubpopulationSize(15)' in e_gen10
        assert 'p1.setSubpopulationSize(asInteger(p2.individualCount * exp(0.100000)))' not in e_gen10
        script.run_slim(verbose=False)

    def test_migration_rate_changes(self):
        scr = utils.basic_two_way_admixture(rho=0)
        scr.add_migration_rate(time=(5, 'late'), rates=[.01,.02])
        # scr.print_script()
        scr.run_slim(verbose=False)  

        scr = utils.basic_two_way_admixture(rho=0)
        scr.add_migration_rate(rates=[.05,.1])
        # scr.print_script()
        scr.run_slim(verbose=False)  

    def test_mass_migration(self):
        scr = utils.basic_two_way_admixture(rho=.001)
        scr.add_mass_migration(prop=[.2,.2], time=(9, 'late'))
        # scr.print_script()
        scr.run_slim(verbose=False)

    def test_mutations(self):
        muts = slime.MutationTypes(mutation_rate=.005, selection_coeffs=[0.2,0.6],
            proportions=[.5,.5], dominance_coeffs=[.7, .45])
        config = msprime.PopulationConfiguration(sample_size=0, initial_size=10)
        script = slime.RecentHistory(final_gen=20, chrom_length=100,
            reference_configs=[config, config], adm_configs=config,
            prop=[0.3,0.7], mutations=muts)
        # script.print_script()
        script.run_slim(verbose=False)
        

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

    def test_recenthistory_defaults(self):
        ref_pops = [msprime.PopulationConfiguration(sample_size=10, initial_size=100),
               msprime.PopulationConfiguration(sample_size=10, initial_size=100)]
        adm_pop = msprime.PopulationConfiguration(sample_size=20, initial_size=50)
        adm_prop = [0.3, 0.7]
        gens = 20 # 100
        script = slime.RecentHistory(final_gen=gens, chrom_length=10,
                recombination=1e-1, reference_configs=ref_pops, adm_configs=adm_pop,
                prop=adm_prop)
        script.run_slim(verbose=False)
        # script.print_script()

    def test_recenthistory_constant_growth(self):
        ref_pops = [msprime.PopulationConfiguration(sample_size=10, initial_size=100),
               msprime.PopulationConfiguration(sample_size=10, initial_size=100)]
        adm_pop = msprime.PopulationConfiguration(sample_size=20, initial_size=50, growth_rate=0.1)
        adm_prop = [0.3, 0.7]
        gens = 20 # 100
        script = slime.RecentHistory(final_gen=gens, chrom_length=10,
                recombination=1e-1, reference_configs=ref_pops, adm_configs=adm_pop,
                prop=adm_prop)
        script.run_slim(verbose=False)

    def test_recenthistory_instant_size_change(self):
        ref_pops = [msprime.PopulationConfiguration(sample_size=10, initial_size=100),
                msprime.PopulationConfiguration(sample_size=10, initial_size=100)]
        adm_pop = msprime.PopulationConfiguration(sample_size=20, initial_size=50)
        adm_prop = [0.3, 0.7]
        gens = 20 # 100
        ch = msprime.PopulationParametersChange(time=11, initial_size=25, population_id=2)
        script = slime.RecentHistory(final_gen=gens, chrom_length=10,
            recombination=1e-1, reference_configs=ref_pops, adm_configs=adm_pop,
            prop=adm_prop)
        script.add_size_change(ch)
        script.run_slim(verbose=False)

    def test_recenthistory_migration_rates(self):
        config = msprime.PopulationConfiguration(sample_size=10, initial_size=100)
        script = slime.RecentHistory(final_gen=20, chrom_length=10, recombination=.001,
            reference_configs=[config, config], adm_configs=config, prop=[0.3,0.7])
        script.add_migration_rate(rates=[.01,.02])
        script.add_migration_rate(rates=[0.0, 0.0], time=(8, 'late'))
        # script.print_script()
        script.run_slim(verbose=False)

    def test_recenthistory_mass_migration(self):
        config = msprime.PopulationConfiguration(sample_size=10, initial_size=100)
        script = slime.RecentHistory(final_gen=20, chrom_length=10, recombination=.001,
            reference_configs=[config, config], adm_configs=config, prop=[0.3,0.7])
        script.add_mass_migration(prop=[.1, .1], time=(9, 'late'))
        script.run_slim(verbose=False)

    def test_mutations(self):
        muts = slime.MutationTypes(mutation_rate=.005, selection_coeffs=[0.2,0.6,.3],
            proportions=[.5,.4, .1], dominance_coeffs=[.7, .45, .57])
        config = msprime.PopulationConfiguration(sample_size=0, initial_size=10)
        script = slime.RecentHistory(final_gen=20, chrom_length=100,
            reference_configs=[config, config], adm_configs=config,
            prop=[0.3,0.7], mutations=muts)
        # script.print_script()
        script.run_slim(verbose=False)

class TestDemography(unittest.TestCase):
    """
    Runs SLiM on randomly-generated scripts under various models of demography.
    This is mainly to ensure that the generated scripts are actually valid.
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

    def test_constantly_growing_populations(self):
        ref_configs = []
        for pop in range(0,3):
            ref_configs.append(msprime.PopulationConfiguration(0, initial_size=10, 
                growth_rate=random.randint(1,10)/100))
        adm_config = msprime.PopulationConfiguration(0, initial_size=10,
            growth_rate = random.randint(1,10)/100)
        script = slime.RecentHistory(final_gen=17, chrom_length=10,
            reference_configs=ref_configs, adm_configs=adm_config,
            prop=[0.2,0.3,0.5])
        self.assertEqual(script.all_generations_and_times(),
            ['1 early', '1 late', '2 early', '2 late', '3:16', '17 early', '17 late'])
        # script.print_script()
        script.run_slim(verbose=False)

    def test_size_changes(self):
        config = msprime.PopulationConfiguration(sample_size=0, initial_size=10)
        config1 = msprime.PopulationConfiguration(sample_size=0, initial_size=10, growth_rate=0.1)
        script = slime.RecentHistory(final_gen=30, chrom_length=10,
            reference_configs=[config, config1], adm_configs=config1,
            prop=[0.3,0.7])
        for pop in [0, 1, 2]:
            times = random.sample(range(3,29), 3)
            new_sizes = random.sample(range(10, 25), 3)
            for i in range(0,3):
                script.add_size_change(msprime.PopulationParametersChange(time=times[i],
                    initial_size=new_sizes[i], population=pop))
        # script.print_script()
        script.run_slim(verbose=False)

    def test_random_recombination_map(self):
        config = msprime.PopulationConfiguration(sample_size=0, initial_size=10)
        utils.random_recombination_map(sequence_length=100, no_rows=20, 
            filename='examples/test.recomb')
        script = slime.RecentHistory(final_gen=20, chrom_length=100,
            reference_configs=[config, config], adm_configs=config,
            prop=[0.3,0.7], recombination='examples/test.recomb')
        script.run_slim(verbose=False)

    def test_migration_rate_changes(self):
        scr = utils.basic_two_way_admixture(rho=.001)
        # Add in some random migrations.
        for m in range(0, 10):
            rates = random.sample(range(1, 100), 2)
            time = random.randint(2, 19)
            rates = [p/200 for p in rates]
            scr.add_migration_rate(rates=rates, time = (time, 'late'))
        # scr.print_script()
        scr.run_slim(verbose=False)

    def test_mass_migrations(self):
        scr = utils.basic_two_way_admixture(rho=.001)
        # Add in some random migrations.
        for m in range(0, 10):
            rates = random.sample(range(1, 100), 2)
            time = random.randint(2, 19)
            rates = [p/200 for p in rates]
            scr.add_mass_migration(prop=rates, time = (time, 'late'))
        # scr.print_script()
        scr.run_slim(verbose=False)
