
import slime
import msprime
import unittest
import re

class TestRecentHistory(unittest.TestCase):

    scr = slime.RecentHistory(final_gen = 5, chrom_length = 10)
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
        scr = slime.RecentHistory(final_gen = 10, chrom_length = 10)
        scr.add_continuous_process((3, 8), 'continuous event')
        # scr.print_script()

    def test_generation_order(self):
        scr = self.scr.dump_script()
        test_regex = re.compile("^\d+", flags=re.MULTILINE)
        generations = test_regex.findall(scr)
        generations = [int(gen) for gen in generations]
        self.assertEqual(generations, sorted(generations))

    def test_time_order(self):
        scr = self.scr.dump_script()
        gen_regex = re.compile("^\d+", flags=re.MULTILINE)
        event_regex = re.compile(r"\d+ \bearly\b|\d+ \blate\b")
        generations = gen_regex.findall(scr)
        generations = [int(gen) for gen in generations]
        events = event_regex.findall(scr)
        for gen in generations:
            early_str = "%i early" % gen
            late_str = "%i late" % gen
            early_regex = re.compile(early_str)
            late_regex = re.compile(late_str)
            early = early_regex.findall(scr)
            late = late_regex.findall(scr)
            self.assertLess(len(early), 2)
            self.assertLess(len(late), 2)
            if len(early) == 1 and len(late) == 1:
                self.assertEqual(events.index(early) + 1, events.index(late))

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
        # scr.print_script()

    def test_add_event_over_time_range(self):
        scr = slime.RecentHistory(final_gen = 5, chrom_length = 10)
        scr.add_event_over_time_range(start = 2, end = 4, event = 'continuous event')  
        # scr.print_script()      


    def test_all_generations_and_times(self):
        scr = slime.RecentHistory(final_gen=5, chrom_length = 10)
        scr.add_event_over_time_range(start = 2, end = 4, event = 'continuous event') 
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '2:4', '5 late'])

    def test_all_early_events(self):
        scr = self.scr
        self.assertEqual(scr.all_early_events(), 
            ['1 early', '4 early'])

    def test_all_late_events(self):
        scr = self.scr
        self.assertEqual(scr.all_late_events(), 
            ['2 late', '3 late', '5 late'])

    def test_all_early_or_late_events(self):
        scr = self.scr
        self.assertEqual(scr.all_early_or_late_events(),
            ['1 early', '2 late', '3 late', '4 early', '5 late'])

    def test_all_continuous_processes(self):
        scr = slime.RecentHistory(final_gen = 10, chrom_length = 10)
        scr.add_event((5, 'late'), 'event in gen 5')
        scr.add_continuous_process((3, 10), event = 'event')
        self.assertEqual(scr.all_continuous_processes(), [(3,4), (6,9)])

    def test_all_events_at_a_given_time(self):
        scr = self.scr
        scr.all_events_at_a_given_time('3 late')
        self.assertEqual(scr.all_events_at_a_given_time('3 late'),
            ['first event', 'second event', 'third event'])

    def test_break_up_range(self):
        scr = slime.RecentHistory(final_gen=5, chrom_length = 10)
        scr.add_continuous_process((2, 4), event = 'growth')
        scr.add_event((2, 'early'), event = 'event')
        events = scr.all_generations_and_times()
        self.assertEqual(events, ['1 early', '2 early', '3:4', '5 late'])
        scr = slime.RecentHistory(final_gen=5, chrom_length = 10)
        scr.add_continuous_process((2, 3), event = 'growth')
        scr.add_event((2, 'early'), event = 'event')
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '2 early', '3 early', '5 late'])
        scr = slime.RecentHistory(final_gen=5, chrom_length = 10)
        scr.add_continuous_process((2, 4), event = 'growth')
        scr.add_event((4, 'early'), event = 'event')
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '2:3', '4 early', '5 late'])
        scr = slime.RecentHistory(final_gen=5, chrom_length = 10)
        scr.add_continuous_process((3, 4), event = 'growth')
        scr.add_event((4, 'early'), event = 'event')
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '3 early', '4 early', '5 late'])
        scr = slime.RecentHistory(final_gen = 10, chrom_length = 10)
        scr.add_continuous_process((2, 8), event = 'continuous')
        scr.add_event((5, 'early'), event = 'middle')
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '2:4', '5 early', '6:8', '10 late'])
        scr = slime.RecentHistory(final_gen = 10, chrom_length = 10)
        scr.add_continuous_process((2, 8), event = 'continuous')
        scr.add_event((3, 'early'), event = 'middle')
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '2 early', '3 early', '4:8', '10 late'])
        scr = slime.RecentHistory(final_gen = 10, chrom_length = 10)
        scr.add_continuous_process((2, 8), event = 'continuous')
        scr.add_event((8, 'early'), event = 'middle')
        self.assertEqual(scr.all_generations_and_times(),
            ['1 early', '2:7', '8 early', '10 late'])
#     # def test_errors(self):
#     #     scr = slime.RecentHistory(final_gen = 5)
#     #     self.assertRaises(ValueError, scr.add_event, (5, 'late', 'an-event'))


class TestDemographyConfig(unittest.TestCase):


    # ref0_config = msprime.PopulationConfiguration(0, 100, growth_rate = 0)
    # ref1_config = msprime.PopulationConfiguration(0, 100, growth_rate = 0)
    # adm_config = msprime.PopulationConfiguration(0, 100, growth_rate = .04)
    # scr = slime.RecentHistory(final_gen = 10, chrom_length = 10)
    # scr.add_reference_population(ref0_config, 'ref0')
    # scr.add_reference_population(ref1_config, 'ref1')
    # scr.add_admixed_population(adm_config, 'adm', proportion,
    #     single_pulse = True)


    # def test_slim_run(self):
    #     scr = slime.RecentHistory(final_gen = 5)        
    #     scr.run_slim()

    def test_demographic_input_type(self):
        scr = slime.RecentHistory(final_gen = 5, chrom_length = 10)
        config = msprime.PopulationConfiguration(0, 100, growth_rate = 0)
        scr.add_reference_population(config, 'pop1')
        # scr.print_script()

    def test_add_reference_population(self):
        scr = slime.RecentHistory(final_gen = 5, chrom_length = 10)
        config = msprime.PopulationConfiguration(0, 100, growth_rate = 0)
        scr.add_reference_population(config, 'ref0') 
        # scr.print_script()

    def test_add_admixed_population(self):
        scr = slime.RecentHistory(final_gen = 5, chrom_length = 50000000)
        rates = "examples/recomb_map_chromlength10.txt"
        scr.initialize_recombination(rate = rates, constant = False)
        ref0_config = msprime.PopulationConfiguration(0, 100, growth_rate = 0)
        ref1_config = msprime.PopulationConfiguration(0, 100, growth_rate = 0)
        config = msprime.PopulationConfiguration(0, 100, growth_rate = .04)
        scr.add_reference_population(ref0_config, 'ref0')
        scr.add_reference_population(ref1_config, 'ref1')
        scr.add_admixed_population(popConfig = config, popLabel = 'adm', proportions = [0.3,0.7], single_pulse = False, migration_rate = .02) 
        scr.print_script() 
        # scr.run_slim()

    def test_initialize_recombination(self):
        # constant recombination.
        scr = slime.RecentHistory(final_gen = 5, chrom_length = 10)
        rates = "examples/recomb_map_chromlength10.txt"
        scr.initialize_recombination(rate = rates, constant = False)
        # scr.print_script()
        scr.initialize_recombination(0.2)
        # scr.print_script()     

    def test_add_demographic_events(self):
        scr = slime.RecentHistory(final_gen = 10, chrom_length = 10)
        config = msprime.PopulationConfiguration(0, 100, growth_rate = 0)
        config1 = msprime.PopulationConfiguration(0, 100, growth_rate = 0.05)
        scr.add_reference_population(config, 'pop0')
        scr.add_reference_population(config, 'pop1')
        scr.add_admixed_population(config1, 'adm', proportions = [0.2, 0.8])
        change1 = msprime.PopulationParametersChange(3, growth_rate =  .5, population_id = 0)
        change2 = msprime.PopulationParametersChange(6, growth_rate =  .6, population_id = 0)
        change3 = msprime.PopulationParametersChange(2, growth_rate =  .3, population_id = 1)
        change4 = msprime.PopulationParametersChange(8, growth_rate =  0, population_id = 1, initial_size = 100)
        scr.add_demographic_events([change3, change1, change2, change4])
        scr.print_script()

    def test_delete_event(self):
        scr = slime.RecentHistory(final_gen = 10, chrom_length = 10)
        scr.delete_event("simulationFinished")
        self.assertEqual(len(scr.all_events_at_a_given_time('10 late')), 1)


    # def test_add_single_pulse_admixture(self):
    #     scr.
