import os
import re
import numpy as np
import pandas as pd
import msprime, pyslim
import tskit # for now, this is my local version of tskit.
import subprocess

class RecentHistory(object):
    """Creates a SLiM script of recent history that the user
    can feed into SLiMe."""
    def __init__(self, final_gen, chrom_length, reference_configs, adm_configs, prop,
        recombination=0, ref_labels=None, adm_label=None, reference_populations=None,
        admixed_population=None, mutations=None, outfile='recent-history.trees',
        model_type='WF', scriptfile='recent-history.slim'):
        self.outfile = outfile
        self.scriptfile = scriptfile
        self.model_type = model_type
        self.final_gen = int(final_gen)
        assert final_gen > 0
        self.chrom_length = chrom_length
        assert chrom_length > 0
        self.sample_sizes = []
        self.initial_sizes = []
        self.growth_rates = []
        self.population_labels = []
        self.populations = [i for i in range(0,len(self.population_labels))]
        self.initial_growth_rates = [0 for i in self.population_labels]

        # Initialize the script.
        self.script = """
initialize(){
    initializeSLiMModelType("%s");
    initializeTreeSeq();
}    

1 early(){      
}""" % self.model_type
        self.reference_populations = reference_populations
        command_to_save_out = "sim.treeSeqOutput(\"%s\")" % self.outfile
        self.add_event((final_gen, "late"), command_to_save_out)
        self.add_event((final_gen, "late"), "sim.simulationFinished()")

        # Add population information into the script.
        if ref_labels is None:
            self.ref_labels = ['pop_ref%i' % i for i in range(0, len(reference_configs))]
        else:
            assert len(ref_labels) == len(reference_configs)
            self.ref_labels = ref_labels
        assert len(reference_configs) > 1
        for pop, label in zip(reference_configs, self.ref_labels):
            assert isinstance(pop, msprime.PopulationConfiguration)
            self.add_reference_population(pop, label)

        if adm_label is None:
            self.adm_label = 'pop_adm'
        else:
            assert len(adm_label) == 1
            self.adm_label = adm_label
        for p in prop: assert p >=0 and p <= 1
        assert sum(prop) == 1
        assert len(reference_configs) == len(prop)
        assert isinstance(adm_configs, msprime.PopulationConfiguration)
        self.add_admixed_population(adm_configs, adm_label, proportions=prop)

        # Add mutations.
        # We need to initialize a genomic element with a mutation type,
        # even if the mutation rate is 0.
        if mutations is None:
            self.mutations = MutationTypes(0, [0], [1.0])
        else:
            assert isinstance(mutations, MutationTypes)
            self.mutations = mutations
        MutationRate, MutationType, GenomicElementType = self.mutations.output()
        self.initialize("initializeMutationRate(%s)" % MutationRate)
        for m in MutationType:
            self.initialize("initializeMutationType(%s)" % m)
        self.initialize("initializeGenomicElementType(%s)" % GenomicElementType)
        # For the moment, assume just one 'genomic element' g1 spanning the whole chromosome.
        self.initialize('initializeGenomicElement(g1, 0, %i)' % self.chrom_length)

        # Add recombination.
        rate = recombination
        if isinstance(rate, (int, float)):
            assert isinstance(recombination, (float, int))
            assert rate >= 0
            assert rate <= 1
            self.initialize('initializeRecombinationRate(%f)' % rate)
        else:
            assert os.path.isfile(rate)
            self.initialize("lines = readFile('%s')" % rate)
            self.initialize("lines = lines[1:(size(lines)-1)]")
            self.initialize("rates = NULL")
            self.initialize("ends = NULL")
            self.initialize("""for (line in lines)
    {
        components = strsplit(line, " ");
        ends = c(ends, asInteger(components[1]));
        rates = c(rates, asFloat(components[2]));
    }
    ends = ends - 1""")
            self.initialize("rates = rates * 1e-8")
            self.initialize("initializeRecombinationRate(rates, ends)")       

    def initialize(self, event, start = False):
        event_ind = self.find_event_index(start = start, initialization = True)
        new_script = self.script[:event_ind] + """
    %s""" % event + ";" + self.script[event_ind:]
        self.script = new_script

    def dump_script(self):
        return(self.script)

    def save_script(self, file = None):
        if file is None:
            file = self.scriptfile
        slim_file = open(file, "w")
        slim_file.writelines(self.script)
        slim_file.close()

    def print_script(self):
        print(self.script)

    def run_slim(self, verbose=True):
        self.save_script(file=self.scriptfile)
        if verbose:
            subprocess.check_call(["slim", "-l", self.scriptfile])
        else:
            subprocess.check_output(["slim", self.scriptfile])

    def time_is_continuous(self, time):
        if len(time) == 2:
            # print(time[0])
            if not isinstance(time[0], (int, np.integer)):
                raise TypeError('Generation number must be an integer.')
            if isinstance(time[1], (int, np.integer)):
                return(True)
            elif time[1] == 'early' or time[1] == 'late':
                return(False)
            else:
                raise ValueError('Time is not in correct format.')
        elif isinstance(time, (int, np.integer)):
            return(False)
        else:
            raise ValueError('Time is not in correct format.')

        
    def time_already_in_script(self, time):
        if not self.time_is_continuous(time):
            gen_time = "%i %s()" % (time[0], time[1])
            return(gen_time in self.script)
        else:
            pass
    
    def find_event_index(self, time = None, start = False, 
        initialization = False):
        """Finds the index of the script at which a new event should
        be inserted."""
        if initialization:
            gen_time = "initialize(){"
        else:
            assert time is not None
            if isinstance(time[1], str):
                gen_time = "%i %s(){" % (time[0], time[1])
            else:
                gen_time = "%i:%i {" % (time[0], time[1])
        gen_location = self.script.find("%s" % gen_time)
        start_loc = gen_location + len(gen_time)
        if start:
            return(start_loc)
        else:
            rest_of_script = self.script[start_loc:]
            end_loc = rest_of_script.find('\n}')
            return(start_loc + end_loc)
               
    def find_time_index(self, time):
        """Finds the index of the script at which a new generation and
        time should be inserted."""
        generation = time[0]
        times = self.all_generations_and_times()
        all_early_or_late_events = self.all_early_or_late_events()

        gen_regex =re.compile(r"\d+")
        if not self.time_is_continuous(time):
            time_regex=re.compile(r"\bearly\b|\blate\b")
            BREAK_TRIGGERED = 0
            for m in times:
                current_gen = int(gen_regex.search(m).group())
                if current_gen > generation:
                    BREAK_TRIGGERED = 1
                    break
                elif current_gen == generation and "late" in m and time[1] == "early":
                    BREAK_TRIGGERED = 1
                    break
            if BREAK_TRIGGERED:
                gen_location = self.script.find("%s" % m)
            else:
                gen_location = len(self.script)
            return(gen_location, BREAK_TRIGGERED)
        else:
            return(len(self.script), 1)

    def add_event(self, time, event, insert_at_start=False, check_continuous_processes = True):
        # check whether there is already an event spanning a range of generations
        # that covers the time of the new event. If so, this range must be broken up.
        generation = time[0]
        if check_continuous_processes:
            for start, end in self.all_continuous_processes():
                event_range = "%i:%i {" % (start, end)
                existing_events = self.all_events_at_a_given_time(event_range[:-2])
                regex = re.compile(event_range)
                BREAK_TRIGGERED = 0
                if generation == start:
                    if start + 1 < end:
                        new_script = regex.sub("%i:%i {" % (start + 1, end), self.script)
                    else:
                        new_script = regex.sub("%i early(){" % end, self.script)
                    self.script = new_script
                    BREAK_TRIGGERED = 1
                elif start < generation and generation < end:
                    if generation + 1 == end:
                        for e in existing_events:
                            self.add_event((end, 'early'), e)
                    else:
                        for e in existing_events:
                            self.add_continuous_process((generation + 1, end), e) 

                    if start + 1 == generation:
                        time_to_put_in = "%i early(){" % start
                    else:
                        time_to_put_in = "%i:%i {" % (start, generation - 1)
                    new_script = regex.sub("%s" % time_to_put_in, self.script)
                    self.script = new_script
                    BREAK_TRIGGERED = 1           
                elif generation == end:
                    if start + 1 < end:
                        new_script = regex.sub("%i:%i {" % (start, end - 1), self.script)
                    else:
                        new_script = regex.sub("%i early(){" % start, self.script)
                    self.script = new_script
                    BREAK_TRIGGERED = 1

                if BREAK_TRIGGERED:
                    for e in existing_events:
                        self.add_event((time[0], time[1]), e, check_continuous_processes = False)  
                    break              

        # Check for existing event at that time.
        if not self.time_already_in_script(time):
            time_ind, NOT_AT_END = self.find_time_index(time)
            if NOT_AT_END:
                new_script = self.script[:time_ind] + """%i %s(){
}
""" % (time[0], time[1]) + "\n" + self.script[time_ind:]
            else:
                new_script = self.script + "\n" + """
%i %s(){
}""" % (time[0], time[1])
            self.script = new_script
        event_ind = self.find_event_index((time[0], time[1]), insert_at_start)
        new_script = self.script[:event_ind] + """
    %s""" % event + ";" + self.script[event_ind:]
        self.script = new_script

    def add_event_over_time_range(self, start, end, event, insert_at_start = False):
        time_ind, NOT_AT_END = self.find_time_index((start, 'early'))
        if not NOT_AT_END:
            raise ValueError("End generation is after the end of the simulation.")
        # See if range is already in script.
        regex = re.compile("%i:%i {" % (start, end))
        ALREADY_IN_SCRIPT = regex.findall(self.script)
        if ALREADY_IN_SCRIPT:
            event_ind = self.find_event_index((start, end), insert_at_start)
            new_script = self.script[:event_ind] + """
    %s""" % event + ";" + self.script[event_ind:]
            self.script = new_script
        else:
            new_script = self.script[:time_ind] + """%i:%i {
    %s;
}
""" % (start, end, event) + "\n" + self.script[time_ind:]
            self.script = new_script 

    def add_continuous_process(self, time, event, early=True):
        """
        Adds an event to each generation in the specified range.
        If there's already a one-off event at a generation within the specified range,
        the event will also be added as a one-off event in that generation. (By default,
        this will happen at the 'early' time period - but this behaviour can be controlled
        by the `early` parameter.).
        """
        start = time[0]
        end = time[1]
        missing_gens = []
        gen_range = [i for i in np.arange(start, end + 1)]
        for gen in gen_range:
            if self.time_already_in_script((gen, 'early')) or self.time_already_in_script((gen, 'late')):
                if early:
                    self.add_event((gen, 'early'), event)
                else:
                    self.add_event((gen, 'late'), event)
            else:
                missing_gens.append(gen)

        # Process the remaining generations. See helpful answer at
        # https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
        if len(missing_gens) == 1:
            gen = missing_gens[0]
            if early:
                self.add_event((gen, 'early'), event)
            else:
                self.add_event((gen, 'late'), event)

        else:
            def gen_groups():
                first = last = missing_gens[0]
                for gen in missing_gens[1:]:
                    if gen == last + 1:
                        last = gen
                    else:
                        yield first, last
                        first = last = gen
                yield first, last
            missing_ranges = list(gen_groups())
            for start, end in missing_ranges:
                time = (start, end)
                self.add_event_over_time_range(time[0], time[1], event)

    def gen_in_range(self, gen, start, end):
        return(start <= gen and gen <= end)

    def all_generations_and_times(self):
        regex = re.compile(r"\d+ \bearly\b|\d+ \blate\b|\d+:\d")
        return(regex.findall(self.script))

    def all_early_events(self):
        regex = re.compile(r"\d+ \bearly\b")
        return(regex.findall(self.script))

    def all_late_events(self):
        regex = re.compile(r"\d+ \blate\b")
        return(regex.findall(self.script))

    def all_early_or_late_events(self):
        regex = re.compile(r"\d+ \bearly\b|\d+ \blate\b")
        return(regex.findall(self.script))

    def all_continuous_processes(self):
        regex = re.compile("\d+:\d")
        ranges = regex.findall(self.script)
        to_return = []
        for range in ranges:
            gens = re.compile("\d")
            pair = gens.findall(range)
            to_return.append((int(pair[0]), int(pair[1])))
        return(to_return)

    def all_events_at_a_given_time(self, time):
        assert time in self.all_generations_and_times()
        ind = self.script.find(time)
        script_shortened = self.script[ind:]
        start_ind = script_shortened.find('{')
        end_ind = script_shortened.find('}')
        all_events_str = script_shortened[start_ind + 1:end_ind]
        all_events_list = all_events_str.split(";")
        all_events_list = [i.strip() for i in all_events_list]
        all_events_list = [i for i in all_events_list if len(i) > 0]
        return(all_events_list)

    def add_reference_population(self, popConfig, popLabel):
        if not isinstance(popConfig, msprime.PopulationConfiguration):
            TypeError("popConfig must be a msprime.PopulationConfiguration object.")
        sample_size = popConfig.sample_size
        initial_size = popConfig.initial_size
        growth_rate = np.exp(popConfig.growth_rate)
        self.sample_sizes.append(sample_size)
        self.initial_sizes.append(initial_size)
        self.growth_rates.append(growth_rate)
        self.population_labels.append(popLabel)
        ind = len(self.population_labels) - 1
        self.populations.append(ind)
        self.add_event((1, 'early'), "sim.addSubpop(\"p%i\", %i)" % (ind, initial_size))
        # if growth_rate != 1:
        ind = self.population_labels.index(popLabel)
        self.initial_growth_rates.append(popConfig.growth_rate)

    def add_admixed_population(self, popConfig, popLabel, proportions, single_pulse = True, migration_rate = None):
        if not isinstance(popConfig, msprime.PopulationConfiguration):
            raise TypeError("popConfig must be a msprime.PopulationConfiguration object.")
        sample_size = popConfig.sample_size
        initial_size = popConfig.initial_size
        growth_rate = np.exp(popConfig.growth_rate)
        self.sample_sizes.append(sample_size)
        self.initial_sizes.append(initial_size)
        self.growth_rates.append(growth_rate)
        self.population_labels.append(popLabel)
        ind = len(self.population_labels) - 1
        self.populations.append(ind)
        self.add_event((1, 'late'), "sim.addSubpop(\"p%i\", %i)" % (ind, initial_size))
        ind = self.population_labels.index(popLabel)
        self.initial_growth_rates.append(popConfig.growth_rate)
            # self.add_continuous_process((2, self.final_gen),
            #     event = "p%i.setSubpopulationSize(asInteger(p%i.individualCount * %f))" % (ind, ind, growth_rate))
        # Add admixture in.
        if not len(proportions) == len(self.population_labels) - 1:
            raise SystemError('A proportion must be allocated to each reference population.')
        if np.sum(proportions) != 1:
            raise ValueError('Admixture proportions must sum to 1.')
        for p in proportions:
            assert p >= 0
            assert p <= 1
        pop_labels = ["p%i" % i for i in range(0, len(self.population_labels) - 1)]
        command_pops = list_to_slim_vector(pop_labels)
        prop_as_string = ["%f" % i for i in proportions]
        command_prop = list_to_slim_vector(prop_as_string)
        command = "p%i.setMigrationRates(%s,%s)" % (len(self.population_labels) - 1, command_pops, command_prop)
        self.add_event((1, 'late'), "%s" % command)
        if single_pulse:
            zeros_as_string = ["0.0" for i in range(0, len(self.population_labels) - 1)]
            final_zeros = list_to_slim_vector(zeros_as_string)
            end_command = "p%i.setMigrationRates(%s,%s)" % (len(self.population_labels) - 1, command_pops, final_zeros)
            self.add_event((2, 'late'), end_command)
        else:
            assert migration_rate >= 0
            assert migration_rate <= 1
            pop_labels = ["p%i" % i for i in range(0, len(self.population_labels))]
            command_pops = list_to_slim_vector(pop_labels)
            scaled_props = ["%f" % (p * migration_rate) for p in proportions]
            scaled_props.append("%f" % (1 - migration_rate))
            command_prop = list_to_slim_vector(scaled_props)
            end_command = "p%i.setMigrationRates(%s,%s)" % (len(self.population_labels) - 1, command_pops, command_prop)
            self.add_event((2, 'late'), end_command)

    def add_demographic_events(self, event_list):
        # Sort events by type.
        param_changes = []
        migr_rate_changes = []
        mass_migrations = []
        for e in event_list:
            if isinstance(e, msprime.PopulationParametersChange):
                param_changes.append(e)
            elif isinstance(e, msprime.MigrationRateChange):
                migr_rate_changes.append(e)
            elif isinstance(e, msprime.MassMigration):
                mass_migrations.append(e)
            else:
                raise TypeError("All items in the event_list must be a PopulationParametersChange, MigrationRateChange or MassMigration object.")

        # Add PopulationParametersChange events.
        self.add_population_parameters_change(param_changes)

    def delete_event(self, string, time = None):
        # Deletes first line of the script containing the given string.
        if time is None:
            event_loc = self.script.find(string)
            script_start = self.script[:event_loc]
            script_start = script_start[:script_start.rfind("\n")]
            script_end = self.script[event_loc:]
            script_end = script_end[script_end.find("\n"):]
            self.script = script_start + script_end
        else:
            time_str = time_as_string(time)
            assert time_str in self.all_generations_and_times()
            assert any(string in s for s in self.all_events_at_a_given_time(time_str))
            time_loc = self.script.find(time_str)
            print(time_loc, "time_loc")
            script_start = self.script[:time_loc]
            script_end = self.script[time_loc:]
            event_loc = script_end.find(string)
            event_start = script_end[:event_loc]
            event_start = event_start[:event_start.rfind("\n")]
            event_end = script_end[event_loc:]
            event_end = event_end[event_end.find("\n"):]
            script_end = event_start + event_end
            self.script = script_start + script_end


    def add_population_parameters_change(self, paramChanges):
        """
        Takes a list of msprime.PopulationParametersChange objects
        and 'translates' them into SLiM commands.
        https://msprime.readthedocs.io/en/stable/api.html#demographic-events
        There are two types of changes: exponential growth rate changes and
        instantaneous population size changes.
        """
        param_changes = paramChanges
        growth_rates = self.initial_growth_rates
        events_to_add = []
        event_times = []
        if any(self.initial_growth_rates) != 0:
            event_times.append(2)
        prev_time = 1
        for e in param_changes:
            if e.time < prev_time:
                raise RuntimeError("At the moment, events must be ordered by forwards-time.")
            if e.growth_rate is not None:
                event_times.append(e.time)
        event_times = list(set(event_times)) # Remove duplicate times.
        event_times.sort()
        # Initialize EventList objects.
        if len(event_times) == 1:
            events_to_add.append(EventList(start_time = event_times[0], 
                end_time = self.final_gen))
        else:
            for first, last in zip(event_times[:-1], event_times[1:]):
                events_to_add.append(EventList(start_time = first, end_time = last))
            events_to_add.append(EventList(start_time = last, end_time = self.final_gen))
        events_to_add = iter(events_to_add)
        current_event = events_to_add.__next__()
        for p in param_changes:
            current_time = p.time
            # 1. Process growth rates.
            pop = p.population
            if p.growth_rate is not None:
                growth_rates[pop] = p.growth_rate           
            # Once there are no more events in an epoch to be processed,
            # add to the relevant event.
            # Note. This assumes that event_list is ordered by time!!!
            if prev_time != current_time:
                while current_event.start_time < current_time:
                    current_event = events_to_add.__next__()
                for pop in self.populations:
                    if growth_rates[pop] != 0:
                        current_event.add_event("p%i.setSubpopulationSize(asInteger(p%i.individualCount * exp(%f)))" % (pop, pop, growth_rates[pop]))
                # Add events into the script.
                for e in current_event.events:
                    if current_event.start_time + 1 == current_event.end_time:
                        self.add_event((current_event.start_time, 'early'), e)
                    else:
                        self.add_event_over_time_range(current_event.start_time, current_event.end_time - 1, e)

            # 2. Process instantaneous population size changes. 
            if p.initial_size is not None:
                command = "p%i.setSubpopulationSize(%i)" % (pop, p.initial_size)
                self.add_event((current_time, 'early'), command)
                # If there is continuous growth at this time, remove it for this generation only.
                popchange = "p%i.setSubpopulationSize(" % pop
                events_at_time = self.all_events_at_a_given_time("%i early" % current_time)
                no_of_popchanges = sum(popchange in e for e in events_at_time)
                if no_of_popchanges > 1:
                    self.delete_event(popchange, (current_time, 'early'))

def time_as_string(time):
    assert len(time) == 2
    assert isinstance(time[0], (int, np.integer))
    if isinstance(time[1], str):
        assert time[1] == 'early' or time[1] == 'late'
        gen_time = "%i %s" % (time[0], time[1])
    else:
        assert isinstance(time[1], (int, np.integer))
        gen_time = "%i:%i" % (time[0], time[1])
    return(gen_time)


class EventList(object):
    """
    Holds a start and and end time, and a list of events that happens
    between those times.
    """
    def __init__(self, start_time, end_time = None):
        self.start_time = start_time
        if end_time is None:
            self.end_time = self.start_time
        else:
            assert end_time >= start_time
            self.end_time = end_time
        self.events = []

    def __lt__(self, other):
        return (self.start_time, self.end_time) < (other.start_time, other.end_time)

    def add_event(self, event):
        assert isinstance(event, str)
        self.events.append(event)

def list_to_slim_vector(list_of_strings):
    # Changes ['s1', 's2', 's3'] into 'c(s1,s2,s3)'
    return("c(" + ",".join(list_of_strings) + ")")
            
class MutationTypes(object):
    """
    Holds information about mutations to put into the simulation.
    """
    def __init__(self, mutation_rate, selection_coeffs, proportions, dominance_coeffs = None):
        self.mutation_rate = mutation_rate
        self.selection_coeffs = selection_coeffs
        assert len(proportions) == len(selection_coeffs)
        assert np.sum(proportions) == 1
        self.proportions = proportions
        if dominance_coeffs is None:
            self.dominance_coeffs = [0.5 for i in range(0, len(selection_coeffs))]
        else:
            self.dominance_coeffs = dominance_coeffs
        
    def output(self):
        MutationRate = self.mutation_rate
        MutationType = []
        for m in range(0, len(self.selection_coeffs)):
            type_to_add = "\"m%i\", %f, \"f\", %f" % (m, self.dominance_coeffs[m], self.selection_coeffs[m])
            MutationType.append(type_to_add)
        mut_labels = "c(" + ",".join(["m%i" % m for m in range(0, len(self.selection_coeffs))]) + ")"
        proportions_str = [str(p) for p in self.proportions]
        proportions = "c(" + ",".join(proportions_str) + ")"
        GenomicElementType = "\"g1\", %s, %s" % (mut_labels, proportions)
        return(MutationRate, MutationType, GenomicElementType)
        
    def print_output(self):
        # Prints the commands to be inserted into the SLiM initialization.
        MutationRate, MutationType, GenomicElementType = self.output()
        print("initializeMutationRate(%s);" % MutationRate)
        for m in MutationType:
            print("initializeMutationType(%s);" % m)
        print("initializeGenomicElementType(%s);" % GenomicElementType)


