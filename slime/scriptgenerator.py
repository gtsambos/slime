import os
import re
import numpy as np
import pandas as pd
import msprime, pyslim
import tskit # for now, this is my local version of tskit.
import subprocess
import itertools

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
        self.initial_prop = prop
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
        ends = c(ends, asInteger(components[0]));
        rates = c(rates, asFloat(components[1]));
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
            gen_time = "\n%i %s()" % (time[0], time[1])
            return(gen_time in self.script)
        else:
            pass # FILL THIS OUT
    
    def find_event_index(self, time = None, start = False, 
        initialization = False):
        """Finds the index of the script at which a new event should
        be inserted."""
        if initialization:
            gen_time = "initialize(){"
        else:
            assert time is not None
            if isinstance(time[1], str):
                gen_time = "\n%i %s(){" % (time[0], time[1])
            else:
                gen_time = "\n%i:%i {" % (time[0], time[1])
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

    def add_event(self, time, event, insert_at_start=False, check_continuous_processes = True,
        early=True):

        # Loop through existing times.
        if early:
            early_or_late = 'early'
        else:
            early_or_late = 'late'
        input_list = self.all_generations_and_times(return_as_strings=False)
        EVENT_ADDED = False

        for ind in range(0, len(input_list)):
            t = input_list[ind]

            # If there is an existing event at this time.
            if t == time:
                self.add_event_to_single_time(time, event, insert_at_start)
                EVENT_ADDED = True
                
            # If there is a continuous range spanning this time.
            elif isinstance(t[1], int) and t[0] <= time[0] and time[0] <= t[1]:

                # Remove existing range.
                existing_events = self.all_events_at_a_given_time(time = "%i:%i" % (t[0], t[1]))
                self.delete_time(time="%i:%i " % (t[0], t[1]))

                if t[0] == time[0] - 1:
                    time_to_add = (t[0], early_or_late)
                    self.add_new_single_time(time=time_to_add)
                    for e in existing_events:
                        self.add_event_to_single_time(time_to_add, e, insert_at_start)
                elif t[0] < time[0] - 1:
                    for e in existing_events:
                        self.add_event_over_time_range(t[0], time[0] - 1, e)
                            
                self.add_new_single_time((time[0], early_or_late))
                if time[1] == 'late':
                    self.add_new_single_time(time)
                for e in existing_events:
                    self.add_event_to_single_time((time[0], early_or_late), e, insert_at_start)
                self.add_event_to_single_time(time, event, insert_at_start)
                            
                if t[1] == time[0] + 1:
                    time_to_add = (t[1], early_or_late)
                    self.add_new_single_time(time=time_to_add)
                    for e in existing_events:
                        self.add_event_to_single_time(time_to_add, e, insert_at_start)
                        
                elif t[1] > time[0] + 1:
                    for e in existing_events:
                        self.add_event_over_time_range(time[0]+1, t[1], e)

                EVENT_ADDED = True

        # If there's no existing event at this time, add one.
        if not EVENT_ADDED:
            self.add_new_single_time(time)
            self.add_event_to_single_time(time, event, insert_at_start)

    def add_event_to_single_time(self, time, event, insert_at_start):
        event_ind = self.find_event_index((time[0], time[1]), insert_at_start)
        new_script = self.script[:event_ind] + """
    %s""" % event + ";" + self.script[event_ind:]
        self.script = new_script        

    def add_new_single_time(self, time):
        """
        Creates a new code block for the supplied single time.
        """
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

    def add_continuous_process(self, new_range, event, early=True):
        """
        Adds an event that occurs over an interval of generations.
        Breaks up existing events if necessary.
        """
        input_list = self.all_generations_and_times(return_as_strings=False)
        if early:
            early_or_late = 'early'
        else:
            early_or_late = 'late'
        new_times = []
        # Loop through times, and print each 'new thing' to be done.
        for ind in range(0, len(input_list)):
            time = input_list[ind]
            
            # Single event ranges
            if isinstance(time[1], str):
                # If there is both an early and late event in this generation,
                # add the event only once.
                if time[0] == input_list[ind - 1][0]:
                    pass
                elif new_range[0] <= time[0] and time[0] <= new_range[1]:
                    self.add_event(time=(time[0], early_or_late), event=event)
                    new_times.append(time[0])
                    
            # Continuous event ranges
            else:
                if time[1] >= new_range[0] and new_range[1] >= time[0]:

                    existing_events = self.all_events_at_a_given_time(time = "%i:%i" % (time[0], time[1]))
                    self.delete_time(time="%i:%i " % (time[0], time[1]))

                    if time[0] == new_range[0] - 1:
                        for e in existing_events:
                            self.add_event(time=(time[0], early_or_late), event=e)
                        new_times.append(time[0])
                    elif time[0] < new_range[0] - 1:
                        for e in existing_events:
                            self.add_event_over_time_range(time[0], new_range[0] - 1, e)
                        for t in range(time[0], new_range[0]):
                            new_times.append(t)

                    if max(time[0], new_range[0]) == min(time[1], new_range[1]):
                        t = (max(time[0], new_range[0]), early_or_late)
                        for e in existing_events:
                            self.add_event(time=t, event=e)
                    else:
                        t = (max(time[0], new_range[0]), min(time[1], new_range[1]))
                        for e in existing_events:
                            self.add_event_over_time_range(t[0], t[1], event=e)
                        self.add_event_over_time_range(t[0], t[1], event=event)

                        for t in range(max(time[0], new_range[0]), min(time[1], new_range[1]) + 1):
                            new_times.append(t)
                    
                    if time[1] == new_range[1] + 1:
                        for e in existing_events:
                            self.add_event(time=(time[1], early_or_late), event=e)
                        new_times.append(time[1])
                        
                    elif time[1] > new_range[1] + 1:
                        for e in existing_events:
                            self.add_event_over_time_range(new_range[1] + 1, time[1], event=e)
                        for t in range(new_range[1] + 1, time[1] + 1):
                            new_times.append(t)
                            
        remaining_times =  set([i for i in range(new_range[0], new_range[1]+1)]) - set(new_times)
        remaining_times = list(remaining_times)
        
        def ranges(i):
            for a, b in itertools.groupby(enumerate(i), lambda args: args[1] - args[0]):
                b = list(b)
                yield b[0][1], b[-1][1]
                
        for r in ranges(remaining_times):
            # single events
            if r[0] == r[1]:
                self.add_event(time=(r[0], early_or_late), event=event)
            else:
                self.add_event_over_time_range(r[0], r[1], event)

    def gen_in_range(self, gen, start, end):
        return(start <= gen and gen <= end)

    def all_generations_and_times(self, return_as_strings=True):
        regex = re.compile(r"\d+ \bearly\b|\d+ \blate\b|\d+:\d+")
        ret = regex.findall(self.script)
        if return_as_strings:
            return(ret)
        else:
            ret_t = []
            for r in ret:
                if ":" in r:
                    gens = re.compile("\d+")
                    pair = gens.findall(r)
                    ret_t.append((int(pair[0]), int(pair[1])))
                else:
                    pair = r.split(" ")
                    ret_t.append((int(pair[0]), str(pair[1])))
            return(ret_t)

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
        regex = re.compile("\d+:\d+")
        ranges = regex.findall(self.script)
        to_return = []
        for range in ranges:
            gens = re.compile("\d+")
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

    def delete_time(self, time):
        time_0 = "\n%s" % time
        start_ind = self.script.find(time_0)
        script_shortened = self.script[start_ind:]
        end_ind = script_shortened.find('}\n')
        new_script = self.script[:start_ind] + script_shortened[end_ind+2:]
        self.script = new_script

    def add_reference_population(self, popConfig, popLabel):
        if not isinstance(popConfig, msprime.PopulationConfiguration):
            TypeError("popConfig must be a msprime.PopulationConfiguration object.")
        sample_size = popConfig.sample_size
        initial_size = popConfig.initial_size
        growth_rate = popConfig.growth_rate
        self.sample_sizes.append(sample_size)
        self.initial_sizes.append(initial_size)
        self.growth_rates.append(growth_rate)
        self.population_labels.append(popLabel)
        ind = len(self.population_labels) - 1
        self.populations.append(ind)
        self.add_event((1, 'early'), "sim.addSubpop(\"p%i\", %i)" % (ind, initial_size))
        # msprime.PopulationParametersChange(3, growth_rate =  .5, population_id = 0)
        if growth_rate != 0:
            self.add_continuous_process((2, self.final_gen),
                event = "p%i.setSubpopulationSize(asInteger(p%i.individualCount * exp(%f)))" % (ind, ind, growth_rate))
        ind = self.population_labels.index(popLabel)
        self.initial_growth_rates.append(popConfig.growth_rate)

    def add_admixed_population(self, popConfig, popLabel, proportions, migration_rate = None):
        if not isinstance(popConfig, msprime.PopulationConfiguration):
            raise TypeError("popConfig must be a msprime.PopulationConfiguration object.")
        sample_size = popConfig.sample_size
        initial_size = popConfig.initial_size
        growth_rate = popConfig.growth_rate
        self.sample_sizes.append(sample_size)
        self.initial_sizes.append(initial_size)
        self.growth_rates.append(growth_rate)
        self.population_labels.append(popLabel)
        ind = len(self.population_labels) - 1
        self.populations.append(ind)
        self.add_event((1, 'late'), "sim.addSubpop(\"p%i\", %i)" % (ind, initial_size))
        if growth_rate != 0:
            self.add_continuous_process((2, self.final_gen),
                event = "p%i.setSubpopulationSize(asInteger(p%i.individualCount * exp(%f)))" % (ind, ind, growth_rate))
        ind = self.population_labels.index(popLabel)
        self.initial_growth_rates.append(popConfig.growth_rate)

        # Add admixture into the script.
        self.add_mass_migration(prop=proportions, time=(1, 'late'))

    def add_migration_rate(self, rates, time=(2, 'late')):
        """
        Specifies rates of migration from references to admixed population.
        """
        assert all(r >= 0 for r in rates) and all(r <=1 for r in rates)
        assert len(rates) == len(self.populations) - 1
        ref_pops = self.populations[:-1]
        adm_pop = self.populations[-1]
        # Add migration proportions.
        ref_pop_vector = list_to_slim_vector(["p%i" % i for i in ref_pops])
        rates_vector = list_to_slim_vector(rates)
        event_to_add = "p%i.setMigrationRates(%s, %s)" % (adm_pop, ref_pop_vector, rates_vector)
        self.add_event(time=time, event=event_to_add)
        # Remove previous migrations at this time.
        time_st = time_as_string(time)
        st = "p%i.setMigrationRates" % adm_pop
        if sum(st in e for e in self.all_events_at_a_given_time(time_st)) == 2:
            self.delete_event(string=st, time=time)

    def add_mass_migration(self, prop, time):
        """
        Specifies a mass migration event from reference populations into admixed 
        population at a specified time.
        This is like changing the migration rate, but the migration rate is
        changed to 0 afterwards. 
        # Later: would be useful if it is changed back to original migration rate?
        """
        self.add_migration_rate(rates=prop, time=time)
        # Change later migration rate back to 0.
        if time[0] != self.final_gen:
            no_refs = len(self.populations) - 1
            self.add_migration_rate(rates=[0 for i in range(0, no_refs)], time=(time[0]+1, 'late'))

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
        # Use with caution...
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
            script_start = self.script[:time_loc]
            script_end = self.script[time_loc:]
            event_loc = script_end.find(string)
            event_start = script_end[:event_loc]
            event_start = event_start[:event_start.rfind("\n")]
            event_end = script_end[event_loc:]
            event_end = event_end[event_end.find("\n"):]
            script_end = event_start + event_end
            self.script = script_start + script_end

    def add_size_change(self, PopulationParametersChange, early=True):
        p = PopulationParametersChange
        if early:
            early_or_late = 'early'
        else:
            early_or_late = 'late'
        if p.population is None:
            raise ValueError("You must supply a population ID.")
        if p.initial_size is not None:
            new_time = (p.time, early_or_late)
            new_time_st = time_as_string(new_time)
            # Add new size.
            self.add_event(time=new_time,
                event="p%i.setSubpopulationSize(%i)" % (p.population, p.initial_size))
            # Delete any prior population size changes at this time.
            st = "p%i.setSubpopulationSize(" % p.population
            if sum(st in e for e in self.all_events_at_a_given_time(new_time_st)) == 2:
                self.delete_event(string=st, time=new_time)

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
    l = [str(s) for s in list_of_strings]
    # Changes ['s1', 's2', 's3'] into 'c(s1,s2,s3)'
    return("c(" + ",".join(l) + ")")
            
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


