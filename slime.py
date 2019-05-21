import os
import re
import numpy as np
import pandas as pd
import msprime, pyslim
import tskit # for now, this is my local version of tskit.

class SLiMError(Exception):
    """ Errors in the SLiM script that are likely to prevent the simulation 
    from completing.
    """
    def __init__(self, message):
        # self.expression = expression
        self.message = message

def simulate_recent_history(file, outFile = "recent-history.trees", logFile = "recent-history.log"):
    """
    Simulates genomic history from the start of admixture until
    the present day with SLiM.
    """
    command_line = "slim" + " " + file + " " + "&> " + logFile
    slim_run = os.system(command_line)
    if slim_run != 0: # Print any errors.
        log_file = open(logFile, "r")
        for line in log_file:
            print(line)

# Functions for sampling from the tree sequence.
class TreeSequenceToSample(object):
    """ Holds a tree sequence corresponding to sample of nodes from the 
    original tree sequence. """
    def __init__(self, ts, sample_nodes = None, populations_to_sample_from = None, sample_sizes = None):
        self.original_ts = ts
        self.sample_nodes = sample_nodes
        self.populations = populations_to_sample_from
        self.sample_sizes = sample_sizes
        if self.populations is not None or self.sample_sizes is not None:
            if len(self.populations) != len(self.sample_sizes):
                raise ValueError("The list of sample sizes must have the same length as the list of populations to sample from")

    def subsample(self):
        """Simplifies tree sequence for a specified subsample of individuals or nodes."""
        if self.sample_nodes is None:
            sample_nodes = self.sample_individuals()
            self.sample_nodes = sample_nodes
        return(self.original_ts.simplify(samples = sample_nodes, keep_unary = True))

    def sample_individuals(self):
        """ Returns a list of sample individual ids to be sampled from.
         Inputs: a tableCollection, two lists of integers."""

        tabs = self.original_ts.dump_tables()
        # Convert the node table into a pandas dataframe.
        data = {'time' : tabs.nodes.time, 
            'population' : tabs.nodes.population, 
            'individual' : tabs.nodes.individual}
        nodes_df = pd.DataFrame(data, columns = ['time', 
                                                 'population', 
                                                 'individual'])
        # Create a list of sample individuals corresponding to each population.
        pop_inds = []
        for pop in self.populations:
            sample_nodes = list(self.get_sample_nodes(pop))
            sample_nodes_df = nodes_df.iloc[sample_nodes]
            relevant_individuals = np.unique(sample_nodes_df[sample_nodes_df.population == pop]['individual'])
            pop_inds.append(relevant_individuals)
        # Take a subsample of individuals and flatten the list.
        sample_inds = []
        for pop in self.populations:
                index = self.populations.index(pop)
                sample_inds.append(np.random.choice(list(pop_inds[index]), 
                    size=self.sample_sizes[index], replace=False))
        sample_inds = np.concatenate(sample_inds)
        # Find the nodes corresponding to these individuals.
        subsample_nodes = []
        for node in self.original_ts.nodes():
            if node.individual in sample_inds:
                subsample_nodes.append(node.id)

        return(subsample_nodes)

    def get_sample_nodes(self, pop):
        """ Finds all present-day samples whose ancestry we wish to obtain.
         Inputs: a treeSequence, an integer. """
        sample_nodes = []
        for node in self.original_ts.nodes():
            if (node.time == 0.0) and (node.population == pop):
                sample_nodes.append(node.id)
        return(sample_nodes)

class AdmixtureSimulation(object):
    def __init__(
        self, slim_script, 
        ancient_recombination_rate,
        ancient_population_configurations,
        ancient_demographic_events,
        neutral_mutation_rate = 0,
        out_file = None, 
        populations_to_sample_from= None,
        sample_sizes = None,
        ):
        self.slim_script = slim_script
        self.out_file = out_file
        self.slim_out = None
        self.populations = populations_to_sample_from
        self.sample_sizes = sample_sizes
        if self.populations is not None:
            self.need_to_subsample = 1
        else:
            self.need_to_subsample = 0
        self.ancient_recombination_rate = ancient_recombination_rate
        self.ancient_population_configurations = ancient_population_configurations
        self.ancient_demographic_events = ancient_demographic_events
        self.neutral_mutation_rate = neutral_mutation_rate

    def go(self):
        """ A wrapper for the admixture simulation."""
        print('Simulating recent history with SLiM...')
        simulate_recent_history(self.slim_script)
        ts = tskit.load(self.slim_out)
        if self.need_to_subsample:
            print('Taking samples from present day populations...')
            ts = TreeSequenceToSample(ts, 
                populations_to_sample_from = self.populations,
                sample_sizes = self.sample_sizes)
            ts = ts.subsample()
        # tabs = ts.tables
        ts = pyslim.SlimTreeSequence.load_tables(ts.tables)
        print('Simulating ancient history with msprime...')
        ts = ts.recapitate(
            recombination_rate = self.ancient_recombination_rate,
            population_configurations = self.ancient_population_configurations,
            demographic_events = self.ancient_demographic_events,
            keep_first_generation = True # needed to get local ancestors
            )
        print('Adding variation...')
        ts = pyslim.SlimTreeSequence(msprime.mutate(ts, 
            rate=self.neutral_mutation_rate, keep=True))
        if self.out_file is not None:
            ts.dump(self.out_file)
        return(ts)

    def debugger(self):
        """ A debugger to run before the simulation. """
        # SLiM debugging
        slim = self.slim_script
        print('\nSLiM input file:', slim)
        # Test 1: is an output file saved?
        with open(slim ,'r') as f:
            lines = f.readlines()
            string_pre = ".treeSeqOutput("
            string_post = ")"
            ind = 0
            for line in lines:
                if string_pre in line and string_post in line:
                    out_file = line.split(string_pre)[1].split(string_post)[0]
                    self.slim_out = out_file.strip('""')
                    print('SLiM output file:', self.slim_out)
                    ind = 1
            if ind == 0:
                print(
    """SLiM error:
    Oh no, your script does not produce a .trees file!
    Please ensure you include a call to 'treeSeqOutput()' at the end of your script.
                    """)
        # Test 2: subsampling
        if self.populations is not None or self.sample_sizes is not None:
            if len(self.populations) != len(self.sample_sizes):
                print(
    """ Subsampling error:
    The list of populations to sample from must have the same length
    as the list of sample sizes."""
    )
            print("We are sampling:")
            for ind in range(len(self.populations)):
                print("-", self.sample_sizes[ind], "individuals from population", 
                    self.populations[ind])
        else:
            "No subsampling will be performed."
        # Test 3: demography debugging in recapitation
        print('Ancient demography:')
        dd = msprime.DemographyDebugger(
            population_configurations=self.ancient_population_configurations,
            demographic_events=self.ancient_demographic_events)
        dd.print_history()
        # Test 4: Adding variation
        # Neutral mutations
        print('Neutral mutation rate:', self.neutral_mutation_rate)

class RecentHistory(object):
    """Creates a SLiM script of recent history that the user
    can feed into SLiMe."""
    def __init__(self, final_gen, outfile='recent-history.slim', model_type='WF'):
        self.outfile = outfile
        self.model_type = model_type
        self.script = """
initialize(){
    initializeSLiMModelType("%s");
    initializeTreeSeq();
}    

1 early(){      
}""" % self.model_type
        self.final_gen = final_gen
        command_to_save_out = "sim.treeSeqOutput(\"%s\")" % self.outfile
        self.add_event((final_gen, "late"), command_to_save_out)
        self.add_event((final_gen, "late"), "sim.simulationFinished()")
        self.number_of_reference_populations = 0
        self.sample_sizes = []
        self.initial_sizes = []
        self.growth_rates = []
        self.population_labels = []

    def dump_script(self):
        return(self.script)
    
    def print_script(self):
        print(self.script)

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
    
    def find_event_index(self, time, INSERT_AT_START):
        """Finds the index of the script at which a new event should
        be inserted."""
        gen_time = "%i %s(){" % (time[0], time[1])
        gen_location = self.script.find("%s" % gen_time)
        start_loc = gen_location + len(gen_time)
        if INSERT_AT_START:
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

        # regex = re.compile(r"\d+ \bearly\b\(\)\{|\d+ \blate\b\(\)\{|\d+:\d+ \{")
        # times = regex.findall(self.script)
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

    def add_event(self, time, event, start=False, check_continuous_processes = True):
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
        event_ind = self.find_event_index((time[0], time[1]), start)
        new_script = self.script[:event_ind] + """
    %s""" % event + ";" + self.script[event_ind:]
        self.script = new_script

    def add_event_over_time_range(self, start, end, event):
        time_ind, NOT_AT_END = self.find_time_index((start, 'early'))
        if not NOT_AT_END:
            raise ValueError("End generation is after the end of the simulation.")
        # See if range is already in script.
        regex = re.compile("%i:%i {" % (start, end))
        ALREADY_IN_SCRIPT = regex.findall(self.script)
        if ALREADY_IN_SCRIPT:
            pass
        else:
            new_script = self.script[:time_ind] + """%i:%i {
    %s;
}
""" % (start, end, event) + "\n" + self.script[time_ind:]
            self.script = new_script 

    def add_continuous_process(self, time, event):
        """
        Adds an event to each generation in the specified range.
        """
        start = time[0]
        end = time[1]
        missing_gens = []
        gen_range = [i for i in np.arange(start, end + 1)]
        for gen in gen_range:
            if self.time_already_in_script((gen, 'early')):
                self.add_event((gen, 'early'), event)
            elif self.time_already_in_script((gen, 'late')):
                self.add_event((gen, 'late'), event, gen == self.final_gen)
            else:
                missing_gens.append(gen)
        # Process the remaining generations. See helpful answer at
        # https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
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
        
#     # why isn't this working - check out later
#     # def add_simulation_end(self, generation):
#     #     command_to_save_out = "sim.treeSeqOutput(\"%s\")" % self.outfile
#     #     self.add_event(generation, "late", command_to_save_out)
#     #     self.add_event(generation, "late", "sim.simulationFinished()")

    def add_reference_population(self, popConfig, popLabel):
        if not isinstance(popConfig, msprime.PopulationConfiguration):
            TypeError("popConfig must be a msprime.PopulationConfiguration object.")
        sample_size = popConfig.sample_size
        initial_size = popConfig.initial_size
        growth_rate = np.exp(popConfig.growth_rate)
        if popConfig.growth_rate != 0:
            raise ValueError('At this stage, only the admixed population can have continuously changing population size.')
        self.sample_sizes.append(sample_size)
        self.initial_sizes.append(initial_size)
        self.growth_rates.append(growth_rate)
        self.population_labels.append(popLabel)
        ind = len(self.population_labels) - 1
        self.add_event((1, 'early'), "sim.addSubpop(\"p%i\", %i)" % (ind, initial_size))
        if growth_rate != 1:
            self.add_continuous_process((1,self.final_gen), 
                event = "newSize = asInteger(p\"%i\".individualCount * %f" % (ind, growth_rate))

    def add_admixed_population(self, popConfig, popLabel):
        if not isinstance(popConfig, msprime.PopulationConfiguration):
            raise TypeError("popConfig must be a msprime.PopulationConfiguration object.")
        sample_size = popConfig.sample_size
        initial_size = popConfig.initial_size
        growth_rate = np.exp(popConfig.growth_rate)
        # if popConfig.growth_rate != 0:
        #     return ValueError('At this stage, only the admixed population can have continuously changing population size.')
        self.sample_sizes.append(sample_size)
        self.initial_sizes.append(initial_size)
        self.growth_rates.append(growth_rate)
        self.population_labels.append(popLabel)
        ind = len(self.population_labels) - 1
        self.add_event((1, 'early'), "sim.addSubpop(\"p%i\", %i)" % (ind, initial_size))
        if growth_rate != 1:
            self.add_continuous_process((2,self.final_gen), 
                event = "newSize = asInteger(p%i.individualCount * %f" % (ind, growth_rate))

    def is_continuous(time):
        return(time.find(":"))


