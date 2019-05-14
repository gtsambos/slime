import os
import numpy as np
import pandas as pd
import msprime, pyslim
import tskit # for now, this is my local version of tskit.

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
        out_directory = "", 
        populations_to_sample_from= None,
        sample_sizes = None,
        ):
        self.slim_script = slim_script
        self.out_directory = out_directory
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
        return(ts)

    def debugger(self):
        """ A debugger to run before the simulation. """
        # SLiM debugging
        slim = self.slim_script
        print('\nSLiM file:', slim)
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






