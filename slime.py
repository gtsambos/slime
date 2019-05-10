import os
import numpy as np
import pandas as pd

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





