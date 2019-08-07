
"""
Functions to assist with tests.
"""
import slime, msprime
import random
import pandas as pd
import numpy as np
import scipy.stats

def random_recombination_map(sequence_length, no_rows, filename=None, beta_pars=(1,5), sep=' '):
	"""
	Generates a random recombination map in the format required
	by the RecentHistory class.
	"""
	bp = list(random.sample(range(1, sequence_length), no_rows - 1))
	bp.sort()
	bp.append(sequence_length + 1)

	rates = scipy.stats.beta.rvs(a=beta_pars[0], b=beta_pars[1], size=no_rows)
	rates = rates.tolist()
	
	# print(type(rates))
	d = {'bp':bp, 'rates':rates}
	df = pd.DataFrame(data=d)
	# print(df)
	df.to_csv(path_or_buf=filename, sep=sep, index=False)


def basic_two_way_admixture(length=10, final_gen=20, prop=[0.3,0.7], rho=.001):
    config = slime.PopulationConfiguration(sample_size=0, initial_size=10)
    script = slime.RecentHistory(final_gen=final_gen, chrom_length=length,
            reference_configs=[config, config], adm_configs=config, recombination=rho,
            prop=prop)
    return(script)