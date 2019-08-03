.. _sec_quickstart:

=============================
Quickstart - a simple example
=============================
 
We'll run through a basic ``slime`` simulation involving a simple two-population
example of admixture.

*******************************************
Step 1: Specify recent history of admixture
*******************************************

Suppose our admixture involves a single-pulse of admixture between two ancestral
populations 10 generations in the past. 
We'll specify configurations for each of these ancestral populations as an
``msprime.PopulationConfiguration`` object.
(Please see the `msprime documentation <`tskit documentation <https://tskit.readthedocs.io/en/stable>`_>`_ for more information about these).

	>>> ref0_config = msprime.PopulationConfiguration(sample_size=0, initial_size=10, growth_rate=0)
	>>> ref1_config = msprime.PopulationConfiguration(sample_size=0, initial_size=15, growth_rate=0)

We'll also specify a configuration for the admixed population, as well as an admixture proportion from each of the reference populations.

	>>> adm_config = msprime.PopulationConfiguration(sample_size=5, initial_size=10, growth_rate = 0)
	>>> adm_props = [0.3, 0.7]

By default, ``slime`` will assume a single-pulse of admixture between all inputted populations, but we can change this later if we wish.
We'll need to decide some other parameters: in particular, a recombination rate (presumed constant at the moment), a sequence length, and the number of generations ago at which the admixture occurred.

	>>> rho = 0.1
	>>> length = 10
	>>> gens = 15

This is the minimal input required to run a ``slime`` simulation.

To generate a SLiM script with this information, we'll create a ``slime.RecentHistory`` object containing all of these parameters as input:

	>>> script = slime.RecentHistory(final_gen=gens, chrom_length=length,
            reference_configs=[ref0_config, ref1_config], adm_configs=adm_config,
            prop=adm_props)