.. _sec_recenthistory:

Specifying admixture history
============================

Don't know how to program in Eidos?
Don't worry!
You can use our script generator to create SLiM scripts for use
with ``slime``! Hooray


.. _sec_recenthistory-defaults:

**********
Quickstart
**********

At minimum, your history of admixture should involve:

 * some number of discrete ancestral populations
 * a present-day admixed population
 * the age of the admixed population, in generations
 * an admixture proportion
 * the number of bases in the genomic segment being simulated
 * a recombination rate

 For each population, we must define a ``msprime.PopulationConfiguration`` object
 specifying the initial population size and desired sample size.

 For example, in the code below, I specify two ancestral populations, each with
 population size 100 at the time of initial admixture.

    >>> ref_pops = [msprime.PopulationConfiguration(sample_size=10, initial_size=100), 
    ...      msprime.PopulationConfiguration(sample_size=10, initial_size=100)]

 Our admixed population has initial population size 50, and is initially created with
 a single-pulse of admixture between the two reference populations that occurred 100
 generations ago.
 These populations contributed migrants to the admixed population with proportions of ``0.3`` and ``0.7`` respectively. 

    >>> adm_pop = msprime.PopulationConfiguration(sample_size=20, initial_size=50)
    >>> adm_prop = [0.3, 0.7]
    >>> gens = 100

 We input all this history into a ``slime.RecentHistory`` object, and
 specify a chromosome of length ``1e7`` and a constant recombination rate of
 ``1e-8``.

    >>> script = slime.RecentHistory(final_gen=gens, chrom_length=1e7,
    ...    recombination=1e-8, reference_configs=ref_pops, adm_configs=adm_pop,
    ...    prop=adm_prop)


Generating a SLiM script
************************

Use the ``dump_script()`` method:

    >>> script.dump_script('recent-history.slim')

You should then have a file called ``recent-history.slim`` in your root directory.

You can also view the SLiM script corresponding to your ``RecentHistory`` object at any
stage using the ``print_script()`` method:

    >>> script.print_script()


.. _sec_recenthistory-growth:

******
Growth
******

To specify a constant rate of population growth for any of the populations, just input
a ``growth`` parameter in the corresponding ``msprime.PopulationConfiguration`` object.

For instance, suppose our admixed population grows at a rate of ``exp(0.1)`` in eachh
new generation:

    >>> adm_pop = msprime.PopulationConfiguration(sample_size=20, initial_size=50, growth_rate=0.1)

