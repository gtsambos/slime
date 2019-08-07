.. _sec_recenthistory_quickstart:

**********
Quickstart
**********

SLiM scripts can be created from ``slime.RecentHistory`` objects.
The following code shows a simple, minimal demographic history of admixture:

    >>> ref_pops = [slime.PopulationConfiguration(initial_size=100), 
    ...      slime.PopulationConfiguration(initial_size=100)]

    >>> adm_pop = slime.PopulationConfiguration(initial_size=50)

    >>> script = slime.RecentHistory(final_gen=100, chrom_length=1e7,
    ...    recombination=1e-8, reference_configs=ref_pops, adm_configs=adm_pop,
    ...    prop=[0.3, 0.7])

For example, in the code below, I specify two ancestral populations, each with
population size 100 at the time of admixture.
Our admixed population has initial population size 50, and is initially created with
a single-pulse of admixture between the two reference populations that occurred 100
generations ago.
Each reference population contributes migrants to the admixed population with proportions
of ``0.3`` and ``0.7`` respectively. 
The simulated chromosome contains ``1e7`` bases, and recombination occurs at a uniform rate
of ``1e-8`` per base per generation.

At minimum, your history of admixture should involve:

 * some number of discrete ancestral populations
 * a present-day admixed population
 * the age of the admixed population, in generations
 * an admixture proportion
 * the number of bases in the genomic segment being simulated
 * a recombination rate

The sections below explain how more complicated admixtures can be represented with ``slime``.


Generating a SLiM script
************************

Use the ``dump_script()`` method:

    >>> script.dump_script('recent-history.slim')

You should then have a file called ``recent-history.slim`` in your root directory.

You can also view the SLiM script corresponding to your ``RecentHistory`` object at any
stage using the ``print_script()`` method:

    >>> script.print_script()