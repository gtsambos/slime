.. _sec_recenthistory_growth:

***********************
Population size changes
***********************

Constant growth
***************

To specify a constant rate of population growth for any of the populations, pass
an exponential ``growth`` parameter to the corresponding ``msprime.PopulationConfiguration`` object.

For instance, suppose our admixed population grows at a rate of ``exp(0.1)`` in each
new generation:

    >>> adm_pop = msprime.PopulationConfiguration(sample_size=20, initial_size=50,
    ...		growth_rate=0.1)


Instant population size changes
*******************************

Create a ``msprime.PopulationParametersChange`` object and add it to your `RecentHistory`
object using the ``add_size_change`` method.

For instance, suppose our admixed population has size 25 from generation 11 onwards:

	>>> size_change = msprime.PopulationParametersChange(time=11, initial_size=25,
	...		population_id=2)
    >>> script = slime.RecentHistory(final_gen=gens, chrom_length=1e7,
    ...    recombination=1e-8, reference_configs=ref_pops, adm_configs=adm_pop,
    ...    prop)
    >>> script.add_size_change(PopulationParametersChange=size_change)

Read more about ``msprime.PopulationParametersChange`` objects here_.

.. _here: https://msprime.readthedocs.io/en/stable/api.html#msprime.PopulationParametersChange


Changes in growth rate
**********************

.. note:: Fill this out.