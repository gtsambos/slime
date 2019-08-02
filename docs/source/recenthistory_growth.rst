.. _sec_recenthistory_growth:

******
Growth
******

To specify a constant rate of population growth for any of the populations, just input
a ``growth`` parameter in the corresponding ``msprime.PopulationConfiguration`` object.

For instance, suppose our admixed population grows at a rate of ``exp(0.1)`` in eachh
new generation:

    >>> adm_pop = msprime.PopulationConfiguration(sample_size=20, initial_size=50, growth_rate=0.1)