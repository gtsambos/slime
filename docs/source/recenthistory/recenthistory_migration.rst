
.. _sec_recenthistory_migration:

*********
Migration
*********

Constant migration
******************

By default, ``slime`` assumes that the admixed populations is created by a single
pulse of migration from the reference populations.
This can be changed by specifying migration rates from each reference population
using the ``add_migration_rate`` method:

	>>> config = msprime.PopulationConfiguration(sample_size=10, initial_size=100)
    >>> script = slime.RecentHistory(final_gen=20, chrom_length=10,
    ...		recombination=.01, reference_configs=[config, config], 
    ...		adm_configs=config, prop=[0.3,0.7])
    >>> script.add_migration_rate(rates=[.01,.02])

To change migration rates at a particular time in the simulation, use the optional
argument ``time``. For example, adding this line to the above code will stop all
migration into the admixed population at the end of the 5th generation:

	>>> script.add_migration_rate(rates=[0, 0], time=(5, 'late'))


Mass migration
**************

You can also add one-off migrations into the admixed population using the
``add_mass_migration`` method.

	>>> script.add_mass_migration(prop=[.1, .1], time=(9, 'late'))

.. note:: At the moment, the ``add_mass_migration`` method will set migration rates
		  in the following generation to 0 (not to what they were in earlier 
		  generations of the simulation). You can set migration rates in generation 
		  n+1 to what they were before using ``add_migration_rate``.
