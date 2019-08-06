
.. sec_recenthistory_mutations:

*********
Mutations
*********

To initialise a ``slime`` simulation with mutations, you make a ``MutationTypes`` object and feed it into your ``RecentHistory`` object via the `mutations` parameter.

.. note:: If you are appending your SLiM simulation to a msprime
         simulation, there's no need to code neutral mutations into your ``RecentHistory`` object, as they are more easily generated later using ``recapitate``. 


Fixed selection coefficients
****************************

Say we wish to specify 2 different types of mutations to occur:

1. Deleterious mutations with selection coefficient 0.2. These will occur with relative proportion 0.2, and have a dominance coefficient of 0.7.
2. Weakly beneficial mutations with selection coefficient 0.6. These will occur with a relative proportion 0.8, and have a dominance coefficient of 0.4.

The following code initialises a `MutationTypes` object that specifies these types of
mutations to occur in our forward-in-time simulation. 

    >>> muts = slime.MutationTypes(mutation_rate=.005,
    ...     selection_coeffs=[0.2,0.6], proportions=[0.2, 0.8], 
    ...     dominance_coeffs=[0.7, 0.4])

We'll then use this ``MutationTypes`` object to initialize a ``RecentHistory`` object via the ``mutations`` parameter.

    >>> config = msprime.PopulationConfiguration(sample_size=0, initial_size=10)
    >>> script = slime.RecentHistory(final_gen=20, chrom_length=100,
    ...     reference_configs=[config, config], adm_configs=config,
    ...     prop=[0.3,0.7], mutations=muts)


Random selection coefficients
*****************************

Coming soon.