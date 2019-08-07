.. _sec_recenthistory_recombination:

*************
Recombination
*************

Recombination rates are specified per base, per generation.
Can be specified for a ``slime.RecentHistory`` object via the ``recombination`` parameter.

Uniform recombination
*********************

Pass the rate to your ``slime.RecentHistory`` object via the ``recombination`` parameter.

For example, the following code indicates a recombination rate of `1e-8` per base per gen:

    >>> script = slime.RecentHistory(final_gen=100, chrom_length=1e7,
    ...    recombination=1e-8, reference_configs=ref_pops, adm_configs=adm_pop,
    ...    prop=[0.3, 0.7])


Non-uniform recombination
*************************

A path to a space-delimited text file containing positions and rates can also be 
passed to the ``recombination`` parameter. 
The first column of the file should contain (ordered) integers representing
base positions. 
The second column should contain recombination rates. The rate in the 'n'th row
is the per-base recombination rate for all bases between the base positions
given in the `n`th (inclusive) and `n+1`th (not inclusive) row of the file.

Special rows: The position in the final row must be the length of the simulated sequence + n, but note that this row isn't actually used (it's a weird quirk of SLiM)...
The first row is a header and will be ignored.

Here is an example of a valid recombination map for a sequence of length 10,
with a recombination hotspot between bases 6 and 10::

	position rate(cM/Mb)
	1 0.1
	6 0.5
	11 0.1

