.. _sec_ancestrytables_intro:


*************************
What are ancestry tables?
*************************

Intuitively, local ancestry might be shown with a matrix that looks something like this:

[make a picture and put it here]

A table that shows the most recent ancestral population from which a given
node has inherited on each segment of their genomes.

Each row of the ancestry table tells something about who has inherited from who
on a particular segment of their genome.
Specifically, each row (L, R, A, P, S) specifies that sample S has inherited from ancestor A
in population P over the genomic interval with coordinates (L, R).

Different ancestry tables can be produced from the same dataset depending on the set of ancestors
that are of interest.
We focus on applications of these tables to studies of local ancestry and IBD. 

An example
**********

Ancestry table with ancestors  from ancestral populations:


Uses
****

Local ancestry
--------------

If the set of ancestors in the table correspond to nodes from specified ancestral
populations, 


Identity-by-descent
-------------------




The API
*******

slime.AncestryTable()

    For a given set of samples and populations or ancestors, 
    MORE LATER

    :ivar left: The array of left coordinates.
    :vartype left: numpy.ndarray, dtype=np.float64
    :ivar right: The array of right coordinates.
    :vartype right: numpy.ndarray, dtype=np.float64