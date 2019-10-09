#!/usr/bin/env python 
import tskit, msprime
import pandas as pd
import numpy as np
import collections


def check_row_lengths(left, right, population, child, ancestor):
    assert len(left) == len(right)
    assert len(right) == len(population)
    assert len(population) == len(child)
    assert len(ancestor) == len(child)


# AncestryTableRow = collections.namedtuple(
#     "AncestryTableRow",
#     ["left", "right", "ancestor", "population", "child"])


class AncestryTable(object):
    """
    A table showing all genomic segments of the specified sample IDs
    that have ancestry with one of the specified populations.
    Each row (L, R, A, P, S) indicates that over the genomic interval
    with coordinates (L, R), the sample node with ID S has inherited
    from the ancestral node with ID A in population P.

    :ivar left: The array of left coordinates.
    :vartype left: numpy.ndarray, dtype=np.float64
    :ivar right: The array of right coordinates.
    :vartype right: numpy.ndarray, dtype=np.float64
    :ivar ancestor: The array of ancestral nodes.
    :vartype ancestor: numpy.ndarray, dtype=np.int32
    :ivar population: The array of population labels.
    :vartype population: numpy.ndarray, dtype=np.int32
    """

    def __init__(self):
        # super().__init__(ll_table, AncestryTableRow)
        self.left = np.array([], dtype=np.float64)
        self.right = np.array([], dtype=np.float64)
        self.ancestor = np.array([], dtype=np.int32)
        self.population = np.array([], dtype=np.int32)
        self.child = np.array([], dtype=np.int32)
        self.num_rows = self.num_rows()

    def __str__(self):
        left = self.left
        right = self.right
        ancestor = self.ancestor
        population = self.population
        child = self.child
        ret = "id\tleft\t\tright\t\tancestor\tpopulation\tchild\n"
        for j in range(self.num_rows):
            ret += "{}\t{:.8f}\t{:.8f}\t{}\t{}\t{}\n".format(
                j, left[j], right[j], ancestor[j], population[j], child[j])
        return ret[:-1]

    def asdict(self):
        """
        Returns a dictionary of table values.
        The keys are the column names, and the values are the
        numpy arrays holding the column values.
        """
        return {
        "left": self.left,
        "right": self.right,
        "ancestor": self.ancestor,
        "population" : self.population,
        "child" : self.child
        }

    def num_rows(self):
        """
        Returns the number of rows in the table.
        """
        check_row_lengths(
            self.left, self.right, self.population, self.child, self.ancestor)
        return len(self.left)

    def add_row(self, left, right, population, child, ancestor = -1):
        """
        Adds a single row with the specified values to the bottom of the table.

        :ivar left: The left coordinate of the segment.
        :vartype left: float
        :ivar right: The right coordinate of the segment.
        :vartype right: float
        :ivar population: The population of the ancestral node.
        :vartype population: int
        :ivar child: The ID of the child node.
        :vartype child: int
        :ivar ancestor: The ID of the ancestral node.
        :vartype ancestor: int
        """
        self.left = np.append(self.left, left)
        self.right = np.append(self.right, right)
        self.ancestor = np.append(self.ancestor, ancestor)
        self.population = np.append(self.population, population)
        self.child = np.append(self.child, child)
        self.num_rows += 1

    def set_columns(self, left, right, population, child, ancestor=None):
        """
        Sets the values in each column of the table.
        This makes it possible to add the information from many rows
        all at once.

        :ivar left: The list of left coordinates.
        :vartype left: list, dtype=np.float64
        :ivar right: The list of right coordinates.
        :vartype right: list, dtype=np.float64
        :ivar ancestor: The list of ancestral nodes.
        :vartype ancestor: list, dtype=np.int32
        :ivar population: The list of population labels.
        :vartype population: list, dtype=np.int32
        """
        if ancestor is None:
            ancestor = np.repeat(-1, len(left))
        check_row_lengths(left, right, population, child, ancestor)
        self.left = np.asarray(left, dtype = np.float64)
        self.right = np.asarray(right, dtype = np.float64)
        self.population = np.asarray(population, dtype=np.int32)
        self.child = np.asarray(child, dtype=np.int32)
        self.ancestor = np.asarray(ancestor, dtype=np.int32)
        self.num_rows = len(left)


def get_ancestry_table(ts, populations, samples=None, keep_ancestors=False):
    """
    Returns an AncestryTable showing local ancestry information for the
    specified set of samples. 

    :ivar ts: The tree sequence containing the dataset.
    :vartype ts: tskit.TreeSequence
    :ivar populations: A list of ancestral population IDs of interest.
    :vartype populations: list, dtype=int
    :ivar samples: A list of sample node IDs of interest. If None, all samples in the inputted tree sequence.
    :vartype samples: list, dtype=int
    :param bool keep_ancestors: If True, ancestral node IDs are retained in the output.
    :return: The ancestry table listing the local ancestry of the genomic segments corresponding to the child nodes.
    :rtype: :class:slime.AncestryTable
    """
    # Extract ancestors with the given population labels.
    # TODO later: make this work with more flexible inputs.
    assert len(populations) > 0
    if samples is None:
        samples = ts.samples()
    ancestors = [u.id for u in ts.nodes() if u.population in populations]
    if len(ancestors) == 0:
        raise ValueError("There are no nodes with the given population IDs.")

    # Apply map_ancestors.
    ancestor_table = ts.tables.map_ancestors(samples=samples, ancestors=ancestors)

    # Copy relevant edges into a new EdgeTable.
    edges_to_squash = tskit.EdgeTable()
    for row in ancestor_table:
        if row.child in samples:
            sample_pop = ts.tables.nodes.population[row.child]
            if sample_pop in populations:
                pop = sample_pop
            else:
                pop = ts.tables.nodes[row.parent].population
            edges_to_squash.add_row(
                left=row.left, right=row.right, parent=pop, child=row.child)

    # squash
    edges_to_squash.squash()

    # sort - would be better without 
    df = pd.DataFrame(
        data = {'left': edges_to_squash.left,
                'right': edges_to_squash.right,
                'population': edges_to_squash.parent,
                'child': edges_to_squash.child
        })
    df.sort_values(by=['child', 'left', 'right', 'population'], inplace=True)

    # Change into an ancestry table.
    ret = AncestryTable()
    ret.set_columns(
        left=np.array(df['left']), right=np.array(df['right']), 
        population=np.array(df['population']), child=np.array(df['child']))

    return(ret)

