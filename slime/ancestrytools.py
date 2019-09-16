#!/usr/bin/env python 
import tskit, msprime
import heapq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import collections


def check_row_lengths(left, right, population, child, ancestor):
    assert len(left) == len(right)
    assert len(right) == len(population)
    assert len(population) == len(child)
    assert len(ancestor) == len(child)


AncestryTableRow = collections.namedtuple(
    "AncestryTableRow",
    ["left", "right", "ancestor", "population", "child"])


class AncestryTable(object):
    """
    A table containing information about local ancestry.
    MORE LATER

    :ivar left: The array of left coordinates.
    :vartype left: numpy.ndarray, dtype=np.float64
    :ivar right: The array of right coordinates.
    :vartype right: numpy.ndarray, dtype=np.float64
    """

    def __init__(self):
        # super().__init__(ll_table, AncestryTableRow)
        self.left = np.array([], dtype=np.float64)
        self.right = np.array([], dtype=np.float64)
        self.ancestor = np.array([], dtype=np.int32)
        self.population = np.array([], dtype=np.int32)
        self.child = np.array([], dtype=np.int32)
        self.num_rows = self.num_rows()
        # self.ll_table = ll_table
        # ll_table = AncestryTable()
        # print(ll_table)

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
        return {
        "left": self.left,
        "right": self.right,
        "ancestor": self.ancestor,
        "population" : self.population,
        "child" : self.child
        }

    def num_rows(self):
        check_row_lengths(
            self.left, self.right, self.population, self.child, self.ancestor)
        return len(self.left)

    # def __eq__(self, other):
    #     ret = False
    #     if type(other) is type(self):
    #         # TODO: Expand on this later
    #         ret = True
    #     return ret

    def add_row(self, left, right, population, child, ancestor = -1):
        # left = left.astype(np.float64)
        # assert type(left) is np.float64
        # assert type(right) == np.float64
        # assert type(ancestor) == np.int32
        # assert type(population) == np.int32
        # assert type(child) == np.int32
        self.left = np.append(self.left, left)
        self.right = np.append(self.right, right)
        self.ancestor = np.append(self.ancestor, ancestor)
        self.population = np.append(self.population, population)
        self.child = np.append(self.child, child)
        self.num_rows += 1

    def set_columns(self, left, right, population, child, ancestor=None):
        if ancestor is None:
            ancestor = np.repeat(-1, len(left))
        check_row_lengths(left, right, population, child, ancestor)
        self.left = np.asarray(left, dtype = np.float64)
        self.right = np.asarray(right, dtype = np.float64)
        self.population = np.asarray(population, dtype=np.int32)
        self.child = np.asarray(child, dtype=np.int32)
        self.ancestor = np.asarray(ancestor, dtype=np.int32)
        self.num_rows = len(left)




    # @property
    # def left(self):
    #     return self.ll_table.left

    # @property
    # def right(self):
    #     return self.ll_table.right

    # @property
    # def ancestor(self):
    #     self.ll_table.ancestor

    # @property
    # def population(self):
    #     self.ll_table.ancestor

    # @property
    # def child(self):
    #     self.ll_table.child

    # def __str__(self):
    #     left = self.left
    #     right = self.right
    #     ancestor = self.ancestor
    #     population = self.population
    #     child = self.child
    #     ret = "id\tleft\t\tright\t\tparent\tchild\n"
    #     for j in range(self.num_rows):
    #         ret += "{}\t{:.8f}\t{:.8f}\t{}\t{}\n".format(
    #             j, left[j], right[j], parent[j], child[j])
    #     return ret[:-1]

    # def add_row(self, left, right, ancestor, population, child):
    #     return self.ll_table.add_row(left, right, ancestor, population, child)


# # Classes.
# class Segment(object):
#     """An ancestral segment mapping a given node ID to a half-open 
#     genomic interval [left, right)."""
    
#     def __init__(self, left, right, node):
#         assert left < right
#         self.left  = left
#         self.right = right
#         self.node  = node

#     def __lt__(self, other):
#         return (self.node, self.left, self.right) < (other.node, other.left, other.right)
    
#     def __repr__(self):
#         return repr((self.left, self.right, self.node))

# class Interval(object):
#     """A genomic interval [left, right]."""
    
#     def __init__(self, left, right):
#         assert left < right
#         self.left  = left
#         self.right = right
        
#     def __lt__(self, other):
#         return (self.left, self.right) < (other.left, other.right)
    
#     def __repr__(self):
#         return repr((self.left, self.right))

# # Classes.
# class IndividualSampleAncestry(object):

#     """Holds two SampleAncestry objects corresponding to the two
#     haplotypes of a given individual."""

#     def __init__(self, sample1, sample2, individual_id = None,
#         population_labels = None):
#         self.sample1 = sample1
#         self.sample2 = sample2
#         self.individual_id = None
#         self.population_labels = population_labels
#         if population_labels:
#             assert len(population_labels) == 2

#     def plot_individual_ancestry(self, plotFile):
#         # print(self.sample1.global_ancestry())
#         # print(self.sample2.global_ancestry())
#         def plot_ancestry_chunk(row, chrom):
#             l = row.left*1e-6
#             r = row.right*1e-6
#             p = row.population
#             if p == 0:
#                 c = 'blue'
#             elif p == 1:
#                 c = 'red'
#             else:
#                 c = 'white'
#             chunk = np.array([[l, 0], [r, 0], [r, 1], [l, 1]])
#             chrom.add_patch(Polygon(xy=chunk, color = c))

#         prop0 = self.sample1.global_ancestry()[0]
#         prop1 = self.sample2.global_ancestry()[0]
#         prop=[]
#         prop.append(round(np.mean([prop0[0], prop1[0]]),3))
#         prop.append(round(np.mean([prop0[1], prop1[1]]),3))

#         consolidated_node0 = self.sample1.sample_ancestry
#         consolidated_node1 = self.sample2.sample_ancestry

#         fig, (chr0, chr1) = plt.subplots(2)

#         fig.set_size_inches(25, 5)
#         if self.individual_id:
#             i = self.individual_id
#         else:
#             i = 0
#         if self.population_labels:
#             pop_labels = self.population_labels
#         else:
#             pop_labels = ['pop0', 'pop1']
#         fig.suptitle('Ancestry in admixed individual %i\nAncestry: %s = %.3f, %s = %.3f' % 
#             (i, pop_labels[0], prop[0], pop_labels[1], prop[1]))
#         fig.frameon=False
#         fig.legend(
#             handles = [Polygon(xy = np.array([[0,0],[0,1],[1,1],[1,0]]), color = 'blue'),
#                       Polygon(xy = np.array([[0,0],[0,1],[1,1],[1,0]]), color = 'red')],
#             labels = pop_labels,
#             loc = 'right'
#         )

#         for row in range(0, len(consolidated_node0)):
#             plot_ancestry_chunk(consolidated_node0.iloc[row], chr0)
#         chr0.set_xticks([10, 20, 30, 40, 50])
#         chr0.set_ylabel('Node 0')
#         chr0.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

#         for row in range(0, len(consolidated_node1)):
#             plot_ancestry_chunk(consolidated_node1.iloc[row], chr1)
#         chr1.set_xticks([10, 20, 30, 40, 50])
#         chr1.set_xlabel('Chromosomal position (Mb)')
#         chr1.set_ylabel('Node 1')
#         chr1.tick_params(left=False, labelleft=False)

#         plt.subplots_adjust(hspace=0)

#         plt.show()

#         if plotFile:
#             fig.savefig(plotFile, dpi=100)

# class SampleAncestry(object):
#     """ Holds local ancestry information for a node."""

#     def __init__(self, node_id, ancestry_table, populations, sequence_length):
#         self.node_id = node_id
#         self.populations = populations
#         self.num_populations = len(populations)
#         self.sequence_length = sequence_length
#         self.sample_ancestry = self.extract_sample_ancestry(ancestry_table)

#     def extract_sample_ancestry(self, ancestry_table):
#         tab = ancestry_table
#         hap = tab[tab.descendant == self.node_id]
#         hap.sort_values(by = 'left', inplace = True)
#         # Consolidate edges.
#         con = pd.DataFrame(
#             columns = ['left', 'right', 'population'])
#         row = 0
#         while row < len(hap):
#             r = hap.iloc[row]
#             right_endpt = r.right.copy()
#             row += 1
#             while row < len(hap) and hap.iloc[row].population == r.population:
#                 right_endpt = hap.iloc[row].right
#                 row += 1
#             new_row = pd.DataFrame([[r.left, right_endpt, r.population]],
#                 columns = ['left', 'right', 'population'])
#             con = con.append(new_row, ignore_index = True)
#         return(con)

#     def global_ancestry(self):
#         r = 0
#         global_ancestry = [0 for pop in self.populations]
#         while r < len(self.sample_ancestry):
#             row = self.sample_ancestry.iloc[r]
#             # index = self.populations.index(int(row.population))
#             # global_ancestry[index] += row.right - row.left
#             global_ancestry[int(row.population)] += row.right - row.left
#             r += 1
#         unassigned_ancestry = self.sequence_length - sum(global_ancestry)
#         global_ancestry = [p/self.sequence_length for p in global_ancestry]
#         unassigned_ancestry = unassigned_ancestry / self.sequence_length
#         return(global_ancestry, unassigned_ancestry)

#     def tract_lengths(self):
#         r = 0
#         tract_lengths = [[] for pop in self.populations]
#         while r < len(self.sample_ancestry):
#             row = self.sample_ancestry.iloc[r]
#             tract_lengths[int(row.population)].append(row.right - row.left)
#             r += 1
#         return(tract_lengths)

# class AncestryTable(object):
#     """A table holding local ancestry information for the sample."""

#     def __init__(self, ts, sample_nodes, populations):
#         self.ts              = ts
#         self.samples         = sample_nodes
#         self.num_samples     = len(sample_nodes)
#         self.populations     = populations
#         self.sequence_length = ts.sequence_length
#         self.ancestral_nodes = []
#         for n in self.ts.nodes():
#             if n.population in self.populations:
#                 self.ancestral_nodes.append(n.id)
#         self.No, self.Eo = self.ancestry_simplify()
#         self.ancestry_table = pd.DataFrame({
#             'left' : self.Eo.left,
#             'right' : self.Eo.right,
#             'ancestor' : self.Eo.parent,
#             'descendant' : self.Eo.child,
#             'population' : tuple(map(self.get_population, self.Eo.parent))
#             })

#     # Methods.
#     def global_ancestries(self):
#         global_ancestries = []
#         for node in self.samples:
#             sampleAnc = SampleAncestry(
#                 node, self.ancestry_table, self.populations, self.sequence_length)
#             global_ancestries.append(sampleAnc.global_ancestry())
#         return(global_ancestries)

#     def tract_lengths(self):
#         tract_lengths = []
#         for node in self.samples:
#             sampleAnc = SampleAncestry(
#                 node, self.ancestry_table, self.populations, self.sequence_length)
#             tract_lengths += sampleAnc.tract_lengths() 
#         return(tract_lengths)

#     def global_ancestry(self):
#         r = 0
#         global_ancestry = [0 for p in self.populations]
#         while r < len(self.ancestry_table):
#             row = self.ancestry_table.iloc[r]
#             global_ancestry[int(row.population)] += row.right - row.left
#             r += 1
#         global_ancestry = [ p/(self.sequence_length * self.num_samples) for p in global_ancestry]
#         unassigned_ancestry = 1 - sum(global_ancestry)
#         return(global_ancestry, unassigned_ancestry)

#     def get_sample_ancestry(self, node):
#         sampAnc = SampleAncestry(node, self.ancestry_table, 
#             self.populations, self.sequence_length)
#         return(sampAnc)

#     # Save output.
#     def save_summary_stats(self, outFile):
#         f = open(outFile, "w+")
#         print("Calculating and writing global ancestry to file...")
#         glob = self.global_ancestry()
#         for pop in self.populations:
#             f.write("Ancestry from population %i: %f\n" % (pop, glob[0][pop]))
#         f.close()

#     # Functions needed in initialisation.
#     def get_population(self, node_id):
#         return(self.No[node_id].population)

#     def ancestry_simplify(self):

#         Ni = self.ts.tables.nodes
#         Ei = self.ts.tables.edges
#         S = self.samples
#         ancestral_nodes = self.ancestral_nodes
#         L = self.sequence_length

#         # Initialising the algorithm.
#         No = tskit.NodeTable()
#         Eo = tskit.EdgeTable()
#         A = [[] for _ in range(len(Ni))]

#         # Create No and Eo.
#         # A: node is a sample
#         for u in S:
#             v    = No.add_row(time=Ni.time[u], flags=1, population=Ni.population[u])
#             A[u] = [Segment(0, L, v)]

#         # B and C: node not a sample
#         Q = []
#         v = -1

#         for ind, e in enumerate(Ei):
#             u = e.parent
#                 # TO DO: change this so that nodes not in
#                 # the final generation can be marked as samples too
#             for x in A[e.child]:
#                 if x.right > e.left and e.right > x.left:
#                     y = Segment(max(x.left, e.left), min(x.right, e.right), x.node)
#                     heapq.heappush(Q, y)

#             if ind % 10000 == 0:
#                 print('Processing edge', ind, 'of', len(Ei))
                    
#             if ind + 1 == len(Ei) or Ei[ind + 1].parent != u:
#                 # Process queue.  
                
#                 # C. node is ancestral.
#                 if u in ancestral_nodes:
                    
#                     Q2 = []
#                     while len(Q) > 0:
#                         x = heapq.heappop(Q)            
#                         if v == -1:
#                             v = No.add_row(time = Ni.time[u], population=Ni.population[u])
#                         child  = x.node
#                         x.node = v
#                         while len(Q) > 0 and Q[0].left <= x.right and Q[0].node == child:
#                             r = max(x.right, Q[0].right)
#                             x.right = r
#                             heapq.heappop(Q)

#                         Eo.add_row(x.left, x.right, v, child)
#                         Q2.append(Interval(x.left, x.right))        

#                     Q2.sort()
#                     while len(Q2) > 0:
#                         x = heapq.heappop(Q2)
#                         while len(Q2) > 0 and Q2[0].left <= x.right:
#                             r = max(x.right, Q2[0].right)
#                             x.right = r
#                             heapq.heappop(Q2)
#                         A[u].append(Segment(x.left, x.right, v))
                
                
#                 # B. node is internal to a branch
#                 else:
#                     while len(Q) > 0:
#                         x = heapq.heappop(Q)
#                         while len(Q) > 0 and Q[0].node == x.node and Q[0].left <= x.right:
#                             r = max(x.right, Q[0].right)
#                             x.right = r
#                             heapq.heappop(Q)             
#                         A[u].append(x)
                  
#                 # Reset node index.
#                 v = -1
#             """Sort the output edges and compact them as much as possible 
#             into the output table. We skip this for the algorithm listing
#             as it's pretty mundane. Could be replaced with calls to
#             squash_edges() and sort_tables()
#             """
            
#         E = list(Eo)
#         Eo.clear()
#         E.sort(key = lambda e: (e.parent, e.child, e.right, e.left))
#         start = 0
#         for j in range(1, len(E)):
#             condition = (
#                 E[j-1].right != E[j].left or
#                 E[j-1].parent != E[j].parent or
#                 E[j-1].child != E[j].child
#             )
#             if condition:
#                 Eo.add_row(E[start].left, E[j-1].right, E[j-1].parent, E[j-1].child)
#                 start = j
#         j = len(E)
#         Eo.add_row(E[start].left, E[j-1].right, E[j-1].parent, E[j-1].child)

#         return (No, Eo)