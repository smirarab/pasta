"""SATe/PASTA - Phylogenetic Tree Container, effectively a wrapper of dendropy.Tree"""

# This file is part of PASTA, and is forked form SATe

# PASTA, like SATe is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Jiaye Yu and Mark Holder, University of Kansas

import copy
from dendropy import Tree
from dendropy import Edge
from dendropy import Node
from dendropy import DataSet as Dataset
from dendropy.datamodel.treemodel import _convert_node_to_root_polytomy as convert_node_to_root_polytomy
from pasta import get_logger

_LOG = get_logger(__name__)

class PhylogeneticTree(object):
    """Data structure to store phylogenetic tree, wrapping dendropy.Tree."""
    def __init__(self, dendropy_tree):
        self._tree = dendropy_tree
        self.n_leaves = self.count_leaves()
        self._tree.seed_node.edge.tail_node = None
        self._tree.seed_node.edge.length = None
# uym2 added: for midpoint decomposition
        self._tree.annotated = False

    def count_leaves(self):
        return len(self._tree.leaf_nodes())

    def count_nodes(self):
        return len(self._tree.nodes())

    def calc_splits(self):
        n = self.count_leaves()
        for i in self._tree.postorder_edge_iter():
            nd = i.head_node
            if nd.is_leaf():
                i.num_leaves_below = 1
            else:
                i.num_leaves_below = sum([j.edge.num_leaves_below for j in nd.child_nodes()])

    def get_centroid_edge(self,spanning=False):
        """Get centroid edge"""
        root = self._tree.seed_node
        root_children = root.child_nodes()
        if root_children and (spanning or not hasattr(root_children[0].edge, "num_leaves_below")):
            self.calc_splits()
            n_leaves = self.count_leaves()
            if spanning and len(root_children) == 1:
                n_leaves += 1
        else:
            if root.edge:
                n_leaves = root.edge.num_leaves_below
            else:
                n_leaves = sum([j.edge.num_leaves_below for j in root_children])
        centroid_edge = None
        centroid_imbalance = n_leaves
        half_taxa = n_leaves/2
        if half_taxa == 0:
            half_taxa = 1
        for edge in self._tree.postorder_edge_iter():
            if edge.tail_node is None:
                continue
            n_descendants = edge.num_leaves_below
            imbalance = abs(half_taxa - n_descendants)
            if (imbalance < centroid_imbalance):
                centroid_edge = edge
                centroid_imbalance = imbalance
        assert centroid_edge is not None
        return centroid_edge

    def get_longest_internal_edge(self):
        longest_internal_edge = None
        longest_len = -1.0
        for edge in self._tree.postorder_edge_iter():
            if edge.tail_node is None:
                continue
            if (edge.length is not None) and edge.is_internal() and edge.length > longest_len:
                longest_internal_edge = edge
                longest_len = edge.length
        assert longest_internal_edge is not None
        return longest_internal_edge

    def get_longest_edge(self):
        longest_edge = None
        longest_len = -1.0
        for edge in self._tree.postorder_edge_iter():
            if edge.tail_node is None:
                continue
            if edge.length is not None and edge.length > longest_len:
                longest_edge = edge
                longest_len = edge.length
        assert longest_edge is not None
        return longest_edge


    #def get_midpoint_edge(self):
	#find_midpoint_edge	

    def get_adjacent_edges(self, e):
        he = [i for i in e.head_node.get_incident_edges() if i is not e]
        te = [i for i in e.tail_node.get_incident_edges() if i is not e]
        he.extend(te)
        return he

    def get_breaking_edge(self, option):
        if option.lower() == 'centroid':
            return self.get_centroid_edge()
        elif option.lower() == 'longest':
            return self.get_longest_edge()
	##### Ziyang & uym2 ################
	#elif option.lower() == 'midpoint': # will have it later
        #    return self.get_midpoint_edge()
	###################################tree1, tree2 = bisect_tree(self.tree, breaking_edge_style=option)		
        else:
            raise ValueError('Unknown break strategy "%s"' % option)

    def bipartition_by_edge(self, e):
        """Prunes the subtree that attached to the head_node of edge e and returns them as a separate tree."""

        t = self._tree
        nr = e.head_node
        assert e.tail_node is not None
        assert e.head_node is not None
        assert nr.parent_node is e.tail_node
        is_valid_tree(t)

        n = self.n_leaves
        potentially_deleted_nd = e.tail_node
        grandparent_nd = potentially_deleted_nd.parent_node
        e.tail_node.remove_child(nr, suppress_unifurcations=True)

        nr.edge.length = None # Length of bisected edge
        nr.parent_node = None
        convert_node_to_root_polytomy(nr)
        t1 = PhylogeneticTree(Tree(seed_node=nr))
        n1 = t1.n_leaves # temp we could speed this up, by telling the Phylogenetic tree how many leaves it has

        if hasattr(e, "num_leaves_below"):
            if grandparent_nd is None:
                old_root = potentially_deleted_nd
                if old_root.edge:
                    old_root.edge.num_leaves_below -= n1
            else:
                if potentially_deleted_nd in grandparent_nd.child_nodes():
                    potentially_deleted_nd.edge.num_leaves_below -= n1
                old_root = grandparent_nd
                if old_root.edge:
                    old_root.edge.num_leaves_below -= n1
                while old_root.parent_node:
                    old_root = old_root.parent_node
                    if old_root.edge:
                        old_root.edge.num_leaves_below -= n1
        else:
            old_root = grandparent_nd or potentially_deleted_nd
            while old_root.parent_node:
                old_root = old_root.parent_node

        t2 = PhylogeneticTree(Tree(seed_node=old_root))

        is_valid_tree(t1._tree)
        is_valid_tree(t2._tree)
        return t1, t2

    def leaf_node_names(self):
        leaves = self._tree.leaf_nodes()
        return [i.taxon.label for i in leaves]

    def __str__(self):
        return self.compose_newick()

    def compose_newick(self,suppress_rooting = True):
        return self._tree.as_string(schema="newick",suppress_rooting = suppress_rooting)
        #return self._tree.compose_newick()

    def read_tree_from_file(self, treefile, file_format):
        dataset = Dataset()
        dataset.read(open(treefile, 'rU'), schema=file_format)
        dendropy_tree = dataset.trees_blocks[0][0]
        self._tree = dendropy_tree
        self.n_leaves = self.count_leaves()

def node_formatter(n):
    return str(id(n))

def edge_formatter(e):
    return "%s %f " % (str(id(e)), e.length)

def is_valid_tree(t):
    assert t and t
    rc = t.seed_node.child_nodes()
    num_children = len(rc)
    if num_children == 0:
        return True
    elif num_children == 1:
        assert len(rc[0].child_nodes()) != 0
        return True
    elif num_children == 2:
        assert((not list(rc[0].child_nodes())) and (not list(rc[1].child_nodes())))
    return True
