#! /usr/bin/env python

import unittest
import datetime
import logging
import os
from optparse import OptionGroup
from optparse import OptionParser
from io import StringIO

from pasta import get_logger
from pasta.tree import PhylogeneticTree
from pasta.test import data_source_path

from dendropy import DataSet as Dataset
from dendropy.treesplit import encode_splits, count_bits, delete_outdegree_one
_LOG = get_logger(__name__)

class PhylogeneticTreeTest(unittest.TestCase):

    def testCentroidBipartition(self):
        treef = data_source_path('diffDecomp.nex')
        pt = self.phylogeneticTreeFromFile(treef, file_format='NEXUS')
        # pt.add_n_leaf_des_attr()
        self.assertEqual(pt.n_leaves, 484)

        e = pt.get_centroid_edge()
        subtree1, subtree2 = pt.bipartition_by_edge(e)

        leaf_num = [subtree1.n_leaves, subtree2.n_leaves]
        leaf_num.sort()
        self.assertEqual(leaf_num, [231, 253])

    def testLongestInternalBipartition(self):
        treef = data_source_path('small.tree')
        pt = self.phylogeneticTreeFromFile(treef, file_format='NEWICK')
        self.assertEqual(pt.n_leaves, 32)

        e = pt.get_longest_internal_edge()
        subtree1, subtree2 = pt.bipartition_by_edge(e)

        leaf_num = [subtree1.n_leaves, subtree2.n_leaves]
        leaf_num.sort()
        self.assertEqual(leaf_num, [2, 30])

    def testLongestBipartition(self):
        treef = data_source_path('small.tree')
        pt = self.phylogeneticTreeFromFile(treef, file_format='NEWICK')
        self.assertEqual(pt.n_leaves, 32)

        e = pt.get_longest_edge()
        subtree1, subtree2 = pt.bipartition_by_edge(e)

        leaf_num = [subtree1.n_leaves, subtree2.n_leaves]
        leaf_num.sort()
        self.assertEqual(leaf_num, [1, 31])

    def phylogeneticTreeFromFile(self, treefile, file_format):
        dataset = Dataset()
        dataset.read(open(treefile, 'rU'), schema=file_format)
        dendropy_tree = dataset.tree_lists[0][0]
        tree = PhylogeneticTree(dendropy_tree)
        tree.calc_splits()
        delete_outdegree_one(tree._tree)
        return tree

if __name__ == "__main__":
    unittest.main()
