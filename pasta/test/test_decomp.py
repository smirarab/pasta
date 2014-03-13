#! /usr/bin/env python

import unittest
import datetime
import logging
import os

from pasta.test import get_testing_configuration, data_source_path, TestLevel, is_test_enabled
from pasta import get_logger
from pasta.tree import PhylogeneticTree
from pasta.alignment import SequenceDataset
from pasta.treeholder import read_and_encode_splits
from pasta.pastaalignerjob import bisect_tree

TOL = 0.000001
_LOG = get_logger(__name__)
config = get_testing_configuration()

class DecompTest(unittest.TestCase):

    def testLongestEdge(self):
        sd = SequenceDataset()
        fp = data_source_path('100T.fasta')
        sd.read(open(fp, 'rU'), file_format='FASTA', datatype='DNA')
        fp = data_source_path('100T.tree')
        tree_list = read_and_encode_splits(sd.dataset, open(fp, "rU"))
        self.assertEqual(len(tree_list), 1)
        t = PhylogeneticTree(tree_list[0])
        self._do_test_longest(t)

    def testCentroidEdge(self):
        sd = SequenceDataset()
        fp = data_source_path('100T.fasta')
        sd.read(open(fp, 'rU'), file_format='FASTA', datatype='DNA')
        fp = data_source_path('100T.tree')
        tree_list = read_and_encode_splits(sd.dataset, open(fp, "rU"))
        self.assertEqual(len(tree_list), 1)
        t = PhylogeneticTree(tree_list[0])
        self._do_test_centroid(t)

    def _do_test_longest(self, t, level="1"):
        if t.n_leaves < 3:
            return
        before_br_len = [e.length for e in t._tree.preorder_edge_iter() if e.length]
        _LOG.debug("code=%s\n before = %s" % (level, t.compose_newick()))
        _LOG.debug(" after len(before_br_len) = %d" % (len(before_br_len)))
        num_real_edges_before = len(before_br_len)
        if len(t._tree.seed_node.child_nodes()) < 3:
            num_real_edges_before -= 1
        t1, t2 = bisect_tree(t, 'longest')
        after_1_br_len = [e.length for e in t1._tree.preorder_edge_iter() if e.length]
        after_2_br_len = [e.length for e in t2._tree.preorder_edge_iter() if e.length]
        num_branches_1 = len(after_1_br_len)
        num_branches_2 = len(after_2_br_len)
        if num_branches_2 == 0:
            num_branches_2 = 1
        if num_branches_1 == 0:
            num_branches_1 = 1
        expected_diff = 3
        if num_branches_2 == 1:
            expected_diff -= 2
        if num_branches_1 == 1:
            expected_diff -= 2
        # cherries are rooted, so they make 1 edge look like 2
        if len(t1._tree.seed_node.child_nodes()) == 2:
            num_branches_1 -= 1
        if len(t2._tree.seed_node.child_nodes()) == 2:
            num_branches_2 -= 1

        _LOG.debug(" after 1 = %s" % (t1.compose_newick()))
        _LOG.debug(" after num_branches_1 = %d" % (num_branches_1))
        _LOG.debug(" after 2 = %s" % (t2.compose_newick()))
        _LOG.debug(" after num_branches_2 = %d" % (num_branches_2))

        #self.assertEqual(len(before_br_len), expected_diff + num_branches_1 + num_branches_2)
        before_br_len.sort(reverse=True)
        before_br_len.pop(0)
        before_sum = sum(before_br_len)
        after_sum = sum(after_1_br_len) + sum(after_2_br_len)
        diff = before_sum - after_sum
        self.assertTrue(abs(diff) < TOL)
        if t1.n_leaves > 2:
            nl = level+ ".1"
            self._do_test_longest(t1, level=nl)
        if t2.n_leaves > 2:
            nl = level+ ".2"
            self._do_test_longest(t2, level=nl)

    def _do_test_centroid(self, t, level="1"):
        if t.n_leaves < 5:
            return

        t.calc_splits()

        t1, t2 = bisect_tree(t, 'centroid')
        assert t1.n_leaves + t2.n_leaves == t.n_leaves
        # indent = level.count(".")

        # print("==============\nInput tree has %s leaf nodes." % t.n_leaves)
        # print("Subtree 1 tree has %s leaf nodes." % t1.n_leaves)
        # print("Subtree 2 tree has %s leaf nodes." % t2.n_leaves)
        # print("==============\n")

        if t1.n_leaves > 2:
            nl = level+ ".1"

            self._do_test_centroid(t1, level=nl)
        if t2.n_leaves > 2:
            nl = level+ ".2"

            self._do_test_centroid(t2, level=nl)

if __name__ == "__main__":
    unittest.main()
