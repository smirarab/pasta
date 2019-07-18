#! /usr/bin/env python

import unittest
import datetime
import logging
import os, sys
import dendropy

from io import StringIO
from pasta.test import data_source_path
from pasta import treeholder
from pasta.errors import TaxaLabelsMismatchError
from pasta import get_logger

_LOG = get_logger(__name__)

class testCheckStartingTreeTaxonLabels(unittest.TestCase):
    def test_match(self):
        ds1 = dendropy.DataSet()
        ds1.read_from_string(
                ('>a\n'
                 'ATCG\n'
                 '>b\n'
                 'ATCG\n'
                 '>c\n'
                 'ATCG\n'),
                schema='fasta',
                data_type='dna')
        ds2 = dendropy.DataSet()
        tree_str = '((a,b),c);'
        ds2.read_from_string(
                tree_str,
                schema='newick')
        extra, missing = treeholder.check_taxon_labels(
                ds2.tree_lists[-1].taxon_set, ds1)
        self.assertEqual(len(missing), 0)
        self.assertEqual(len(extra), 0)

    def test_extra_tips(self):
        ds1 = dendropy.DataSet()
        ds1.read_from_string(
                ('>a\n'
                 'ATCG\n'
                 '>b\n'
                 'ATCG\n'
                 '>c\n'
                 'ATCG\n'),
                schema='fasta',
                data_type='dna')
        ds2 = dendropy.DataSet()
        tree_str = '((a,b),(c,d));'
        ds2.read_from_string(
                tree_str,
                schema='newick')
        extra, missing = treeholder.check_taxon_labels(
                ds2.tree_lists[-1].taxon_set, ds1)
        self.assertEqual(len(missing), 0)
        self.assertEqual(len(extra), 1)
        self.assertEqual(extra[0], 'd')
        tree_stream = StringIO()
        tree_stream.write(tree_str)
        self.assertRaises(TaxaLabelsMismatchError,
                treeholder.read_trees_into_dataset,
                ds1, tree_stream, starting_tree=True)

    def test_missing_tips(self):
        ds1 = dendropy.DataSet()
        ds1.read_from_string(
                ('>a\n'
                 'ATCG\n'
                 '>b\n'
                 'ATCG\n'
                 '>c\n'
                 'ATCG\n'
                 '>d\n'
                 'ATCG\n'),
                schema='fasta',
                data_type='dna')
        ds2 = dendropy.DataSet()
        tree_str = '((a,b),c);'
        ds2.read_from_string(
                tree_str,
                schema='newick')
        extra, missing = treeholder.check_taxon_labels(
                ds2.tree_lists[-1].taxon_set, ds1)
        self.assertEqual(len(missing), 1)
        self.assertEqual(len(extra), 0)
        self.assertEqual(missing[0], 'd')
        tree_stream = StringIO()
        tree_stream.write(tree_str)
        self.assertRaises(TaxaLabelsMismatchError,
                treeholder.read_trees_into_dataset,
                ds1, tree_stream, starting_tree=True)


if __name__ == "__main__":
    unittest.main()
