#! /usr/bin/env python

import unittest
import datetime
import logging
import os

from cStringIO import StringIO
from sate import get_logger
from sate.alignment import Alignment, SequenceDataset
from sate.treeholder import read_and_encode_splits

from sate.test import get_testing_configuration, data_source_path, TestLevel, is_test_enabled

from dendropy import treesplit
import dendropy

_LOG = get_logger(__name__)

class AlignmentTest(unittest.TestCase):
    def testAlignment(self):
        filename = data_source_path('small.fasta')
        alignment = Alignment()
        alignment.read_filepath(filename, 'FASTA')
        num_taxa = alignment.get_num_taxa()
        self.assertEqual(num_taxa, 32)
        alignment.write_filepath(filename+'.phy', 'PHYLIP')
        alignment.write_unaligned_fasta(filename+'.raw')
        alignment.sub_alignment( alignment.keys()[0:2] ).write_unaligned_fasta(filename+'.partial.raw')

    def testConcatenateAlignments(self):
        filename1 = data_source_path('small.fasta')
        filename2 = data_source_path('small.fasta')
        a = Alignment()
        b = Alignment()
        a.datatype = "DNA"
        b.datatype = "DNA"
        a.read_filepath(filename1, 'FASTA')
        b.read_filepath(filename2, 'FASTA')

class SeqDatasetTest(unittest.TestCase):

    #def testTaxonRelabeling(self):
    #    sd = SequenceDataset()
    #    fp = data_source_path('bad_names.fasta')
    #    sd.read(open(fp, 'rU'), file_format='FASTA', datatype='DNA')
    #    fp = data_source_path('bad_names.tree')
    #    tree_list = read_and_encode_splits(sd.dataset, open(fp, 'rU'))
    #    self.assertEqual(len(tree_list), 1)
    #    tree = tree_list[0]

    #    real_tree_string = """('a Bad name','a  nothe++-_=+r',('an!@#"$o^&*()}{_ther',(another,'a Badn,ame')))"""
    #    safe_tree_string = """(abadname,another1,(another2,(another,abadname1)))"""
    #    self.assertEqual(tree.compose_newick(), real_tree_string)

    #    alignment = sd.relabel_for_sate()

    #    k = list(alignment.keys())
    #    k.sort()
    #    self.assertEqual(k, ['abadname', 'abadname1', 'another', 'another1', 'another2'])
    #    _LOG.debug(tree.description(9))
    #    self.assertEqual(tree.compose_newick(), safe_tree_string)

    #    alt_safe_tree_string = """((abadname1,another1),another,(another2,abadname))"""
    #    alt_real_tree_string = """(('a Badn,ame','a  nothe++-_=+r'),another,('an!@#"$o^&*()}{_ther','a Bad name'))"""
    #    alt_tree = read_and_encode_splits(sd.dataset, StringIO(alt_safe_tree_string))[0]
    #    self.assertEqual(alt_tree.compose_newick(), alt_safe_tree_string)

    #    sd.restore_taxon_names()

    #    self.assertEqual(tree.compose_newick(), real_tree_string)
    #    self.assertEqual(alt_tree.compose_newick(), alt_real_tree_string)

    def test100T(self):
        sd = SequenceDataset()
        fp = data_source_path('100T.fasta')
        sd.read(open(fp, 'rU'), file_format='FASTA', datatype='DNA')
        fp = data_source_path('100T.tree')
        tree_list = read_and_encode_splits(sd.dataset, open(fp, "rU"))
        self.assertEqual(len(tree_list), 1)

    def test1000T(self):
        sd = SequenceDataset()
        fp = data_source_path('1000T.fasta')
        sd.read(open(fp, 'rU'), file_format='FASTA', datatype='DNA')
        fp = data_source_path('1000T.tree')
        tree_list = read_and_encode_splits(sd.dataset, open(fp, "rU"))
        self.assertEqual(len(tree_list), 1)

    def testDNAFasta(self):
        sd = SequenceDataset()
        fp = data_source_path('anolis.fasta')
        sd.read(open(fp, 'rU'), file_format='FASTA', datatype='DNA')

    #def testDeleteMissing(self):
    #    sd = SequenceDataset()
    #    fp = data_source_path('missing.fasta')
    #    sd.read(open(fp, 'rU'), file_format='FASTA', datatype='DNA')
    #    self.assertFalse(sd.sequences_are_valid())

    #    self.assertTrue(sd.sequences_are_valid(True, None))
    #    aln = sd.relabel_for_sate()

    #    fp = data_source_path('delmissing.fasta')
    #    sd.read(open(fp, 'rU'), file_format='FASTA', datatype='DNA')
    #    del_aln = sd.relabel_for_sate()
    #    for k, v in aln.iteritems():
    #        self.assertEqual(v.upper(), del_aln[k].upper())

    #def testMissingToAmbig(self):
    #    sd = SequenceDataset()
    #    fp = data_source_path('missing.fasta')
    #    sd.read(open(fp, 'rU'), file_format='FASTA', datatype='DNA')
    #    self.assertFalse(sd.sequences_are_valid())

    #    self.assertTrue(sd.sequences_are_valid(True, sd.alphabet.any_residue.symbol))
    #    aln = sd.relabel_for_sate()

    #    fp = data_source_path('ambigmissing.fasta')
    #    sd.read(open(fp, 'rU'), file_format='FASTA', datatype='DNA')
    #    ambig_aln = sd.relabel_for_sate()
    #    for k, v in aln.iteritems():
    #        self.assertEqual(v.upper(), ambig_aln[k].upper())

if __name__ == "__main__":
    unittest.main()
