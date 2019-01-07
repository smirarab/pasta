#! /usr/bin/env python

import unittest
import datetime
import logging
import os

from io import StringIO
from pasta import get_logger
from pasta.alignment import Alignment, SequenceDataset, MultiLocusDataset,\
    merge_in
from pasta.treeholder import read_and_encode_splits

from pasta.test import get_testing_configuration, data_source_path, TestLevel, is_test_enabled
from pasta.test.support.sate_test_case import SateTestCase

from dendropy import treesplit
import dendropy
import sys

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
        alignment.sub_alignment( list(alignment.keys())[0:2] ).write_unaligned_fasta(filename+'.partial.raw')

    def testConcatenateAlignments(self):
        filename1 = data_source_path('small.fasta')
        filename2 = data_source_path('small.fasta')
        a = Alignment()
        b = Alignment()
        a.datatype = "DNA"
        b.datatype = "DNA"
        a.read_filepath(filename1, 'FASTA')
        b.read_filepath(filename2, 'FASTA')

    def testMaxSequenceLength(self):
        a = Alignment()
        a['1'] = 'A--CG--T'
        a['2'] = 'AC----GT'
        a['3'] = 'A-C-G-T-'
        a['4'] = 'ACGT---T'
        self.assertEqual(a.max_sequence_length(), 5)
            
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

class MultiLocusDatasetTest(SateTestCase):
    def setUp(self):
        self.set_up()
        self.tmp_sub_dir = self.ts.create_temp_subdir(
                parent=self.ts.top_level_temp,
                prefix='MultiLocusDatasetTest')
        self.data_path = os.path.join(self.tmp_sub_dir,
                self.job_name + '_test.fasta')
        self.mlds = MultiLocusDataset()
    
    def tearDown(self):
        self.register_files()
        self.ts.remove_dir(self.tmp_sub_dir)
        self.tear_down()

    def _create_seq_file(self, seq_str):
        out = open(self.data_path, 'w')
        out.write(seq_str)
        out.close()

    def _parse_seq_dataset(self, sd):
        d = {}
        for taxon, char_vec in list(sd.dataset.char_matrices[0].items()):
            d[taxon.label] = ''.join([i for i in char_vec])
        return d

    def testRNAConversion(self):
        sf = StringIO()
        sf.write('>a\nAUGCAUGC\n')
        sf.write('>b\nAUGCAUGC\n')
        self._create_seq_file(sf.getvalue())
        sf.seek(0)
        seqs = self.parse_fasta_file(sf)
        self.mlds.read_files([self.data_path],
                datatype='RNA',
                file_format='FASTA')
        self.assertEqual(len(self.mlds), 1)
        self.assertSameDataSet([seqs, self._parse_seq_dataset(self.mlds[0])])
        self.mlds.convert_rna_to_dna()
        seqs = self.convert_rna_to_dna(seqs)
        for k, v in list(seqs.items()):
            self.assertTrue('T' in v)
        self.assertSameDataSet([seqs,
                self._parse_seq_dataset(self.mlds[0])])
        self.mlds.convert_rna_to_dna()
        self.assertSameDataSet([seqs,
                self._parse_seq_dataset(self.mlds[0])])
        self.mlds.convert_dna_to_rna()
        seqs = self.convert_rna_to_dna(seqs, reverse=True)
        for k, v in list(seqs.items()):
            self.assertFalse('T' in v)
        self.assertSameDataSet([seqs,
                self._parse_seq_dataset(self.mlds[0])])
        self.mlds.convert_rna_to_dna()
        seqs = self.convert_rna_to_dna(seqs)
        for k, v in list(seqs.items()):
            self.assertTrue('T' in v)
        self.assertSameDataSet([seqs,
                self._parse_seq_dataset(self.mlds[0])])

if __name__ == "__main__":
    unittest.main()
