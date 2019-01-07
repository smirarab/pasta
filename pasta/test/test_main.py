#! /usr/bin/env python

import unittest
import logging
import os, sys

from pasta.configure import get_configuration
from pasta.test import is_test_enabled, TestLevel, data_source_path
from pasta.test.support.sate_test_case import SateTestCase
from pasta.errors import TaxaLabelsMismatchError
from pasta import get_logger

_LOG = get_logger(__name__)

class MainTest(SateTestCase):
    def setUp(self):
        self.set_up()

    def tearDown(self):
        self.tear_down()

    def testBasic(self):
       self._main_execution(['--hep'], rc=2)
       self._main_execution([], rc=1)
       self._main_execution(['--help'], rc=0)

    def testMulti(self):
        if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG,
                module_name=".".join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            self._main_execution(['-m',
                    '-i', data_source_path('testmulti'),
                    '-o', self.ts.top_level_temp,
                    '--temporaries=%s' % self.ts.top_level_temp,
                    '-j', self.job_name,
                    '--iter-limit=1'])

class TestTaxonLabelMismatch(SateTestCase):
    def setUp(self):
        self.set_up()
        self.data = data_source_path('tiny.fasta')
        self.tree = data_source_path('tiny_name_mismatch.tre')

    def tearDown(self):
        self.tear_down()

    def testStartingTreeLabelMismatch(self):
        cmd = ['-i', self.data,
               '-t', self.tree,
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-j', self.job_name,
               '--iter-limit=1']
        self.assertRaises(TaxaLabelsMismatchError, self._exe, cmd)

class TestLowerCaseCharacters(SateTestCase):
    def setUp(self):
        self.set_up()
        self.data = data_source_path('tiny.lowercase.fasta')

    def tearDown(self):
        self.tear_down()

    def testLowerCaseData(self):
        cmd = ['-i', self.data,
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-j', self.job_name,
               '--max-subproblem-size=2',
               '--iter-limit=1']
        self._main_execution(cmd, rc=0)

class TestUnicodePathCharacters(SateTestCase):
    def setUp(self):
        self.set_up()
        data_file = data_source_path('tiny.fasta')
        unicode_name = 'm\xe9ss\xfdp\xe4th'
        self.tmp_sub_dir = self.get_subdir(unicode_name)
        self.data_path = self.get_path(
                name=unicode_name + '.fasta',
                parent_dir=self.tmp_sub_dir)
        src = open(data_file, 'rU')
        out = open(self.data_path, 'w')
        for line in src:
            out.write(line)
        src.close()
        out.close()
    
    def tearDown(self):
        self.tear_down()

    def testUnicodePath(self):
        cmd = ['-i', self.data_path,
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-j', self.job_name,
               '--max-subproblem-size=2',
               '--iter-limit=1']
        self._exe_run_sate(cmd, rc=0)

class TestSpacesInPath(SateTestCase):
    def setUp(self):
        self.set_up()
        data_file = data_source_path('tiny.fasta')
        space_name = 'a path with a lot of spaces'
        self.tmp_sub_dir = self.get_subdir(space_name)
        self.data_path = self.get_path(
                name=space_name + '.fasta',
                parent_dir=self.tmp_sub_dir)
        src = open(data_file, 'rU')
        out = open(self.data_path, 'w')
        for line in src:
            out.write(line)
        src.close()
        out.close()
    
    def tearDown(self):
        self.tear_down()

    def testSpaces(self):
        cmd = ['-i', self.data_path,
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-j', self.job_name,
               '--max-subproblem-size=2',
               '--iter-limit=1']
        self._exe_run_sate(cmd, rc=0)

class TestRnaData(SateTestCase):
    def setUp(self):
        self.set_up()
        self.tiny_rna = data_source_path('tinyrna.fasta')
        self.small_rna = data_source_path('smallrna.fasta')
        self.small_tree = data_source_path('small.tree')
        self.tiny_aln_path = self.get_path(
                '.marker001.tinyrna.aln') 
        self.small_aln_path = self.get_path(
                '.marker001.smallrna.aln') 
        self.init_aln_path = self.get_path(
                '_temp_iteration_initialsearch_seq_alignment.txt')
        self.iter_aln_path = self.get_path(
                '_temp_iteration_0_seq_alignment.txt')
        self.cfg_path = self.get_path(
                '_temp_sate_config.txt')

    def tearDown(self):
        self.tear_down()

    def testDefaultError(self):
        cmd = ['-i', self.tiny_rna,
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-j', self.job_name,
               '--iter-limit=1']
        self.assertRaises(Exception, self._exe, cmd)

    def testDnaTypeError(self):
        cmd = ['-i', self.tiny_rna,
               '-d', 'dna',
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-j', self.job_name,
               '--iter-limit=1']
        self.assertRaises(Exception, self._exe, cmd)

    def testProteinTypeError(self):
        cmd = ['-i', self.tiny_rna,
               '-d', 'protein',
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-j', self.job_name,
               '--iter-limit=1']
        self.assertRaises(Exception, self._exe, cmd)

    def testTinyRna(self):
        cmd = ['-i', self.tiny_rna,
               '-d', 'rna',
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-j', self.job_name,
               '--keeptemp',
               '--max-subproblem-size=2',
               '--iter-limit=1']
        self._exe(cmd)
        self.assert_is_nuc(self.tiny_rna, 'RNA')
        self.assert_is_nuc(self.tiny_aln_path, 'RNA')
        self.assertSameInputOutputSequenceData(
                [self.tiny_rna],
                [self.tiny_aln_path])
        self.assertNoGapColumns([self.tiny_aln_path,
                self.init_aln_path,
                self.iter_aln_path])
        self.assert_is_nuc(self.init_aln_path, 'DNA')
        self.assert_is_nuc(self.iter_aln_path, 'DNA')
        self.assertSameSequences([
                self.tiny_rna,
                self.tiny_aln_path,
                self.convert_rna_to_dna(self.init_aln_path, reverse=True),
                self.convert_rna_to_dna(self.iter_aln_path, reverse=True)])
        cfg = get_configuration(self.cfg_path)
        self.assertEqual(cfg.commandline.datatype.upper(), 'RNA')

    def testSmallRna(self):
        if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG,
                module_name=".".join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            cmd = ['-i', self.small_rna,
                   '-t', self.small_tree,
                   '-d', 'rna',
                   '-o', self.ts.top_level_temp,
                   '--temporaries=%s' % self.ts.top_level_temp,
                   '-j', self.job_name,
                   '--keeptemp',
                   '--iter-limit=1']
            self._exe(cmd)
            self.assert_is_nuc(self.small_rna, 'RNA')
            self.assert_is_nuc(self.small_aln_path, 'RNA')
            self.assertSameInputOutputSequenceData(
                    [self.small_rna],
                    [self.small_aln_path])
            self.assertNoGapColumns([self.small_aln_path,
                    self.iter_aln_path])
            self.assert_is_nuc(self.iter_aln_path, 'DNA')
            self.assertSameSequences([
                    self.small_rna,
                    self.small_aln_path,
                    self.convert_rna_to_dna(self.iter_aln_path, reverse=True)])
            cfg = get_configuration(self.cfg_path)
            self.assertEqual(cfg.commandline.datatype.upper(), 'RNA')

    def testTinyRnaAutoError(self):
        cmd = ['-i', self.tiny_rna,
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-j', self.job_name,
               '--auto']
        self.assertRaises(Exception, self._exe, cmd)

    def testTinyRnaAuto(self):
        cmd = ['-i', self.tiny_rna,
               '-d', 'rna',
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-j', self.job_name,
               '--keeptemp',
               '--max-subproblem-size=2',
               '--auto']
        self._exe(cmd)
        self.assert_is_nuc(self.tiny_rna, 'RNA')
        self.assert_is_nuc(self.tiny_aln_path, 'RNA')
        self.assertSameInputOutputSequenceData(
                [self.tiny_rna],
                [self.tiny_aln_path])
        self.assertNoGapColumns([self.tiny_aln_path,
                self.init_aln_path,
                self.iter_aln_path])
        self.assert_is_nuc(self.init_aln_path, 'DNA')
        self.assert_is_nuc(self.iter_aln_path, 'DNA')
        self.assertSameSequences([
                self.tiny_rna,
                self.tiny_aln_path,
                self.convert_rna_to_dna(self.init_aln_path, reverse=True),
                self.convert_rna_to_dna(self.iter_aln_path, reverse=True)])
        cfg = get_configuration(self.cfg_path)
        self.assertEqual(cfg.commandline.datatype.upper(), 'RNA')

class TestRnaMultiLocus(SateTestCase):
    def setUp(self):
        self.set_up()
        self.multi_dir = data_source_path('testmulti/')
        self.multi_rna_dir = os.path.join(self.multi_dir, 'tinyrna')
        self.in_path1 = os.path.join(self.multi_rna_dir, 'tinyrna1.fasta')
        self.in_path2 = os.path.join(self.multi_rna_dir, 'tinyrna2.fasta')
        self.aln_path1 = self.get_path(
                '.marker001.tinyrna1.aln') 
        self.aln_path2 = self.get_path(
                '.marker002.tinyrna2.aln') 
        self.cfg_path = self.get_path(
                '_temp_sate_config.txt')
        self.concat_path = self.get_path(
                '_temp_iteration_0_seq_alignment.txt')

    def tearDown(self):
        self.tear_down()

    def testMultiDefaultError(self):
        cmd = ['-i', self.multi_rna_dir,
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-m',
               '-j', self.job_name,
               '--iter-limit=1']
        self.assertRaises(Exception, self._exe, cmd)

    def testMultiDnaTypeError(self):
        cmd = ['-i', self.multi_rna_dir,
               '-d', 'dna',
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-m',
               '-j', self.job_name,
               '--iter-limit=1']
        self.assertRaises(Exception, self._exe, cmd)

    def testMultiProteinTypeError(self):
        cmd = ['-i', self.multi_rna_dir,
               '-d', 'protein',
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-m',
               '-j', self.job_name,
               '--iter-limit=1']
        self.assertRaises(Exception, self._exe, cmd)

    def testMultiAutoError(self):
        cmd = ['-i', self.multi_rna_dir,
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-m',
               '-j', self.job_name,
               '--auto']
        self.assertRaises(Exception, self._exe, cmd)

    def testMultiTinyRna(self):
        cmd = ['-i', self.multi_rna_dir,
               '-d', 'rna',
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-m',
               '-j', self.job_name,
               '--max-subproblem-size=2',
               '--iter-limit=1']
        self._exe(cmd)
        self.assert_is_nuc(self.in_path1, 'RNA')
        self.assert_is_nuc(self.in_path2, 'RNA')
        self.assert_is_nuc(self.aln_path1, 'RNA')
        self.assert_is_nuc(self.aln_path2, 'RNA')
        self.assert_is_nuc(self.concat_path, 'DNA')
        self.assertSameInputOutputSequenceData(
                [self.in_path1, self.in_path2],
                [self.aln_path1, self.aln_path2])
        self.assertNoGapColumns([
                self.aln_path1,
                self.aln_path2,
                self.concat_path])
        self.assertSameConcatenatedSequences(
                concatenated_data=self.convert_rna_to_dna(
                        self.concat_path,
                        reverse=True),
                seq_data_list=[self.in_path1, self.in_path2])
        cfg = get_configuration(self.cfg_path)
        self.assertEqual(cfg.commandline.datatype.upper(), 'RNA')

class TestRnaDnaIdentity(SateTestCase):
    def setUp(self):
        self.set_up()
        self.dna = data_source_path('small.fasta')
        self.rna = data_source_path('smallrna.fasta')
        self.tree = data_source_path('small.tree')
        self.dna_tmp = self.get_subdir('dna')
        self.rna_tmp = self.get_subdir('rna')
        self.dna_aln = self.get_path(
                name='.marker001.small.aln',
                parent_dir=self.dna_tmp)
        self.dna_tree = self.get_path(
                name='.tre',
                parent_dir=self.dna_tmp)
        self.rna_aln = self.get_path(
                name='.marker001.smallrna.aln',
                parent_dir=self.rna_tmp)
        self.rna_tree = self.get_path(
                name='.tre',
                parent_dir=self.rna_tmp)
        self.dna_score = self.get_path(
                name='.score.txt',
                parent_dir=self.dna_tmp)
        self.rna_score = self.get_path(
                name='.score.txt',
                parent_dir=self.rna_tmp)
        self.dna_tmp_aln = self.get_path(
                name='_temp_iteration_0_seq_alignment.txt',
                parent_dir=self.dna_tmp)
        self.rna_tmp_aln = self.get_path(
                name='_temp_iteration_0_seq_alignment.txt',
                parent_dir=self.rna_tmp)

    def tearDown(self):
        self.tear_down()

    def testRnaDnaIdentity(self):
        if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG,
                module_name=".".join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            dna_cmd = ['-i', self.dna,
                   '-d', 'dna',
                   '-t', self.tree,
                   '-o', self.dna_tmp,
                   '--temporaries=%s' % self.dna_tmp,
                   '-j', self.job_name,
                   '--aligner=mafft',
                   '--merger=muscle',
                   '--tree-estimator=fasttree',
                   '--start-tree-search-from-current',
                   '--tree-estimator-model=-gtr -gamma -seed 1111',
                   '--iter-limit=1']
            self._exe_run_sate(dna_cmd)
            rna_cmd = ['-i', self.rna,
                   '-d', 'rna',
                   '-t', self.tree,
                   '-o', self.rna_tmp,
                   '--temporaries=%s' % self.rna_tmp,
                   '-j', self.job_name,
                   '--aligner=mafft',
                   '--merger=muscle',
                   '--tree-estimator=fasttree',
                   '--start-tree-search-from-current',
                   '--tree-estimator-model=-gtr -gamma -seed 1111',
                   '--iter-limit=1']
            self._exe_run_sate(rna_cmd)
            self.assert_is_nuc(self.dna_aln, 'DNA')
            self.assert_is_nuc(self.rna_aln, 'RNA')
            self.assertNoGapColumns([
                    self.dna_aln,
                    self.rna_aln])
            self.assertSameDataSet([
                    self.rna,
                    self.rna_aln,
                    self.convert_rna_to_dna(self.dna, reverse=True),
                    self.convert_rna_to_dna(self.dna_aln, reverse=True)])
            # self.assertSameScores([self.dna_score, self.rna_score])
            # self.assertSameTrees([self.dna_tree, self.rna_tree])
            self.assertSameFiles([self.dna_tmp_aln, self.rna_tmp_aln])

class TestMixedMultiLocus(SateTestCase):
    def setUp(self):
        self.set_up()
        self.multi_dir = data_source_path('testmulti/')
        self.multi_mixed_dir = os.path.join(self.multi_dir, 'mixed')
        self.in_path1 = os.path.join(self.multi_mixed_dir, 'tinydna.fasta')
        self.in_path2 = os.path.join(self.multi_mixed_dir, 'tinyrna.fasta')
        self.aln_path1 = self.get_path(
                '.marker001.tinydna.aln') 
        self.aln_path2 = self.get_path(
                '.marker002.tinyrna.aln') 
        self.cfg_path = self.get_path(
                '_temp_sate_config.txt')
        self.concat_path = self.get_path(
                '_temp_iteration_0_seq_alignment.txt')

    def tearDown(self):
        self.tear_down()

    def testMixedDefaultError(self):
        cmd = ['-i', self.multi_mixed_dir,
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-m',
               '-j', self.job_name,
               '--iter-limit=1']
        self.assertRaises(Exception, self._exe, cmd)

    def testMixedDnaTypeError(self):
        cmd = ['-i', self.multi_mixed_dir,
               '-d', 'dna',
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-m',
               '-j', self.job_name,
               '--iter-limit=1']
        self.assertRaises(Exception, self._exe, cmd)

    def testMixedRnaTypeError(self):
        cmd = ['-i', self.multi_mixed_dir,
               '-d', 'rna',
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-m',
               '-j', self.job_name,
               '--iter-limit=1']
        self.assertRaises(Exception, self._exe, cmd)

    def testMixedProteinTypeError(self):
        cmd = ['-i', self.multi_mixed_dir,
               '-d', 'protein',
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-m',
               '-j', self.job_name,
               '--iter-limit=1']
        self.assertRaises(Exception, self._exe, cmd)

    def testMixedAutoError(self):
        cmd = ['-i', self.multi_mixed_dir,
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-m',
               '-j', self.job_name,
               '--auto']
        self.assertRaises(Exception, self._exe, cmd)

if __name__ == "__main__":
    unittest.main()

