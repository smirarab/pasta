#! /usr/bin/env python

import unittest
import os
import sys
import subprocess
import logging

from pasta.test import (data_source_path, TestLevel,
        is_test_enabled)
from pasta.test.support.sate_test_case import SateTestCase
from pasta import get_logger

_LOG = get_logger(__name__)

class RunSateTest(SateTestCase):
    def setUp(self):
        self.set_up()
        self.anolis_file = data_source_path('anolis.fasta')
        self.caenophidia_file = data_source_path('caenophidia_mos.fasta')
        self.multi_dir = data_source_path('testmulti/')
        self.multi_aa_dir = os.path.join(self.multi_dir, 'caenophidia')
        self.figwasp_dir = os.path.join(self.multi_dir, 'figwasps')
        self.hummingbird_dir = os.path.join(self.multi_dir, 'hummingbirds')
        self.ambig_dna = data_source_path('small.ambiguities.fasta')
        self.ambig_dna_tree = data_source_path('small.tree')
        self.ambig_aa = data_source_path('caenophidia_mos.ambiguities.fasta')
        self.ambig_aa_tree = data_source_path('caenophidia_mos.tre')

    def tearDown(self):
        self.tear_down()

    def testBasic(self):
        self._exe_run_sate(['--blah'], rc=2)
        self._exe_run_sate([], rc=1)
        self._exe_run_sate(['-h'], rc=0)

    def testSingleDnaLocusRun(self):
        if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG,
                module_name=".".join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            arg_list = ['-d', 'dna',
                        '--temporaries=%s' % self.ts.top_level_temp,
                        '--iter-limit=1',
                        '-j', self.job_name,
                        '-o', self.ts.top_level_temp,
                        '-i', self.anolis_file,]
            self._exe_run_sate(arg_list, rc=0)
            self.assertSameInputOutputSequenceData(
                    [self.anolis_file],
                    [os.path.join(self.ts.top_level_temp,
                            self.job_name + '.marker001.anolis.aln')])
            self.assertNoGapColumns([os.path.join(self.ts.top_level_temp,
                    self.job_name + '.marker001.anolis.aln')])

    def testSingleAminoAcidLocusRun(self):
        if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG,
                module_name=".".join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            arg_list = ['-d', 'protein',
                        '--temporaries=%s' % self.ts.top_level_temp,
                        '--iter-limit=1',
                        '-j', self.job_name,
                        '-o', self.ts.top_level_temp,
                        '-i', self.caenophidia_file,]
            self._exe_run_sate(arg_list, rc=0)
            self.assertSameInputOutputSequenceData(
                    [self.caenophidia_file],
                    [os.path.join(self.ts.top_level_temp,
                            self.job_name + '.marker001.caenophidia_mos.aln')])
            self.assertNoGapColumns([os.path.join(self.ts.top_level_temp,
                    self.job_name + '.marker001.caenophidia_mos.aln')])

    def testMultiDnaLocusRun(self):
        if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG,
                module_name=".".join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            arg_list = ['-d', 'dna',
                        '--temporaries=%s' % self.ts.top_level_temp,
                        '--iter-limit=1',
                        '-m',
                        '-j', self.job_name,
                        '-o', self.ts.top_level_temp,
                        '-i', self.multi_dir,]
            self._exe_run_sate(arg_list, rc=0)
            seqs_in1_path = os.path.join(self.multi_dir, '1.fasta')
            seqs_in2_path = os.path.join(self.multi_dir, '2.fasta')
            seqs_out1_path = os.path.join(self.ts.top_level_temp,
                    self.job_name + '.marker001.1.aln')
            seqs_out2_path = os.path.join(self.ts.top_level_temp,
                    self.job_name + '.marker002.2.aln')
            self.assertSameInputOutputSequenceData(
                    [seqs_in1_path, seqs_in2_path],
                    [seqs_out1_path, seqs_out2_path])
            concat_out = os.path.join(self.ts.top_level_temp,
                    self.job_name + '_temp_iteration_0_seq_alignment.txt')
            self.assertSameConcatenatedSequences(
                    concatenated_data=concat_out,
                    seq_data_list=[seqs_in1_path, seqs_in2_path])
            self.assertNoGapColumns([seqs_out1_path, seqs_out2_path,
                    concat_out])

    def testMultiAminoAcidLocusRun(self):
        if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG,
                module_name=".".join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            arg_list = ['-d', 'protein',
                        '--temporaries=%s' % self.ts.top_level_temp,
                        '--iter-limit=1',
                        '-m',
                        '-j', self.job_name,
                        '-o', self.ts.top_level_temp,
                        '-i', self.multi_aa_dir,]
            self._exe_run_sate(arg_list, rc=0)
            seqs_in1_path = os.path.join(self.multi_aa_dir,
                    'caenophidia_mos.fasta')
            seqs_in2_path = os.path.join(self.multi_aa_dir,
                    'caenophidia_mos2.fasta')
            seqs_out1_path = os.path.join(self.ts.top_level_temp,
                    self.job_name + '.marker001.caenophidia_mos.aln')
            seqs_out2_path = os.path.join(self.ts.top_level_temp,
                    self.job_name + '.marker002.caenophidia_mos2.aln')
            self.assertSameInputOutputSequenceData(
                    [seqs_in1_path, seqs_in2_path],
                    [seqs_out1_path, seqs_out2_path])
            concat_out = os.path.join(self.ts.top_level_temp,
                    self.job_name + '_temp_iteration_0_seq_alignment.txt')
            self.assertSameConcatenatedSequences(
                    concatenated_data=concat_out,
                    seq_data_list=[seqs_in1_path, seqs_in2_path])
            self.assertNoGapColumns([seqs_out1_path, seqs_out2_path,
                    concat_out])

    def testFigWaspDataRun(self):
         if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG,
                module_name=".".join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            arg_list = ['-d', 'dna',
                        '--temporaries=%s' % self.ts.top_level_temp,
                        '--iter-limit=1',
                        '--start-tree-search-from-current',
                        '--treefile=%s' % os.path.join(
                                self.figwasp_dir,
                                'starting.tre'),
                        '--merger=muscle',
                        '--tree-estimator=fasttree',
                        '-m',
                        '-j', self.job_name,
                        '-o', self.ts.top_level_temp,
                        '-i', self.figwasp_dir,]
            self._exe_run_sate(arg_list, rc=0)
            seqs_in1_path = os.path.join(self.figwasp_dir,
                    'M1504.fasta')
            seqs_in2_path = os.path.join(self.figwasp_dir,
                    'M1505.fasta')
            seqs_out1_path = os.path.join(self.ts.top_level_temp,
                    self.job_name + '.marker001.M1504.aln')
            seqs_out2_path = os.path.join(self.ts.top_level_temp,
                    self.job_name + '.marker002.M1505.aln')
            self.assertSameInputOutputSequenceData(
                    [seqs_in1_path, seqs_in2_path],
                    [seqs_out1_path, seqs_out2_path])
            concat_out = os.path.join(self.ts.top_level_temp,
                    self.job_name + '_temp_iteration_0_seq_alignment.txt')
            self.assertSameConcatenatedSequences(
                    concatenated_data=concat_out,
                    seq_data_list=[seqs_in1_path, seqs_in2_path])
            self.assertNoGapColumns([seqs_out1_path, seqs_out2_path,
                    concat_out])       

    def testHummingBirdDataRun(self):
         if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG,
                module_name=".".join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            arg_list = ['-d', 'dna',
                        '--temporaries=%s' % self.ts.top_level_temp,
                        '--iter-limit=1',
                        '--start-tree-search-from-current',
                        '--treefile=%s' % os.path.join(
                                self.hummingbird_dir,
                                'starting.tre'),
                        '--merger=muscle',
                        '--tree-estimator=fasttree',
                        '-m',
                        '-j', self.job_name,
                        '-o', self.ts.top_level_temp,
                        '-i', self.hummingbird_dir,]
            self._exe_run_sate(arg_list, rc=0)
            seqs_in1_path = os.path.join(self.hummingbird_dir,
                    'AK1.fasta')
            seqs_in2_path = os.path.join(self.hummingbird_dir,
                    'bfib.fasta')
            seqs_in3_path = os.path.join(self.hummingbird_dir,
                    'nd2.fasta')
            seqs_in4_path = os.path.join(self.hummingbird_dir,
                    'nd4.fasta')
            seqs_out1_path = os.path.join(self.ts.top_level_temp,
                    self.job_name + '.marker001.AK1.aln')
            seqs_out2_path = os.path.join(self.ts.top_level_temp,
                    self.job_name + '.marker002.bfib.aln')
            seqs_out3_path = os.path.join(self.ts.top_level_temp,
                    self.job_name + '.marker003.nd2.aln')
            seqs_out4_path = os.path.join(self.ts.top_level_temp,
                    self.job_name + '.marker004.nd4.aln')
            self.assertSameInputOutputSequenceData(
                    [seqs_in1_path, seqs_in2_path,
                     seqs_in3_path, seqs_in4_path],
                    [seqs_out1_path, seqs_out2_path,
                     seqs_out3_path, seqs_out4_path])
            concat_out = os.path.join(self.ts.top_level_temp,
                    self.job_name + '_temp_iteration_0_seq_alignment.txt')
            self.assertSameConcatenatedSequences(
                    concatenated_data=concat_out,
                    seq_data_list=[seqs_in1_path, seqs_in2_path,
                            seqs_in3_path, seqs_in4_path])
            self.assertNoGapColumns([seqs_out1_path, seqs_out2_path,
                    seqs_out3_path, seqs_out4_path, concat_out])       

    def testDnaAmbiguousCharactersMafftFasttreeTrusted(self):
        if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG,
                module_name=".".join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            arg_list = ['-d', 'dna',
                        '--temporaries=%s' % self.ts.top_level_temp,
                        '--iter-limit=1',
                        '-j', self.job_name,
                        '-o', self.ts.top_level_temp,
                        '-i', self.ambig_dna,
                        '-t', self.ambig_dna_tree,
                        '--aligner=mafft',
                        '--merger=muscle',
                        '--tree-estimator=fasttree',]
            self._exe_run_sate(arg_list, rc=0)
            self.assertSameInputOutputSequenceData(
                    [self.ambig_dna],
                    [os.path.join(self.ts.top_level_temp,
                        self.job_name + '.marker001.small.ambiguities.aln')])
            self.assertNoGapColumns([os.path.join(self.ts.top_level_temp,
                    self.job_name + '.marker001.small.ambiguities.aln')])

    def testDnaAmbiguousCharactersClustalRaxmlUntrusted(self):
        if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG,
                module_name=".".join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            arg_list = ['-d', 'dna',
                        '--temporaries=%s' % self.ts.top_level_temp,
                        '--iter-limit=1',
                        '-j', self.job_name,
                        '-o', self.ts.top_level_temp,
                        '-i', self.ambig_dna,
                        '-t', self.ambig_dna_tree,
                        '--aligner=clustalw2',
                        '--merger=muscle',
                        '--tree-estimator=raxml',
                        '--untrusted',]
            self._exe_run_sate(arg_list, rc=0)
            self.assertSameInputOutputSequenceData(
                    [self.ambig_dna],
                    [os.path.join(self.ts.top_level_temp,
                        self.job_name + '.marker001.small.ambiguities.aln')])
            self.assertNoGapColumns([os.path.join(self.ts.top_level_temp,
                    self.job_name + '.marker001.small.ambiguities.aln')])

    def testProteinAmbiguousCharactersMafftFasttreeTrusted(self):
        if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG,
                module_name=".".join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            arg_list = ['-d', 'protein',
                        '--temporaries=%s' % self.ts.top_level_temp,
                        '--iter-limit=1',
                        '-j', self.job_name,
                        '-o', self.ts.top_level_temp,
                        '-i', self.ambig_aa,
                        '-t', self.ambig_aa_tree,
                        '--aligner=mafft',
                        '--merger=muscle',
                        '--tree-estimator=fasttree',]
            self._exe_run_sate(arg_list, rc=0)
            self.assertSameInputOutputSequenceData(
                    [self.ambig_aa],
                    [os.path.join(self.ts.top_level_temp,
                        self.job_name + '.marker001.caenophidia_mos.ambiguities.aln')])
            self.assertNoGapColumns([os.path.join(self.ts.top_level_temp,
                    self.job_name + '.marker001.caenophidia_mos.ambiguities.aln')])

    def testProteinAmbiguousCharactersClustalRaxmlUntrusted(self):
        if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG,
                module_name=".".join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            arg_list = ['-d', 'protein',
                        '--temporaries=%s' % self.ts.top_level_temp,
                        '--iter-limit=1',
                        '-j', self.job_name,
                        '-o', self.ts.top_level_temp,
                        '-i', self.ambig_aa,
                        '-t', self.ambig_aa_tree,
                        '--aligner=clustalw2',
                        '--merger=muscle',
                        '--tree-estimator=raxml',
                        '--untrusted',]
            self._exe_run_sate(arg_list, rc=0)
            self.assertSameInputOutputSequenceData(
                    [self.ambig_aa],
                    [os.path.join(self.ts.top_level_temp,
                        self.job_name + '.marker001.caenophidia_mos.ambiguities.aln')])
            self.assertNoGapColumns([os.path.join(self.ts.top_level_temp,
                    self.job_name + '.marker001.caenophidia_mos.ambiguities.aln')])

if __name__ == "__main__":
    unittest.main()

