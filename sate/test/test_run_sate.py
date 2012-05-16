#! /usr/bin/env python

import unittest
import os
import sys
import subprocess
import logging

import sate
from sate.test import (get_testing_configuration, data_source_path, TestLevel,
        is_test_enabled, TESTS_DIR)
from sate.test.support.sate_test_case import SateTestCase
from sate import get_logger
from sate.filemgr import TempFS

_LOG = get_logger(__name__)
config = get_testing_configuration()

class RunSateTest(SateTestCase):
    def setUp(self):
        self.script_path = os.path.join(sate.sate_home_dir(), 'run_sate.py')
        self.ts = TempFS()
        self.ts.create_top_level_temp(prefix='runSateTest',
                parent=TESTS_DIR)
        self.anolis_file = data_source_path('anolis.fasta')
        self.caenophidia_file = data_source_path('caenophidia_mos.fasta')
        self.multi_dir = data_source_path('testmulti/')
        self.multi_aa_dir = os.path.join(self.multi_dir, 'caenophidia')
        self.figwasp_dir = os.path.join(self.multi_dir, 'figwasps')
        self.hummingbird_dir = os.path.join(self.multi_dir, 'hummingbirds')

    def tearDown(self):
        dir_list = self.ts.get_remaining_directories()
        for dir in dir_list:
            self.ts.remove_dir(dir)

    def _exe_run_sate(self, args, return_code=0, stdout=None, stderr=None):
        if isinstance(args, str):
            arg_list = args.split()
        else:
            arg_list = args
        cmd = ['python', self.script_path] + arg_list
        _LOG.debug("Command:\n\t" + " ".join(cmd))
        p = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
        o, e = p.communicate()
        exit_code = p.wait()
        if exit_code != return_code:
            _LOG.error("exit code (%s) did not match %s" % (exit_code,
                    return_code))
            _LOG.error("here is the stdout:\n%s" % o)
            _LOG.error("here is the stderr:\n%s" % e)
        self.assertEquals(exit_code, return_code)
        if stdout != None:
            self.assertEquals(o, stdout)
        if stderr != None:
            self.assertEquals(e, stderr)

    def testBasic(self):
        self._exe_run_sate(['--blah'], return_code=2)
        self._exe_run_sate([], return_code=1)
        self._exe_run_sate(['-h'], return_code=0)

    def testSingleDnaLocusRun(self):
        if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG,
                module_name=".".join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            arg_list = ['-d', 'dna',
                        '-j', 'satejob',
                        '--keepalignmenttemps',
                        '--keeptemp',
                        '--temporaries=%s' % self.ts.top_level_temp,
                        '--iter-limit=1',
                        '-o', self.ts.top_level_temp,
                        '-i', self.anolis_file,]
            self._exe_run_sate(arg_list, return_code=0)
            self.assertSameInputOutputSequenceData(
                    [self.anolis_file],
                    [os.path.join(self.ts.top_level_temp,
                            'satejob.marker001.anolis.aln')])
            self.assertNoGapColumns([os.path.join(self.ts.top_level_temp,
                    'satejob.marker001.anolis.aln')])

    def testSingleAminoAcidLocusRun(self):
        if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG,
                module_name=".".join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            arg_list = ['-d', 'protein',
                        '-j', 'satejob',
                        '--keepalignmenttemps',
                        '--keeptemp',
                        '--temporaries=%s' % self.ts.top_level_temp,
                        '--iter-limit=1',
                        '-o', self.ts.top_level_temp,
                        '-i', self.caenophidia_file,]
            self._exe_run_sate(arg_list, return_code=0)
            self.assertSameInputOutputSequenceData(
                    [self.caenophidia_file],
                    [os.path.join(self.ts.top_level_temp,
                            'satejob.marker001.caenophidia_mos.aln')])
            self.assertNoGapColumns([os.path.join(self.ts.top_level_temp,
                    'satejob.marker001.caenophidia_mos.aln')])

    def testMultiDnaLocusRun(self):
        if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG,
                module_name=".".join([self.__class__.__name__,
                        sys._getframe().f_code.co_name])):
            arg_list = ['-d', 'dna',
                        '-j', 'satejob',
                        '--keepalignmenttemps',
                        '--keeptemp',
                        '--temporaries=%s' % self.ts.top_level_temp,
                        '--iter-limit=1',
                        '-m',
                        '-o', self.ts.top_level_temp,
                        '-i', self.multi_dir,]
            self._exe_run_sate(arg_list, return_code=0)
            seqs_in1_path = os.path.join(self.multi_dir, '1.fasta')
            seqs_in2_path = os.path.join(self.multi_dir, '2.fasta')
            seqs_out1_path = os.path.join(self.ts.top_level_temp,
                    'satejob.marker001.1.aln')
            seqs_out2_path = os.path.join(self.ts.top_level_temp,
                    'satejob.marker002.2.aln')
            self.assertSameInputOutputSequenceData(
                    [seqs_in1_path, seqs_in2_path],
                    [seqs_out1_path, seqs_out2_path])
            concat_out = os.path.join(self.ts.top_level_temp,
                    'satejob_temp_iteration_0_seq_alignment.txt')
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
                        '-j', 'satejob',
                        '--keepalignmenttemps',
                        '--keeptemp',
                        '--temporaries=%s' % self.ts.top_level_temp,
                        '--iter-limit=1',
                        '-m',
                        '-o', self.ts.top_level_temp,
                        '-i', self.multi_aa_dir,]
            self._exe_run_sate(arg_list, return_code=0)
            seqs_in1_path = os.path.join(self.multi_aa_dir,
                    'caenophidia_mos.fasta')
            seqs_in2_path = os.path.join(self.multi_aa_dir,
                    'caenophidia_mos2.fasta')
            seqs_out1_path = os.path.join(self.ts.top_level_temp,
                    'satejob.marker001.caenophidia_mos.aln')
            seqs_out2_path = os.path.join(self.ts.top_level_temp,
                    'satejob.marker002.caenophidia_mos2.aln')
            self.assertSameInputOutputSequenceData(
                    [seqs_in1_path, seqs_in2_path],
                    [seqs_out1_path, seqs_out2_path])
            concat_out = os.path.join(self.ts.top_level_temp,
                    'satejob_temp_iteration_0_seq_alignment.txt')
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
                        '-j', 'satejob',
                        '--keepalignmenttemps',
                        '--keeptemp',
                        '--temporaries=%s' % self.ts.top_level_temp,
                        '--iter-limit=1',
                        '--start-tree-search-from-current',
                        '--treefile=%s' % os.path.join(
                                self.figwasp_dir,
                                'starting.tre'),
                        '--merger=muscle',
                        '--tree-estimator=fasttree',
                        '-m',
                        '-o', self.ts.top_level_temp,
                        '-i', self.figwasp_dir,]
            self._exe_run_sate(arg_list, return_code=0)
            seqs_in1_path = os.path.join(self.figwasp_dir,
                    'M1504.fasta')
            seqs_in2_path = os.path.join(self.figwasp_dir,
                    'M1505.fasta')
            seqs_out1_path = os.path.join(self.ts.top_level_temp,
                    'satejob.marker001.M1504.aln')
            seqs_out2_path = os.path.join(self.ts.top_level_temp,
                    'satejob.marker002.M1505.aln')
            self.assertSameInputOutputSequenceData(
                    [seqs_in1_path, seqs_in2_path],
                    [seqs_out1_path, seqs_out2_path])
            concat_out = os.path.join(self.ts.top_level_temp,
                    'satejob_temp_iteration_0_seq_alignment.txt')
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
                        '-j', 'satejob',
                        '--keepalignmenttemps',
                        '--keeptemp',
                        '--temporaries=%s' % self.ts.top_level_temp,
                        '--iter-limit=1',
                        '--start-tree-search-from-current',
                        '--treefile=%s' % os.path.join(
                                self.hummingbird_dir,
                                'starting.tre'),
                        '--merger=muscle',
                        '--tree-estimator=fasttree',
                        '-m',
                        '-o', self.ts.top_level_temp,
                        '-i', self.hummingbird_dir,]
            self._exe_run_sate(arg_list, return_code=0)
            seqs_in1_path = os.path.join(self.hummingbird_dir,
                    'AK1.fasta')
            seqs_in2_path = os.path.join(self.hummingbird_dir,
                    'bfib.fasta')
            seqs_in3_path = os.path.join(self.hummingbird_dir,
                    'nd2.fasta')
            seqs_in4_path = os.path.join(self.hummingbird_dir,
                    'nd4.fasta')
            seqs_out1_path = os.path.join(self.ts.top_level_temp,
                    'satejob.marker001.AK1.aln')
            seqs_out2_path = os.path.join(self.ts.top_level_temp,
                    'satejob.marker002.bfib.aln')
            seqs_out3_path = os.path.join(self.ts.top_level_temp,
                    'satejob.marker003.nd2.aln')
            seqs_out4_path = os.path.join(self.ts.top_level_temp,
                    'satejob.marker004.nd4.aln')
            self.assertSameInputOutputSequenceData(
                    [seqs_in1_path, seqs_in2_path,
                     seqs_in3_path, seqs_in4_path],
                    [seqs_out1_path, seqs_out2_path,
                     seqs_out3_path, seqs_out4_path])
            concat_out = os.path.join(self.ts.top_level_temp,
                    'satejob_temp_iteration_0_seq_alignment.txt')
            self.assertSameConcatenatedSequences(
                    concatenated_data=concat_out,
                    seq_data_list=[seqs_in1_path, seqs_in2_path,
                            seqs_in3_path, seqs_in4_path])
            self.assertNoGapColumns([seqs_out1_path, seqs_out2_path,
                    seqs_out3_path, seqs_out4_path, concat_out])       

if __name__ == "__main__":
    unittest.main()
