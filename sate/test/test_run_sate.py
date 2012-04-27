#! /usr/bin/env python

import unittest
import os
import sys
import subprocess
import logging

import sate
from sate.test import get_testing_configuration, data_source_path, TestLevel, is_test_enabled, TESTS_DIR
from sate import get_logger
from sate.filemgr import TempFS

_LOG = get_logger(__name__)
config = get_testing_configuration()

class RunSateTest(unittest.TestCase):
    def setUp(self):
        self.script_path = os.path.join(sate.sate_home_dir(), 'run_sate.py')
        self.ts = TempFS()
        self.ts.create_top_level_temp(prefix='runSateTest',
                parent=TESTS_DIR)
        self.anolis_file = data_source_path('anolis.fasta')
        self.multi_dir = data_source_path('testmulti')

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
        self.assertEquals(exit_code, return_code)
        if stdout != None:
            self.assertEquals(o, stdout)
        if stderr != None:
            self.assertEquals(e, stderr)

    def testBasic(self):
        self._exe_run_sate(['--blah'], return_code=2)
        self._exe_run_sate([], return_code=1)
        self._exe_run_sate(['-h'], return_code=0)

    def testSingleLocusRun(self):
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

    def testMultLocuRun(self):
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

if __name__ == "__main__":
    unittest.main()
