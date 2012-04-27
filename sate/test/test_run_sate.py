#! /usr/bin/env python

import unittest
import os
import sys
import subprocess
import logging

from sate.test import get_testing_configuration, data_source_path, TestLevel, is_test_enabled
from sate import get_logger
import sate

_LOG = get_logger(__name__)
config = get_testing_configuration()

class RunSateTest(unittest.TestCase):
    def setUp(self):
        self.script_path = os.path.join(sate.sate_home_dir(), 'run_sate.py')

    def tearDown(self):
        pass

    def _exe_run_sate(self, args, return_code=0, stdout=None, stderr=None):
        if isinstance(args, str):
            arg_list = args.split()
        else:
            arg_list = args
        cmd = ['python', self.script_path] + arg_list
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

if __name__ == "__main__":
    unittest.main()
