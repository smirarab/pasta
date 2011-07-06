#! /usr/bin/env python

import unittest
import subprocess
import logging
import os, sys
from cStringIO import StringIO

from sate.test import get_testing_configuration, data_source_path, TestLevel, is_test_enabled
from sate import get_logger
from sate.mainsate import sate_main

_LOG = get_logger(__name__)
config = get_testing_configuration()

class MainTest(unittest.TestCase):
    def _main_execution(self, argv, stderr=None, stdout=None, rc=0):
        try:
            cmd = "import sys; from sate.mainsate import main_sate; main_sate(%s)[0] or sys.exit(1)" % repr(argv)
            invoc = [sys.executable, '-c', cmd]
            _LOG.debug("Command:\n\tpython -c " + repr(cmd))
            p = subprocess.Popen(invoc,
                                 stderr=subprocess.PIPE,
                                 stdout=subprocess.PIPE)
            (o, e) = p.communicate()
            r = p.wait()
            self.assertEquals(r, rc)
            if stderr is not None:
                self.assertEquals(e, stderr)
            if stdout is not None:
                self.assertEquals(o, stdout)
        except Exception, v:
            #self.assertEquals(str(v), 5)
            raise

    def testBasic(self):
        if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG):
            self._main_execution(['--hep'], rc=2)
            self._main_execution([], rc=1)
            self._main_execution(['--help'], rc=0)
_LOG.warn('SKIPPING multitest')
class A:
    def testMulti(self):
        if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG):
            self._main_execution(['-m', '-i', data_source_path('testmulti'), '--iter-limit=1'])



if __name__ == "__main__":
    unittest.main()
