#! /usr/bin/env python

import unittest
import logging
import os, sys

from sate.test import is_test_enabled, TestLevel, data_source_path
from sate.test.support.sate_test_case import SateTestCase
from sate import get_logger

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

if __name__ == "__main__":
    unittest.main()
