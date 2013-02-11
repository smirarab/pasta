#! /usr/bin/env python

import unittest
import logging
import os, sys

from sate.test import is_test_enabled, TestLevel, data_source_path
from sate.test.support.sate_test_case import SateTestCase
from sate.errors import TaxaLabelsMismatchError
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
               '--iter-limit=1']
        self._main_execution(cmd, rc=0)

class TestUnicodePathCharacters(SateTestCase):
    def setUp(self):
        self.set_up()
        data_file = data_source_path('tiny.fasta')
        unicode_name = u'm\xe9ss\xfdp\xe4th'
        self.tmp_sub_dir = self.ts.create_temp_subdir(
                parent=self.ts.top_level_temp,
                prefix=unicode_name)
        self.data_path = os.path.join(self.tmp_sub_dir,
                self.job_name + unicode_name + '.fasta')
        src = open(data_file, 'rU')
        out = open(self.data_path, 'w')
        for line in src:
            out.write(line)
        src.close()
        out.close()
    
    def tearDown(self):
        self.register_files()
        self.ts.remove_dir(self.tmp_sub_dir)
        self.tear_down()

    def testUnicodePath(self):
        cmd = ['-i', self.data_path,
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-j', self.job_name,
               '--iter-limit=1']
        self._exe_run_sate(cmd, rc=0)

class TestSpacesInPath(SateTestCase):
    def setUp(self):
        self.set_up()
        data_file = data_source_path('tiny.fasta')
        space_name = 'a path with a lot of spaces'
        self.tmp_sub_dir = self.ts.create_temp_subdir(
                parent=self.ts.top_level_temp,
                prefix=space_name)
        self.data_path = os.path.join(self.tmp_sub_dir,
                self.job_name + space_name + '.fasta')
        src = open(data_file, 'rU')
        out = open(self.data_path, 'w')
        for line in src:
            out.write(line)
        src.close()
        out.close()
    
    def tearDown(self):
        self.register_files()
        self.ts.remove_dir(self.tmp_sub_dir)
        self.tear_down()

    def testUnicodePath(self):
        cmd = ['-i', self.data_path,
               '-o', self.ts.top_level_temp,
               '--temporaries=%s' % self.ts.top_level_temp,
               '-j', self.job_name,
               '--iter-limit=1']
        self._exe_run_sate(cmd, rc=0)

if __name__ == "__main__":
    unittest.main()
