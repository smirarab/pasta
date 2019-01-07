#! /usr/bin/env python

import unittest
import datetime
import logging
import os

from pasta.test import get_testing_configuration, data_source_path, TestLevel, is_test_enabled

from pasta import get_logger
from pasta.alignment import Alignment
from pasta.scheduler import jobq, start_worker
from pasta.filemgr import TempFS

_LOG = get_logger(__name__)

config = get_testing_configuration()

start_worker(1)

class MergerTest(unittest.TestCase):
    def setUp(self):
        self.ts = TempFS()
        self.ts.create_top_level_temp(prefix='mergerTesting', parent=os.curdir)

    def tearDown(self):
        dir_list = self.ts.get_remaining_directories()
        for dir in dir_list:
            self.ts.remove_dir(dir)

    def get_merger(self, name):
        try:
            return config.create_merger(name=name, temp_fs=self.ts)
        except RuntimeError:
            _LOG.warn("""Could not create an merger of type %s !

This could indicate a bug in create_merger_using_config() or could mean that
your installation is not configured to run this tool.
""" % name)
            return None

    def testOpal(self):
        if is_test_enabled(TestLevel.SLOW, _LOG):
            self._impl_test_merger('opal')

    def testMuscle(self):
        if is_test_enabled(TestLevel.SLOW, _LOG):
            self._impl_test_merger('muscle')

    def _impl_test_merger(self, name):
        filename = data_source_path('merger1.fasta')
        alignment1 = Alignment()
        alignment1.read_filepath(filename, 'FASTA')
        filename = data_source_path('merger2.fasta')
        alignment2 = Alignment()
        alignment2.read_filepath(filename, 'FASTA')

        aln = self.get_merger('%s merger' % name)
        if aln is None:
            _LOG.warn("test%s skipped" % name)
            return
        a = aln.run(alignment1,
                    alignment2,
                    tmp_dir_par=self.ts.top_level_temp,
                    delete_temps=True)

        reference_fn = data_source_path('merger_result.fasta')
        reference_aln = Alignment()
        reference_aln.read_filepath(reference_fn, 'FASTA')
        self.assertEqual(reference_aln, a)

if __name__ == "__main__":
    unittest.main()
