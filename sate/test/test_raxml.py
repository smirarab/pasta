#! /usr/bin/env python

import unittest
import datetime
import logging
import os

from sate.test import get_testing_configuration, data_source_path, TestLevel, is_test_enabled
from sate import get_logger, log_exception
from sate.alignment import Alignment
from sate.scheduler import jobq, start_worker
from sate.filemgr import TempFS

_LOG = get_logger(__name__)

config = get_testing_configuration()
start_worker(1)

class TreeEstimatorTest(unittest.TestCase):
    def setUp(self):
        self.ts = TempFS()
        self.ts.create_top_level_temp(prefix='treeEstimatorTest', parent=os.curdir)

    def tearDown(self):
        dir_list = self.ts.get_remaining_directories()
        for dir in dir_list:
            self.ts.remove_dir(dir)
    def get_tree_estimator(self, name):
        try:
            return config.create_tree_estimator(name=name, temp_fs=self.ts)
        except RuntimeError:
            log_exception(_LOG)
            _LOG.warn("""Could not create an aligner of type %s !

This could indicate a bug in create_tree_estimator_using_config() or could mean that
your installation is not configured to run this tool.
""" % name)
            return None

    def _impl_test_tree_estimator(self, name, datatype, partitions):
        filename = data_source_path('anolis.fasta')
        alignment = Alignment()
        alignment.read_filepath(filename, 'FASTA')
        te = self.get_tree_estimator(name)
        if te is None:
            _LOG.warn("test%s skipped" % name)
            return
        alignment.datatype = datatype
        a = te.run(alignment=alignment,
                   partitions=partitions,
                   tmp_dir_par=self.ts.top_level_temp,
                   delete_temps=True)

    def testRaxml(self):
        filename = data_source_path('mafft.anolis.fasta')
        alignment = Alignment()
        alignment.read_filepath(filename, 'FASTA')

        if is_test_enabled(TestLevel.SLOW, _LOG):
            self._impl_test_tree_estimator('raxml', datatype="DNA", partitions=[("DNA", 1, 1456)])

if __name__ == "__main__":
    unittest.main()
