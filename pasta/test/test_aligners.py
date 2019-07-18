#! /usr/bin/env python

import unittest
import datetime
import logging
import sys
import os

from pasta.test import get_testing_configuration, data_source_path, TestLevel, is_test_enabled

from pasta import get_logger
from pasta.alignment import Alignment
from pasta.scheduler import jobq, start_worker
from pasta.filemgr import TempFS

_LOG = get_logger(__name__)

config = get_testing_configuration()

start_worker(1)

class AlignerTest(unittest.TestCase):
    def setUp(self):
        self.ts = TempFS()
        self.ts.create_top_level_temp(prefix='alignerTesting', parent=os.curdir)

    def tearDown(self):
        dir_list = self.ts.get_remaining_directories()
        for dir in dir_list:
            self.ts.remove_dir(dir)
    def get_aligner(self, name):
        try:
            return config.create_aligner(name=name, temp_fs=self.ts)
        except RuntimeError:
            _LOG.warn("""Could not create an aligner of type %s !

This could indicate a bug in create_aligner_using_config() or could mean that
your installation is not configured to run this tool.
""" % name)
            return None

    def _impl_test_aligner(self, name, fn):
        filename = data_source_path(fn)
        alignment = Alignment()
        alignment.read_filepath(filename, 'FASTA')

        aln = self.get_aligner('%s' % name)
        if aln is None:
            _LOG.warn("test%s skipped" % name)
            return
        a = aln.run(alignment,
                    tmp_dir_par=self.ts.top_level_temp,
                    delete_temps=True)

        reference_fn = data_source_path('%s.%s' % (name, fn))
        reference_aln = Alignment()
        reference_aln.read_filepath(reference_fn, 'FASTA')
        _LOG.debug('Checking results from %s against %s' % (name, reference_fn))
        if reference_aln != a:
            i = 1
            while True:
                nrfn  = reference_fn + '.' + str(i)
                if os.path.exists(nrfn):
                    reference_aln = Alignment()
                    reference_aln.read_filepath(nrfn, 'FASTA')
                    _LOG.debug('Checking results from %s against %s' % (name, nrfn))
                    if reference_aln == a:
                        self.assertEqual(reference_aln, a)
                        return True
                    i += 1
                else:
                    self.assertEqual(reference_aln, a)


    def testClustalW2(self):
        if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG):
            self._impl_test_aligner('clustalw2', 'anolis.fasta')

    def testOpal(self):
        if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG):
            self._impl_test_aligner('opal', 'anolis.fasta')

    def testMafft(self):
        if is_test_enabled(TestLevel.EXHAUSTIVE, _LOG):
            self._impl_test_aligner('mafft', 'anolis.fasta')




if __name__ == "__main__":
    unittest.main()
