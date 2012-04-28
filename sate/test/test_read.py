#! /usr/bin/env python

import unittest
import datetime
import logging
import os, sys

from cStringIO import StringIO
from sate import get_logger
from sate.alignment import summary_stats_from_parse
from sate.test import data_source_path

from dendropy import treesplit
import dendropy

_LOG = get_logger(__name__)

class DiagnoseDatatypeTest(unittest.TestCase):
    def testDiagnoseDNA(self):
        fp = data_source_path('small.fasta')
        print fp
        s = summary_stats_from_parse([fp], ["DNA", "RNA", "PROTEIN"], careful_parse=False)
        self.assertEqual(s[0], "DNA")
        self.assertEqual(s[1], [(32, 1650)])
        self.assertEqual(s[2], 32)
        s = summary_stats_from_parse([fp], ["DNA", "RNA", "PROTEIN"], careful_parse=True)
        self.assertEqual(s[0], "DNA")
        self.assertEqual(s[1], [(32, 1650)])
        self.assertEqual(s[2], 32)
        self.assertRaises(Exception, summary_stats_from_parse, [fp], ["RNA"], careful_parse=False)
        self.assertRaises(Exception, summary_stats_from_parse, [fp], ["RNA"], careful_parse=True)
    def testDiagnoseRNA(self):
        fp = data_source_path('smallrna.fasta')
        print fp
        s = summary_stats_from_parse([fp], ["DNA", "RNA", "PROTEIN"], careful_parse=False)
        self.assertEqual(s[0], "RNA")
        self.assertEqual(s[1], [(32, 1650)])
        self.assertEqual(s[2], 32)
        s = summary_stats_from_parse([fp], ["DNA", "RNA", "PROTEIN"], careful_parse=True)
        self.assertEqual(s[0], "RNA")
        self.assertEqual(s[1], [(32, 1650)])
        self.assertEqual(s[2], 32)
        self.assertRaises(Exception, summary_stats_from_parse, [fp], ["DNA", "PROTEIN"], careful_parse=False)
        _LOG.warn("WARNING: summary_stats_from_parse does not distinguish between RNA and DNA in 'careful' mode") 
        #self.assertRaises(Exception, summary_stats_from_parse, [fp], ["DNA", "PROTEIN"], careful_parse=True)
    def testDiagnoseProt(self):
        fp = data_source_path('caenophidia_mos.fasta')
        print fp
        s = summary_stats_from_parse([fp], ["DNA", "RNA", "PROTEIN"], careful_parse=False)
        self.assertEqual(s[0], "PROTEIN")
        self.assertEqual(s[1], [(114, 189)])
        self.assertEqual(s[2], 114)
        s = summary_stats_from_parse([fp], ["DNA", "RNA", "PROTEIN"], careful_parse=True)
        self.assertEqual(s[0], "PROTEIN")
        self.assertEqual(s[1], [(114, 189)])
        self.assertEqual(s[2], 114)
        self.assertRaises(Exception, summary_stats_from_parse, [fp], ["DNA", "RNA"], careful_parse=False)
        self.assertRaises(Exception, summary_stats_from_parse, [fp], ["DNA", "RNA"], careful_parse=True)
    def testDiagnoseBogus(self):
        fp = data_source_path('caenophidia_mos_bogus.fasta')
        self.assertRaises(Exception, summary_stats_from_parse, [fp], ["DNA", "RNA", "PROTEIN"], careful_parse=False)
        _LOG.warn("WARNING: summary_stats_from_parse does not distinguish between all bogus sequences in 'careful' mode") 
        #self.assertRaises(Exception, summary_stats_from_parse, [fp], ["DNA", "RNA", "PROTEIN"], careful_parse=True)

if __name__ == "__main__":
    unittest.main()
