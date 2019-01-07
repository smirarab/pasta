#! /usr/bin/env python

import unittest
import datetime
import logging
import os, sys

from io import StringIO
from pasta import get_logger
from pasta.alignment import summary_stats_from_parse
from pasta.test import data_source_path

from dendropy import treesplit
import dendropy

_LOG = get_logger(__name__)

class DiagnoseDatatypeTest(unittest.TestCase):
    def testDiagnoseDNA(self):
        fp = data_source_path('small.fasta')
        print(fp)
        s = summary_stats_from_parse([fp], ["DNA", "RNA", "PROTEIN"], careful_parse=False)
        self.assertEqual(s[0], "DNA")
        self.assertEqual(s[1], [(32, 1650)])
        self.assertEqual(s[2], 32)
        self.assertEqual(s[3], True)
        s = summary_stats_from_parse([fp], ["DNA", "RNA", "PROTEIN"], careful_parse=True)
        self.assertEqual(s[0], "DNA")
        self.assertEqual(s[1], [(32, 1650)])
        self.assertEqual(s[2], 32)
        self.assertEqual(s[3], True)
        self.assertRaises(Exception, summary_stats_from_parse, [fp], ["RNA"], careful_parse=False)
        self.assertRaises(Exception, summary_stats_from_parse, [fp], ["RNA"], careful_parse=True)
    def testDiagnoseRNA(self):
        fp = data_source_path('smallrna.fasta')
        print(fp)
        s = summary_stats_from_parse([fp], ["DNA", "RNA", "PROTEIN"], careful_parse=False)
        self.assertEqual(s[0], "RNA")
        self.assertEqual(s[1], [(32, 1650)])
        self.assertEqual(s[2], 32)
        self.assertEqual(s[3], True)
        s = summary_stats_from_parse([fp], ["DNA", "RNA", "PROTEIN"], careful_parse=True)
        self.assertEqual(s[0], "RNA")
        self.assertEqual(s[1], [(32, 1650)])
        self.assertEqual(s[2], 32)
        self.assertEqual(s[3], True)
        self.assertRaises(Exception, summary_stats_from_parse, [fp], ["DNA", "PROTEIN"], careful_parse=False)
        _LOG.warn("WARNING: summary_stats_from_parse does not distinguish between RNA and DNA in 'careful' mode") 
        #self.assertRaises(Exception, summary_stats_from_parse, [fp], ["DNA", "PROTEIN"], careful_parse=True)
    def testDiagnoseProt(self):
        fp = data_source_path('caenophidia_mos.fasta')
        print(fp)
        s = summary_stats_from_parse([fp], ["DNA", "RNA", "PROTEIN"], careful_parse=False)
        self.assertEqual(s[0], "PROTEIN")
        self.assertEqual(s[1], [(114, 189)])
        self.assertEqual(s[2], 114)
        self.assertEqual(s[3], False)
        s = summary_stats_from_parse([fp], ["DNA", "RNA", "PROTEIN"], careful_parse=True)
        self.assertEqual(s[0], "PROTEIN")
        self.assertEqual(s[1], [(114, 189)])
        self.assertEqual(s[2], 114)
        self.assertEqual(s[3], False)
        self.assertRaises(Exception, summary_stats_from_parse, [fp], ["DNA", "RNA"], careful_parse=False)
        self.assertRaises(Exception, summary_stats_from_parse, [fp], ["DNA", "RNA"], careful_parse=True)
    def testDiagnoseBogus(self):
        fp = data_source_path('caenophidia_mos_bogus.fasta')
        self.assertRaises(Exception, summary_stats_from_parse, [fp], ["DNA", "RNA", "PROTEIN"], careful_parse=False)
        _LOG.warn("WARNING: summary_stats_from_parse does not distinguish between all bogus sequences in 'careful' mode") 
        #self.assertRaises(Exception, summary_stats_from_parse, [fp], ["DNA", "RNA", "PROTEIN"], careful_parse=True)

    def testDiagnoseMulti(self):
        multi_dir = data_source_path('testmulti/caenophidia')
        fp = os.path.join(multi_dir,'caenophidia_mos.fasta')
        fp2 = os.path.join(multi_dir,'caenophidia_mos2.fasta')
        s = summary_stats_from_parse([fp, fp2], ["DNA", "RNA", "PROTEIN"], careful_parse=False)
        self.assertEqual(s[0], "PROTEIN")
        self.assertEqual(s[1], [(114, 189), (109, 202)])
        self.assertEqual(s[2], 116) # two taxa names were changed and 5 were deleted, so the union is 116
        self.assertEqual(s[3], False)

        fp3 = data_source_path('smallrna.fasta')
        s = summary_stats_from_parse([fp3, fp3], ["DNA", "RNA", "PROTEIN"], careful_parse=False)
        self.assertEqual(s[0], "RNA")
        self.assertEqual(s[1], [(32, 1650),(32, 1650)])
        self.assertEqual(s[2], 32)
        self.assertEqual(s[3], True)
        self.assertRaises(Exception, summary_stats_from_parse, [fp, fp3], ["DNA", "RNA", "PROTEIN"], careful_parse=False)
        _LOG.warn("WARNING: summary_stats_from_parse will read multi with dna and protein as entirely protein. MIXED data type support is needed!") 


        fp4 = data_source_path('small.fasta')
        fp5 = data_source_path('smallunaligned.fasta')
        s = summary_stats_from_parse([fp4, fp4], ["DNA", "RNA", "PROTEIN"], careful_parse=False)
        self.assertEqual(s[0], "DNA")
        self.assertEqual(s[1], [(32, 1650),(32, 1650)])
        self.assertEqual(s[2], 32) 
        self.assertEqual(s[3], True)
        self.assertRaises(Exception, summary_stats_from_parse, [fp, fp3], ["DNA", "RNA", "PROTEIN"], careful_parse=False)
        _LOG.warn("WARNING: summary_stats_from_parse will read multi with dna and protein as entirely protein. MIXED data type support is needed!") 

        fp4 = data_source_path('small.fasta')
        fp5 = data_source_path('smallunaligned.fasta')
        s = summary_stats_from_parse([fp4, fp5], ["DNA", "RNA", "PROTEIN"], careful_parse=False)
        self.assertEqual(s[0], "DNA")
        self.assertEqual(s[1], [(32, 1650),(32, 1650)])
        self.assertEqual(s[2], 32) 
        self.assertEqual(s[3], False)
        self.assertRaises(Exception, summary_stats_from_parse, [fp, fp3], ["DNA", "RNA", "PROTEIN"], careful_parse=False)
        _LOG.warn("WARNING: summary_stats_from_parse will read multi with dna and protein as entirely protein. MIXED data type support is needed!") 

if __name__ == "__main__":
    unittest.main()
