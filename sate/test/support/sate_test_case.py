#! /usr/bin/env python

import os
import unittest
import itertools
from cStringIO import StringIO

class SateTestCase(unittest.TestCase):

    def parse_fasta_file(self, file):
        if isinstance(file, str):
            file_stream = open(file, 'rU')
        else:
            file_stream = file
        line_iter = iter(file_stream)
        data = {}
        seq = StringIO()
        name = None
        for i, line in enumerate(line_iter):
            l = line.strip()
            if l.startswith('>'):
                if name:
                    data[name] = seq.getvalue().replace('-','')
                name = l[1:]
                seq = StringIO()
            else:
                seq.write(l)
        if name:
            data[name] = seq.getvalue().replace('-','')
        file_stream.close()
        return data

    def assertSameTaxa(self, sequence_dict1, sequence_dict2):
        self.assertEqual(sequence_dict1.keys().sort(), 
                         sequence_dict2.keys().sort())

    def assertSameSequences(self, sequence_dict1, sequence_dict2):
        self.assertEqual(sequence_dict1.values().sort(), 
                         sequence_dict2.values().sort())

    def assertSameDateSet(self, sequence_dict1, sequence_dict2):
        self.assertSameTaxa(sequence_dict1, sequence_dict2)
        self.assertSameSequences(sequence_dict1, sequence_dict2)
        for name, seq in sequence_dict1.iteritems():
            self.assertEqual(seq, sequence_dict2[name])
