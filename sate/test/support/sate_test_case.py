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
                    data[name] = seq.getvalue().replace('-','').upper()
                name = l[1:]
                seq = StringIO()
            else:
                seq.write(l)
        if name:
            data[name] = seq.getvalue().replace('-','').upper()
        file_stream.close()
        return data

    def concatenate_sequences(self, file_path_list):
        taxa = set()
        data_sets = []
        file_path_list.sort()
        for f in file_path_list:
            seqs = self.parse_fasta_file(f)
            taxa.update(seqs.keys())
            data_sets.append(seqs)
        data = {}
        for t in taxa:
            data[t] = ''
        for ds in data_sets:
            for name in taxa:
                data[name] += ds.get(name, '')
        return data

    def assertSameTaxa(self, sequence_dict1, sequence_dict2):
        self.assertEqual(sequence_dict1.keys().sort(), 
                         sequence_dict2.keys().sort())

    def assertSameSequences(self, sequence_dict1, sequence_dict2):
        self.assertEqual(sequence_dict1.values().sort(), 
                         sequence_dict2.values().sort())

    def assertSameDataSet(self, sequence_dict1, sequence_dict2):
        self.assertSameTaxa(sequence_dict1, sequence_dict2)
        self.assertSameSequences(sequence_dict1, sequence_dict2)
        for name, seq in sequence_dict1.iteritems():
            self.assertEqual(seq, sequence_dict2[name])

    def assertSameInputOutputSequenceData(self, 
            file_path_list1, file_path_list2):
        for i in range(len(file_path_list1)):
            seqs1 = self.parse_fasta_file(file_path_list1[i])
            seqs2 = self.parse_fasta_file(file_path_list2[i])
            self.assertSameDataSet(seqs1, seqs2)

    def assertSameConcatenatedSequences(self, 
            concatenated_file_path, file_path_list):
        concat_in = self.concatenate_sequences(file_path_list)
        concat_out = self.parse_fasta_file(concatenated_file_path)
        self.assertSameSequences(concat_in, concat_out)
