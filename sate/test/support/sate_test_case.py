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
                    data[name] = seq.getvalue().upper()
                name = l[1:]
                seq = StringIO()
            else:
                seq.write(l)
        if name:
            data[name] = seq.getvalue().upper()
        file_stream.close()
        return data

    def parseSequenceArg(self, seq_arg):
        if isinstance(seq_arg, dict):
            return seq_arg
        else:
            return self.parse_fasta_file(seq_arg)

    def remove_gaps(self, sequence_dict):
        sd = self.parseSequenceArg(sequence_dict)
        new_seqs = dict(sd)
        for name, seq in sd.iteritems():
            new_seqs[name] = seq.replace('-','')
        return new_seqs

    def concatenate_sequences(self, seq_data_list):
        taxa = set()
        data_sets = []
        for f in seq_data_list:
            seqs = self.parseSequenceArg(f)
            taxa.update(seqs.keys())
            data_sets.append(seqs)
        data = {}
        for t in taxa:
            data[t] = ''
        for ds in data_sets:
            for name in taxa:
                data[name] += ds.get(name, '')
        return data

    def assertSameTaxa(self, seq_data_list):
        if len(seq_data_list) < 2:
            return
        seqs1 = self.parseSequenceArg(seq_data_list[0])
        for i in range(1, len(seq_data_list)):
            seqs2 = self.parseSequenceArg(seq_data_list[i])
            self.assertEqual(sorted(seqs1.keys()), 
                             sorted(seqs2.keys()))

    def assertSameSequences(self, seq_data_list):
        seqs1 = self.parseSequenceArg(seq_data_list[0])
        sd1 = self.remove_gaps(seqs1)
        for i in range(1, len(seq_data_list)):
            seqs2 = self.parseSequenceArg(seq_data_list[i])
            sd2 = self.remove_gaps(seqs2)
            self.assertEqual(sorted(sd1.values()), 
                             sorted(sd2.values()))

    def assertSameDataSet(self, seq_data_list):
        seqs1 = self.parseSequenceArg(seq_data_list[0])
        sd1 = self.remove_gaps(seqs1)
        for i in range(1, len(seq_data_list)):
            seqs2 = self.parseSequenceArg(seq_data_list[i])
            sd2 = self.remove_gaps(seqs2)
            self.assertSameTaxa([sd1, sd2])
            self.assertSameSequences([sd1, sd2])
            for name, seq in sd1.iteritems():
                self.assertEqual(seq, sd2[name])

    def assertSameInputOutputSequenceData(self, 
            seq_data_list1, seq_data_list2):
        for i in range(len(seq_data_list1)):
            seqs1 = self.parseSequenceArg(seq_data_list1[i])
            seqs2 = self.parseSequenceArg(seq_data_list2[i])
            sd1 = self.remove_gaps(seqs1)
            sd2 = self.remove_gaps(seqs2)
            self.assertSameDataSet([sd1, sd2])

    def assertSameConcatenatedSequences(self, 
            concatenated_data, seq_data_list):
        concat_in = self.concatenate_sequences(sorted(seq_data_list))
        concat_out = self.parseSequenceArg(concatenated_data)
        sd_in = self.remove_gaps(concat_in)
        sd_out = self.remove_gaps(concat_out)
        self.assertSameSequences([sd_in, sd_out])

    def assertNoGapColumns(self, seq_data_list):
        for seq_data in seq_data_list:
            sd = self.parseSequenceArg(seq_data)
            columns_to_taxa = {}
            for name, seq in sd.iteritems():
                for column_index, residue in enumerate(seq):
                    if residue == '-':
                        if column_index not in columns_to_taxa.keys():
                            columns_to_taxa[column_index] = [name]
                        else:
                            columns_to_taxa[column_index].append(name)
            self.assertEqual(len(columns_to_taxa.keys()), len(set(columns_to_taxa.keys())))
            for col, name_list in columns_to_taxa.iteritems():
                self.assertEqual(len(name_list), len(set(name_list)))
                self.assertNotEqual(len(name_list), len(sd.keys()))

