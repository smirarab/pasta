#! /usr/bin/env python

import os
import sys
import unittest
import itertools
import re
import subprocess
import random
import string
from io import StringIO

import dendropy

import pasta
from pasta import get_logger
from pasta.test import TESTS_DIR
from pasta.filemgr import TempFS
from pasta.mainpasta import pasta_main

_LOG = get_logger(__name__)

class SateTestCase(unittest.TestCase):

    def set_up(self):
        self.ts = TempFS()
        self.ts.create_top_level_temp(prefix='runSateTest',
                parent=TESTS_DIR)
        self.job_name = 'satejob' + self.random_id(8)
        self.dirs = set([self.ts.top_level_temp])
        self.paths = set()

    def tear_down(self):
        self.register_files()
        self.remove_dir()

    def _main_execution(self, args, stdout=None, stderr=None, rc=0):
        try:
            cmd = "import sys; from pasta.mainpasta import pasta_main; pasta_main(%s)[0] or sys.exit(1)" % repr(args)
            invoc = [sys.executable, '-c', cmd]
            _LOG.debug("Command:\n\tpython -c " + repr(cmd))
            p = subprocess.Popen(invoc,
                                 stderr=subprocess.PIPE,
                                 stdout=subprocess.PIPE)
            (o, e) = p.communicate()
            r = p.wait()
            if r != rc:
                _LOG.error("exit code (%s) did not match %s" % (r,
                        rc))
                _LOG.error("here is the stdout:\n%s" % o)
                _LOG.error("here is the stderr:\n%s" % e)
            self.assertEqual(r, rc)
            if stderr is not None:
                self.assertEqual(e, stderr)
            if stdout is not None:
                self.assertEqual(o, stdout)
        except Exception as v:
            #self.assertEquals(str(v), 5)
            raise

    def _exe_run_sate(self, args, stdout=None, stderr=None, rc=0):
        script_path = os.path.join(pasta.pasta_home_dir(), 'run_sate.py')
        if isinstance(args, str):
            arg_list = args.split()
        else:
            arg_list = args
        cmd = ['python', script_path] + arg_list
        _LOG.debug("Command:\n\t" + " ".join(cmd))
        p = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
        o, e = p.communicate()
        exit_code = p.wait()
        if exit_code != rc:
            _LOG.error("exit code (%s) did not match %s" % (exit_code,
                    rc))
            _LOG.error("here is the stdout:\n%s" % o)
            _LOG.error("here is the stderr:\n%s" % e)
        self.assertEqual(exit_code, rc)
        if stdout != None:
            self.assertEqual(o, stdout)
        if stderr != None:
            self.assertEqual(e, stderr)

    def _exe(self, args):
        return pasta_main(args)

    def parse_fasta_file(self, file):
        if isinstance(file, str):
            _LOG.info('parsing fasta file {0!r}...'.format(file))
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
                seq.write(l.replace(' ', ''))
        if name:
            data[name] = seq.getvalue().upper()
        file_stream.close()
        return data
    
    def parse_score_file(self, file_obj):
        file_stream = file_obj
        if isinstance(file_obj, str):
            _LOG.info('parsing score file {0!r}...'.format(file_obj))
            file_stream = open(file_obj, 'rU')
        return(float(file_stream.read().strip()))

    def parse_score_arg(self, arg):
        if isinstance(arg, float):
            return arg
        return self.parse_score_file(arg)

    def parse_tree_arg(self, arg):
        if isinstance(arg, dendropy.Tree):
            return arg
        return self.parse_tree_file(arg)

    def parse_tree_file(self, file_obj):
        file_stream = file_obj
        if isinstance(file_obj, str):
            _LOG.info('parsing tree file {0!r}...'.format(file_obj))
            file_stream = open(file_obj, 'rU')
        t = dendropy.Tree()
        t.read_from_stream(file_stream, schema='newick')
        file_stream.close()
        return t

    def assertSameTrees(self, tree_list, percent_tol=1e-6):
        if len(tree_list) < 2:
            return
        tree1 = self.parse_tree_arg(tree_list.pop(0))
        for t in tree_list:
            tree2 = self.parse_tree_arg(t)
            self.assertEqualTreePair(tree1, tree2, percent_tol)

    def assertEqualTreePair(self, tree1, tree2, percent_tol=1e-6):
        self.assertEqual(sorted(tree1.taxon_set.labels()),
                         sorted(tree2.taxon_set.labels()))
        nodes1 = [n for n in tree1.postorder_node_iter()]
        nodes2 = [n for n in tree2.postorder_node_iter()]
        self.assertEqual(len(nodes1), len(nodes2))
        for i, n1 in enumerate(nodes1):
            n2 = nodes2[i]
            if n1.taxon is not None:
                self.assertTrue(n2.taxon is not None)
                self.assertEqual(n1.taxon.label, n2.taxon.label)
            else:
                self.assertEqual(n2.taxon, None)
            if n1.edge.length is not None:
                self.assertTrue(n2.edge.length is not None)
                self.assertApproxEqual(n1.edge.length, n2.edge.length,
                        percent_tol)
            else:
                self.assertEqual(n2.edge.length, None)

    def assertSameScores(self, score_list, percent_tol=1e-6):
        if len(score_list) < 2:
            return
        score1 = self.parse_score_arg(score_list.pop(0))
        for s in score_list:
            score2 = self.parse_score_arg(s)
            self.assertApproxEqual(score1, score2, percent_tol=percent_tol)

    def assertApproxEqual(self, x, y, percent_tol=1e-6):
        self.assertTrue(
                ((abs(x-y) / ((abs(x)+abs(y))/2))*100) < percent_tol)
        
    def parseSequenceArg(self, seq_arg):
        if isinstance(seq_arg, dict):
            return seq_arg
        else:
            return self.parse_fasta_file(seq_arg)

    def remove_gaps(self, sequence_dict):
        sd = self.parseSequenceArg(sequence_dict)
        new_sd = {}
        for name, seq in list(sd.items()):
            new_seq = re.sub(r'[-?]', '', seq)
            if new_seq != '':
                new_sd[name] = new_seq
        return new_sd

    def concatenate_sequences(self, seq_data_list):
        taxa = set()
        data_sets = []
        for f in seq_data_list:
            seqs = self.parseSequenceArg(f)
            taxa.update(list(seqs.keys()))
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
            for name, seq in list(sd1.items()):
                self.assertEqual(seq, sd2[name])

    def assertSameFiles(self, files):
        all_equal = True
        f1 = files.pop(0)
        if isinstance(f1, str):
            f1 = open(f1, 'rU')
        s1 = f1.read()
        for f2 in files:
            if isinstance(f2, str):
                f2 = open(f2, 'rU')
            s2 = f2.read()
            if not s1 == s2:
                all_equal = False
                _LOG.error('files {0!r} and {1!r} are different!'.format(
                        f1.name, f2.name))
        self.assertTrue(all_equal)

    def assertSameInputOutputSequenceData(self, 
            seq_data_list1, seq_data_list2):
        for i in range(len(seq_data_list1)):
            _LOG.debug("comparing %s to %s" % (seq_data_list1[i],
                    seq_data_list2[i]))
            seqs1 = self.parseSequenceArg(seq_data_list1[i])
            seqs2 = self.parseSequenceArg(seq_data_list2[i])
            self.assertSameDataSet([seqs1, seqs2])

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
            for name, seq in list(sd.items()):
                for column_index, residue in enumerate(seq):
                    if residue == '-':
                        if column_index not in list(columns_to_taxa.keys()):
                            columns_to_taxa[column_index] = [name]
                        else:
                            columns_to_taxa[column_index].append(name)
            self.assertEqual(len(list(columns_to_taxa.keys())), len(set(columns_to_taxa.keys())))
            for col, name_list in list(columns_to_taxa.items()):
                self.assertEqual(len(name_list), len(set(name_list)))
                self.assertNotEqual(len(name_list), len(list(sd.keys())))

    def random_id(self, length=8,
            char_pool=string.ascii_letters + string.digits):
        return ''.join(random.choice(char_pool) for i in range(length))

    def get_subdir(self, name, parent_dir=None):
        if not parent_dir:
            parent_dir = self.ts.top_level_temp
        if not parent_dir in self.dirs:
            raise Exception('{0!r} is not a registered test dir. '
                    'you should not create a test dir outside of the '
                    'unit test temp file system'.format(parent_dir))
        d = self.ts.create_temp_subdir(
                parent=parent_dir,
                prefix=self.job_name + name)
        self.dirs.add(d)
        self.register_dir(d)
        return d

    def get_path(self, name, parent_dir=None):
        if not parent_dir:
            parent_dir = self.ts.top_level_temp
        if not parent_dir in self.dirs:
            raise Exception('{0!r} is not a registered test dir. '
                    'you should not create a test file outside of the '
                    'unit test temp file system'.format(parent_dir))
        p = os.path.join(parent_dir, self.job_name + name)
        self.paths.add(p)
        return p

    def register_file(self, path):
        shared = os.path.commonprefix([self.ts.top_level_temp, path])
        if not shared == self.ts.top_level_temp:
            raise Exception('cannot register file outside of the unit test '
                    'temp file system')
        self.paths.add(path)

    def register_dir(self, d):
        self.ts._directories_created_lock.acquire()
        self.ts._directories_created.add(d)
        self.ts._directories_created_lock.release()

    def dir_registered(self, d):
        self.ts._directories_created_lock.acquire()
        b = d in self.ts._directories_created
        self.ts._directories_created_lock.release()
        return b

    def register_files(self):
        _LOG.debug('registering temp file system...')
        self.ts.run_generated_filenames.extend(self.paths)
        for d in self.dirs:
            self.register_dir(d)
        self.register_dir(os.path.join(
                self.ts.top_level_temp, self.job_name))
        for path, dirs, files in os.walk(self.ts.top_level_temp):
            for f in files:
                if f.startswith(self.job_name):
                    self.ts.run_generated_filenames.append(
                            os.path.join(path, f))
                elif os.path.basename(f) in self.ts.run_generated_filenames:
                    self.ts.run_generated_filenames.append(
                            os.path.join(path, f))
                else:
                    if not f in self.ts.run_generated_filenames:
                        _LOG.warning('could not register {0!r}'.format(
                                os.path.join(path, f)))
            for d in dirs:
                self.register_dir(os.path.join(path, d))

    def remove_dir(self, d=None):
        if not d:
            d = self.ts.top_level_temp
        _LOG.debug('removing {0!r}...'.format(d))
        self.ts.remove_dir(d)

    def convert_rna_to_dna(self, seqs, reverse=False):
        seq_dict = self.parseSequenceArg(seqs)
        d = {}
        for taxon, seq in list(seq_dict.items()):
            if reverse:
                d[taxon] = seq.replace('T', 'U')
            else:
                d[taxon] = seq.replace('U', 'T')
        return d

    def assert_is_nuc(self, seqs, datatype):
        seq_dict = self.parseSequenceArg(seqs)
        has_u = False
        has_t = False
        for taxon, seq in list(seq_dict.items()):
            if 'U' in seq:
                has_u = True
            if 'T' in seq:
                has_t = True
        if datatype.upper() == 'DNA':
            self.assertTrue(has_t)
            self.assertFalse(has_u)
        elif datatype.upper() == 'RNA':
            self.assertTrue(has_u)
            self.assertFalse(has_t)


