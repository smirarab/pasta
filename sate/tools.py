#!/usr/bin/env python

"""Interface to external tools (alignment and tree inference programs)
"""

# This file is part of SATe

# SATe is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Jiaye Yu and Mark Holder, University of Kansas

import os
import random
import sys
import time
import platform

from alignment import Alignment
from sate import get_logger, GLOBAL_DEBUG, SATE_SYSTEM_PATHS_CFGFILE, DEFAULT_MAX_MB
from sate.filemgr import open_with_intermediates
from sate.scheduler import jobq, start_worker, DispatchableJob, FakeJob

from cStringIO import StringIO
from alignment import Alignment

from tree import PhylogeneticTree

_LOG = get_logger(__name__)

def is_file_checker(p):
    if not p:
        return False, "Expecting the path to an executable, got an empty string"
    if not os.path.isfile(p):
        return False, "%s is not a valid file" % str(p)
    return True, p

def is_executable_checker(p):
    r, msg = is_file_checker(p)
    if not r:
        return (r, msg)
    try:
        import stat
    except:
        return True, p
    if (os.stat(p)[stat.ST_MODE] & (stat.S_IXUSR|stat.S_IXGRP|stat.S_IXOTH)) == 0:
        return False, "%s does not appear to be an executable file" % p
    else:
        return True, p

def read_internal_alignment(fn,
                            file_format='FASTA',
                            datatype=None,
                            dirs_to_delete=(),
                            temp_fs=None):
    alignment = Alignment()
    alignment.datatype = datatype
    alignment.read_filepath(fn, file_format=file_format)
    if len(alignment) >= 1:
        if dirs_to_delete:
            assert(temp_fs)
            for d in dirs_to_delete:
                time.sleep(1) #TODO: not sure why this is here!
                temp_fs.remove_dir(d)
        return alignment
    else:
        raise ValueError("The alignment file has no sequences. SATe quits." % fn)

def read_raxml_results(dir, dirs_to_delete, temp_fs):
    flist = os.listdir(dir)
    id = None
    for f in flist:
        if f.startswith('RAxML_log'):
            id = f.split('.')[1]
            break
    raxml_log = os.path.join(dir, 'RAxML_log.%s' % id)
    logsc = filter(lambda x : x.find("Final GAMMA-based Score of best tree")!=-1 , open(raxml_log, 'rU').readlines())
    if logsc:
        score = float(logsc[0].split(" ")[-1])
    else:
        score = float(open(raxml_log, 'rU').readlines()[-1].split()[1])
    raxml_result = os.path.join(dir, 'RAxML_result.%s' % id)
    tree_str = open(raxml_result, 'rU').read().strip()
    for d in dirs_to_delete:
        temp_fs.remove_dir(d)
    return score, tree_str

def read_fasttree_results(toclose, dir, fasttree_restults_file, log, delete_dir=False):
        toclose.close()
        tree_str = open(fasttree_restults_file, 'rU').read().strip()
        score = None
        for line in reversed(open(log, 'rU').readlines()):
            if (line.split()[0] == 'Gamma20LogLk'):
                score = float(line.split()[1])
                break                                                                      
            if (line.split()[0] == 'TreeLogLk'):
                score = float(line.split()[2])
                break
        if score is None:
            message = "FastTree did not report a log-likelhood for the data: the data set is probably too weird"
            raise Exception(message)
        return score, tree_str

class ExternalTool (object):

    is_bundled_tool = False

    def __init__(self, name, temp_fs, **kwargs):
        self.name = name
        self.temp_fs = temp_fs
        self.exe = kwargs['path']
        if not os.path.exists(self.exe):
            if hasattr(self, "is_bundled_tool") and self.is_bundled_tool:
                msg = "The path '%s' does not exist. Please check the installation and try again." % self.exe
            else:
                msg = """'%s' not found. Please install '%s' and/or configure its location correctly in '%s'""" % (self.exe, self.name, SATE_SYSTEM_PATHS_CFGFILE)
            raise ValueError(msg)

        self.delete_temps = kwargs.get('delete_temps', True)

    @staticmethod
    def exists(self):
        return is_executable_checker(self.exe)[0]

    def make_temp_workdir(self, tmp_dir_par):
        pref = 'temp' + self.name
        scratch_dir = self.temp_fs.create_temp_subdir(parent=tmp_dir_par, prefix=pref)
        return scratch_dir

    def run(self, *args, **kwargs):
        start_worker(1)
        job = self.create_job(*args, **kwargs)
        jobq.put(job)
        return job.get_results()

    def create_job(self, *args, **kwargs):
        raise NotImplementedError('Abstract ExternalTool class cannot spawn jobs.')

class Aligner(ExternalTool):
    def __init__(self, name, temp_fs, **kwargs):
        ExternalTool.__init__(self, name, temp_fs, **kwargs)
        self.user_opts = kwargs.get('args', ' ').split()

    def _prepare_input(self, alignment, **kwargs):
        """Wraps up the writing of raw fasta, creation of temp dir, ... for common aligners.
        Returns directory, input filename, output filename."""
        tdp = kwargs.get('tmp_dir_par')
        if not tdp:
            raise AssertionError('The tmp_dir_par must be specified when calling create_job or _prepare_input')
        scratch_dir = self.make_temp_workdir(tmp_dir_par=tdp)
        seqfn = os.path.join(scratch_dir, "input.fasta")
        alignment.write_unaligned_fasta(seqfn)
        alignedfn = os.path.join(scratch_dir, 'input.aligned')
        return scratch_dir, seqfn, alignedfn

    def create_job(self, *args, **kwargs):
        raise NotImplementedError('Abstract Aligner class cannot spawn jobs.')

    def _finish_standard_job(self, alignedfn, datatype, invoc, scratch_dir, job_id, delete_temps):
        dirs_to_delete = []
        if delete_temps:
            dirs_to_delete = [scratch_dir]
        # create a results processor to read the alignment file
        rpc = lambda : read_internal_alignment(alignedfn,
                                               datatype=datatype,
                                               dirs_to_delete=dirs_to_delete,
                                               temp_fs=self.temp_fs)
        job = DispatchableJob(invoc, result_processor=rpc,  cwd=scratch_dir, context_str=job_id)
        return job

class CustomAligner(Aligner):
    section_name = 'custom aligner'

    def __init__(self, name, temp_fs, **kwargs):
        Aligner.__init__(self, name, temp_fs, **kwargs)

    def create_job(self, alignment, guide_tree=None):
        raise NotImplementedError('User-provided Aligner NOT supported yet.')


class MafftAligner(Aligner):
    section_name = 'mafft aligner'
    url = 'http://align.bmr.kyushu-u.ac.jp/mafft/software'
    is_bundled = True

    def __init__(self, temp_fs, **kwargs):
        Aligner.__init__(self, 'mafft', temp_fs, **kwargs)

    def create_job(self, alignment, guide_tree=None, **kwargs):
        job_id = kwargs.get('context_str', '') + '_mafft'
        if alignment.get_num_taxa() == 1:
            return FakeJob(alignment, context_str=job_id)

        scratch_dir, seqfn, alignedfn = self._prepare_input(alignment, **kwargs)
        aligned_fileobj = open_with_intermediates(alignedfn, 'w')

        invoc = []
        if platform.system() == "Windows":
            invoc.append(self.exe)
        else:
            invoc.extend([sys.executable, self.exe])
        if len(alignment) <= 200 and alignment.max_sequence_length() < 10000:
            invoc.extend(['--localpair', '--maxiterate', '1000'])
        if '--ep' not in self.user_opts:
            invoc.extend(['--ep', '0.123'])
        invoc.extend(['--quiet', seqfn])
        invoc.extend(self.user_opts)

        # The MAFFT job creation is slightly different from the other
        #   aligners because we redirect and read standard output.

        dirs_to_delete = []
        if kwargs.get('delete_temps', self.delete_temps):
            dirs_to_delete.append(scratch_dir)

        def mafft_result_processor(to_close=aligned_fileobj,
                                    fn=alignedfn,
                                    datatype=alignment.datatype,
                                    dirs_to_delete=dirs_to_delete,
                                    temp_fs=self.temp_fs):
            to_close.close()
            return read_internal_alignment(fn=alignedfn,
                                           datatype=datatype,
                                           dirs_to_delete=dirs_to_delete,
                                           temp_fs=temp_fs)

        job = DispatchableJob(invoc,
                              result_processor=mafft_result_processor,
                              cwd=scratch_dir,
                              stdout=aligned_fileobj,
                              context_str=job_id)
        return job


class OpalAligner(Aligner):
    section_name = 'opal aligner'
    url = "http://opal.cs.arizona.edu"
    is_bundled = True

    @staticmethod
    def checker(p, config):
        return is_file_checker(p)

    def __init__(self, temp_fs, **kwargs):
        Aligner.__init__(self, 'opal', temp_fs, **kwargs)
        self.max_mem_mb = kwargs.get("max_mem_mb", DEFAULT_MAX_MB)

    def create_job(self, alignment, guide_tree=None, **kwargs):
        job_id = kwargs.get('context_str', '') + '_opal'
        if alignment.get_num_taxa() == 1:
            return FakeJob(alignment, context_str=job_id)
        scratch_dir, seqfn, alignedfn = self._prepare_input(alignment, **kwargs)

        invoc = ['java', '-Xmx%dm' % self.max_mem_mb, '-jar', self.exe, '--in', seqfn, '--out', alignedfn, '--quiet']
        invoc.extend(self.user_opts)

        return self._finish_standard_job(alignedfn=alignedfn,
                                        datatype=alignment.datatype,
                                        invoc=invoc,
                                        scratch_dir=scratch_dir,
                                        job_id=job_id,
                                        delete_temps=kwargs.get('delete_temps', self.delete_temps))


class Clustalw2Aligner(Aligner):
    section_name = 'clustalw2 aligner'
    url = 'http://www.ebi.ac.uk/Tools/clustalw2/index.html'
    is_bundled = True

    def __init__(self, temp_fs, **kwargs):
        Aligner.__init__(self, 'clustalw2', temp_fs, **kwargs)

    def create_job(self, alignment, guide_tree=None, **kwargs):
        job_id = kwargs.get('context_str', '') + '_clustalw2'
        if alignment.get_num_taxa() == 1:
            return FakeJob(alignment, context_str=job_id)
        scratch_dir, seqfn, alignedfn = self._prepare_input(alignment, **kwargs)

        invoc = [self.exe, '-align', '-infile=%s' % seqfn, '-outfile=%s' % alignedfn, '-output=fasta']
        invoc.extend(self.user_opts)

        return self._finish_standard_job(alignedfn=alignedfn,
                                        datatype=alignment.datatype,
                                        invoc=invoc,
                                        scratch_dir=scratch_dir,
                                        job_id=job_id,
                                        delete_temps=kwargs.get('delete_temps', self.delete_temps))

class ProbalignAligner(Aligner):
    section_name = 'probalign'
    url = 'http://probalign.njit.edu'
    is_bundled_tool = True

    def __init__(self, temp_fs, **kwargs):
        Aligner.__init__(self, 'probalign', temp_fs, **kwargs)

    def create_job(self, alignment, guide_tree=None, **kwargs):
        job_id = kwargs.get('context_str', '') + '_probalign'
        if alignment.get_num_taxa() == 1:
            return FakeJob(alignment, context_str=job_id)
        scratch_dir, seqfn, alignedfn = self._prepare_input(alignment, **kwargs)

        invoc = [self.exe, '-nuc', "-o", alignedfn, seqfn]
        invoc.extend(self.user_opts)

        return self._finish_standard_job(alignedfn=alignedfn,
                                        datatype=alignment.datatype,
                                        invoc=invoc,
                                        scratch_dir=scratch_dir,
                                        job_id=job_id,
                                        delete_temps=kwargs.get('delete_temps', self.delete_temps))


class PrankAligner(Aligner):
    section_name = 'prank aligner'
    url = 'http://www.ebi.ac.uk/goldman-srv/prank/prank'
    is_bundled = True

    def __init__(self, temp_fs, **kwargs):
        Aligner.__init__(self, 'prank', temp_fs, **kwargs)

    def create_job(self, alignment, guide_tree=None, **kwargs):
        job_id = kwargs.get('context_str', '') + '_prank'
        if alignment.get_num_taxa() == 1:
            return FakeJob(alignment, context_str=job_id)
        scratch_dir, seqfn, alignedfn = self._prepare_input(alignment, **kwargs)

        invoc = [self.exe, '-once', '-noxml', '-notree', '-nopost', '+F', '-quiet', '-matinitsize=5', '-uselogs', '-d=%s' % seqfn, '-o=%s' % alignedfn]
        alignedfn = alignedfn + '.1.fas'

        return self._finish_standard_job(alignedfn=alignedfn,
                                        datatype=alignment.datatype,
                                        invoc=invoc,
                                        scratch_dir=scratch_dir,
                                        job_id=job_id,
                                        delete_temps=kwargs.get('delete_temps', self.delete_temps))


class FakeAligner(Aligner):
    "Simply returns the input data -- I hope that it is aligned!"
    section_name = 'fakealigner'
    url = ''

    def __init__(self, temp_fs, **kwargs):
        Aligner.__init__(self, 'fakealigner', temp_fs, **kwargs)

    def create_job(self, alignment, guide_tree=None, **kwargs):
        job_id = kwargs.get('context_str', '') + '_fakealigner'
        return FakeJob(alignment, context_str=job_id)

class PadAligner(Aligner):
    section_name = 'padaligner'
    url = ''
    is_bundled = True

    def __init__(self, temp_fs, **kwargs):
        Aligner.__init__(self, 'padaligner', temp_fs, **kwargs)

    def create_job(self, alignment, guide_tree=None, **kwargs):
        job_id = kwargs.get('context_str', '') + '_padaligner'
        if alignment.get_num_taxa() == 1:
            return FakeJob(alignment, context_str=job_id)
        scratch_dir, seqfn, alignedfn = self._prepare_input(alignment, **kwargs)

        invoc = [sys.executable, self.exe, alignment.datatype, seqfn, alignedfn]

        return self._finish_standard_job(alignedfn=alignedfn,
                                        datatype=alignment.datatype,
                                        invoc=invoc,
                                        scratch_dir=scratch_dir,
                                        job_id=job_id,
                                        delete_temps=kwargs.get('delete_temps', self.delete_temps))



class Merger(ExternalTool):
    def __init__(self, name, temp_fs, **kwargs):
        ExternalTool.__init__(self, name, temp_fs, **kwargs)

    def _prepare_input(self, alignment1, alignment2, **kwargs):
        scratch_dir = self.make_temp_workdir(tmp_dir_par=kwargs['tmp_dir_par'])
        seqfn1 = os.path.join(scratch_dir, "1.fasta")
        seqfn2 = os.path.join(scratch_dir, "2.fasta")
        alignment1.write_filepath(seqfn1, 'FASTA')
        alignment2.write_filepath(seqfn2, 'FASTA')
        outfn = os.path.join(scratch_dir, 'out.fasta')
        return scratch_dir, seqfn1, seqfn2, outfn

    def _finish_standard_job(self, alignedfn, datatype, invoc, scratch_dir, job_id, delete_temps):
        dirs_to_delete = []
        if delete_temps:
            dirs_to_delete = [scratch_dir]
        # create a results processor to read the alignment file
        rpc = lambda : read_internal_alignment(alignedfn,
                                               datatype=datatype,
                                               dirs_to_delete=dirs_to_delete,
                                               temp_fs=self.temp_fs)
        job = DispatchableJob(invoc, result_processor=rpc,  cwd=scratch_dir, context_str=job_id)
        return job


class CustomMerger(Merger):
    section_name = 'custom aligner'
    url = ''

    def __init__(self, name, temp_fs, **kwargs):
        Merger.__init__(self, name, temp_fs, **kwargs)

    def create_job(self, alignment, guide_tree=None, **kwargs):
        raise NotImplementedError('User-provided Merger NOT supported yet.')

class FakeMerger(Merger):
    section_name = 'fakealigner'
    url = ''
    is_bundled = True

    def __init__(self, temp_fs, **kwargs):
        Merger.__init__(self, 'fakealigner', temp_fs, **kwargs)

    def create_job(self, alignment1, alignment2, **kwargs):
        alignment1.update(alignment2)
        job_id = kwargs.get('context_str', '') + '_fakealigner'
        return FakeJob(alignment1, context_str=job_id)


class PadMerger(Merger):
    section_name = 'padaligner'
    url = ''
    is_bundled = True

    def __init__(self, temp_fs, **kwargs):
        Merger.__init__(self, 'padaligner', temp_fs, **kwargs)

    def create_job(self, alignment1, alignment2, **kwargs):
        scratch_dir, seqfn1, seqfn2, outfn = self._prepare_input(alignment1, alignment2, **kwargs)

        invoc = [sys.executable, self.exe, alignment1.datatype, seqfn1, seqfn2, outfn]

        job_id = kwargs.get('context_str', '') + '_padaligner'

        return self._finish_standard_job(alignedfn=outfn,
                                         datatype=alignment1.datatype,
                                         invoc=invoc,
                                         scratch_dir=scratch_dir,
                                         job_id=job_id,
                                         delete_temps=kwargs.get('delete_temps', self.delete_temps))

class MuscleMerger (Merger):
    section_name = 'muscle merger'
    url = "http://www.drive5.com/muscle"
    is_bundled = True

    def __init__(self, temp_fs, **kwargs):
        Merger.__init__(self, 'muscle', temp_fs, **kwargs)

    def create_job(self, alignment1, alignment2, **kwargs):
        scratch_dir, seqfn1, seqfn2, outfn = self._prepare_input(alignment1, alignment2, **kwargs)

        invoc = [self.exe, '-in1', seqfn1, '-in2', seqfn2, '-out', outfn, '-quiet', '-profile']

        job_id = kwargs.get('context_str', '') + '_muscle'

        return self._finish_standard_job(alignedfn=outfn,
                                         datatype=alignment1.datatype,
                                         invoc=invoc,
                                         scratch_dir=scratch_dir,
                                         job_id=job_id,
                                         delete_temps=kwargs.get('delete_temps', self.delete_temps))

class OpalMerger (Merger):
    section_name = "opal merger"
    url = "http://opal.cs.arizona.edu"
    is_bundled = True

    @staticmethod
    def checker(p, config):
        return is_file_checker(p)

    def __init__(self, temp_fs, **kwargs):
        Merger.__init__(self, 'opal', temp_fs, **kwargs)
        self.max_mem_mb = kwargs.get("max_mem_mb", DEFAULT_MAX_MB)

    def create_job(self, alignment1, alignment2, **kwargs):
        scratch_dir, seqfn1, seqfn2, outfn = self._prepare_input(alignment1, alignment2, **kwargs)
        assert(alignment1.datatype == alignment2.datatype)

        invoc = ['java', '-Xmx%dm' % self.max_mem_mb, '-jar', self.exe, '--in', seqfn1, '--in2', seqfn2, '--out', outfn, '--align_method', 'profile']
        job_id = kwargs.get('context_str', '') + '_opal'

        return self._finish_standard_job(alignedfn=outfn,
                                         datatype=alignment1.datatype,
                                         invoc=invoc,
                                         scratch_dir=scratch_dir,
                                         job_id=job_id,
                                         delete_temps=kwargs.get('delete_temps', self.delete_temps))

class TreeEstimator(ExternalTool):
    def __init__(self, name, temp_fs, **kwargs):
        ExternalTool.__init__(self, name, temp_fs, **kwargs)
        self.model = kwargs.get('model')

    def _prepare_input(self, alignment, **kwargs):
        raise NotImplmentedError('Abstract TreeEstimator class!')

    @staticmethod
    def _read_results(fn):
        raise NotImplmentedError('Abstract TreeEstimator class!')

class CustomTreeEstimator(TreeEstimator):
    section_name = 'custom tree_estimator'
    url = ''

    def __init__(self, name, temp_fs, **kwargs):
        TreeEstimator.__init__(self, name, temp_fs, **kwargs)

    def create_job(self, alignment, starting_tree=None, **kwargs):
        raise NotImplementedError('User-provided Merger NOT supported yet.')

class Randtree(TreeEstimator):
    section_name = 'randtree tree_estimator'
    url = 'http://phylo.bio.ku.edu/software/sate-exe'
    is_bundled = True

    def _prepare_input(self, alignment, **kwargs):
        scratch_dir = self.make_temp_workdir(tmp_dir_par=kwargs['tmp_dir_par'])
        seqfn = os.path.join(scratch_dir, "input.fasta")
        alignment.write_filepath(seqfn, 'FASTA')
        score_fn = os.path.join(scratch_dir, 'scorefile')
        return scratch_dir, seqfn, alignment.datatype, score_fn

    def __init__(self, temp_fs, **kwargs):
        TreeEstimator.__init__(self, 'randtree', temp_fs, **kwargs)

    def create_job(self, alignment, starting_tree=None, name='default', **kwargs):
        scratch_dir, seqfn, dt, score_fn = self._prepare_input(alignment, **kwargs)
        invoc = [sys.executable,
                self.exe,
                seqfn,
                dt,
                os.path.join(scratch_dir, 'output.tre'),
                ]
        score_fileobj = open_with_intermediates(score_fn, 'w')

        dirs_to_delete = []
        if kwargs.get('delete_temps', self.delete_temps):
            dirs_to_delete.append(scratch_dir)

        def randtree_result_processor(dir=scratch_dir,
                                      to_close=score_fileobj,
                                      score_fn=score_fn,
                                      fn=os.path.join(scratch_dir, 'output.tre'),
                                      dirs_to_delete=dirs_to_delete,
                                      temp_fs=self.temp_fs):
            to_close.close()
            score = float(open(score_fn, 'rU').read().strip())
            tree_str = open(fn, 'rU').read().strip()
            for d in dirs_to_delete:
                temp_fs.remove_dir(d)
            return (score, tree_str)

        job_id = kwargs.get('context_str', '') + '_randtree'
        job = DispatchableJob(invoc,
                              result_processor=randtree_result_processor,
                              cwd=scratch_dir,
                              context_str=job_id,
                              stdout=score_fileobj)
        return job

class FakeTreeEstimator(TreeEstimator):
    "Must be sent an starting tree.  It simply returns this tree"
    section_name = 'faketree tree_estimator'
    url = ''
    is_bundled = True

    def __init__(self, temp_fs, **kwargs):
        TreeEstimator.__init__(self, 'faketree', temp_fs, **kwargs)

    def create_job(self, alignment, starting_tree=None, name='default', **kwargs):
        assert(starting_tree)
        job_id = kwargs.get('context_str', '') + '_fake'
        if isinstance(starting_tree, str):
            tree_str = starting_tree
        else:
            tree_str = starting_tree.compose_newick()
        score = hash(tree_str)/10000.0
        blob = (score, tree_str)
        return FakeJob(blob, context_str=job_id)

class FastTree(TreeEstimator):
    section_name = 'fasttree tree_estimator'
    url = 'http://www.microbesonline.org/fasttree/'
    is_bundled_tool = False
    user_opts = []

    def __init__(self, **kwargs):
        TreeEstimator.__init__(self, 'fasttree', **kwargs)
        self.options = kwargs.get('options').split()

    def _prepare_input(self, multilocus_dataset, **kwargs):
        curdir = self.make_temp_workdir(tmp_dir_par=kwargs.get('tmp_dir_par'))
        seqfn = os.path.join(curdir, "input.fasta")

        # TODO: @mth: I added this line following the RAxML tool; is it correct?
        alignment, partitions = multilocus_dataset.concatenate_alignments()
        alignment.write_filepath(seqfn, 'FASTA')

        if alignment.datatype == 'DNA':
            datatype = '-nt'
        elif alignment.datatype == 'PROTEIN':
            datatype = ''
        else:
            raise ValueError('Datatype "%s" not recognized by FastTree' % str(alignment.datatype))

        options = self.options if self.options is not None else ''

        return curdir, seqfn, datatype, options

    def create_job(self, alignment, starting_tree=None, **kwargs):
        scratch_dir, seqfn, datatype, options = self._prepare_input(alignment, **kwargs)
        num_cpus = kwargs.get('num_cpus')
        log_file = os.path.join(scratch_dir, 'log');

        invoc = [self.exe, '-quiet']
        if datatype != '':
            invoc.extend([datatype])

        model = self.model  if self.model is not None else ''
        if model != '':
            model = model.split(' ')
            invoc.extend(model)
        if options is not None and len(options) >=1 :
            invoc.extend(options)

        fasttree_result = os.path.join(scratch_dir, 'results')
        results_fileobj = open_with_intermediates(fasttree_result, 'w')

        if starting_tree is not None:
            if isinstance(starting_tree, str):
                tree_str = starting_tree
            else:
                tree_str = starting_tree.compose_newick()
            tree_fn = os.path.join(os.path.abspath(scratch_dir), "start.tre")
            tree_file_obj = open(tree_fn, "w")
            tree_file_obj.write("%s;\n" % tree_str)
            tree_file_obj.close()
            invoc.extend(['-intree', tree_fn])

        invoc.extend(['-log', log_file,    seqfn ])

        if num_cpus > 1:
            os.putenv("OMP_NUM_THREADS", str(num_cpus))
            invoc[0] += 'MP'
#            if platform.system() == 'Windows':
#                x = invoc[0].split('.')
#                x[-2] += 'p'
#                invoc[0] = '.'.join(x)
#            else:
#                invoc[0] += 'p'

        dirs_to_delete = []
        if kwargs.get('delete_temps', self.delete_temps):
            dirs_to_delete.append(scratch_dir)

        rpc = lambda : read_fasttree_results(results_fileobj, scratch_dir, fasttree_result , log_file, delete_dir=kwargs.get('delete_temps', self.delete_temps))
        job_id = kwargs.get('context_str', '') + '_fasttree'
        job = DispatchableJob(invoc, result_processor=rpc, cwd=scratch_dir, stdout=results_fileobj, context_str=job_id)
        return job

class Raxml(TreeEstimator):
    section_name = 'raxml tree_estimator'
    url = 'http://icwww.epfl.ch/~stamatak'

    def _write_partition_filepath(self, parfn, partitions, model):
        # partition --- list of tuples, [("DNA", 1, 30), ("DNA", 31, 60), ("PROTEIN", 61, 100)]
        file_obj = open_with_intermediates(parfn,'w')
        count = 0
        for item in partitions:
            key = ""
            count += 1
            if item[0] == "DNA":
                key = "DNA"
            elif item[0] == "PROTEIN":
                if model.startswith("PROTGAMMA"):
                    key = model[len("PROTGAMMAI"):] if model.startswith("PROTGAMMAI") else model[len("PROTGAMMA"):]
                if model.startswith("PROTCAT"):
                    key = model[len("PROTCATI"):] if model.startswith("PROTCATI") else model[len("PROTCAT"):]
            file_obj.write("%s, p%s=%s-%s\n" % (key, count, item[1], item[2]) )
        file_obj.close()

    def _prepare_input(self, multilocus_dataset, **kwargs):
        scratch_dir = self.make_temp_workdir(tmp_dir_par=kwargs['tmp_dir_par'])
        seqfn = os.path.join(scratch_dir, "input.phy")
        model = self.model
        alignment, partitions = multilocus_dataset.concatenate_alignments()
        alignment.write_filepath(seqfn, 'PHYLIP')

        if alignment.datatype == 'DNA':
            model = self.model  if self.model is not None and self.model != '' else 'GTRCAT'
        elif alignment.datatype == 'PROTEIN':
            model = self.model  if self.model is not None and self.model != '' else 'PROTCATWAGF'
        else:
            raise ValueError('Datatype "%s" not recognized by RAxML' % str(alignment.datatype))
        parfn = os.path.join(scratch_dir, "partition.txt")
        self._write_partition_filepath(parfn, partitions, model)
        return scratch_dir, seqfn, parfn, model

    def __init__(self, temp_fs, **kwargs):
        TreeEstimator.__init__(self, 'raxml', temp_fs, **kwargs)

    def create_job(self, alignment, starting_tree=None, name='default', **kwargs):
        scratch_dir, seqfn, parfn, model = self._prepare_input(alignment, **kwargs)
        num_cpus = kwargs.get('num_cpus')
        invoc = [self.exe,
                '-m', model,
                '-n', name,
                '-q', parfn,
                '-s', seqfn,
                # '-M', # Branch length estimates per partition
                ]
        x = open(parfn).readlines()
        npar = [i.count(',') for i in x]

        # if npar > 1:
        #   invoc.extend('-M')

        if starting_tree is not None:
            if isinstance(starting_tree, str):
                tree_str = starting_tree
            else:
                tree_str = starting_tree.compose_newick()
            tree_fn = os.path.join(os.path.abspath(scratch_dir), "start.tre")
            tree_file_obj = open(tree_fn, "w")
            tree_file_obj.write("%s;\n" % tree_str)
            tree_file_obj.close()
            invoc.extend(['-t', tree_fn])
        if num_cpus > 1:
            invoc.extend(['-T', str(num_cpus)])
            if platform.system() == 'Windows':
                x = invoc[0].split('.')
                x[-2] += 'p'
                invoc[0] = '.'.join(x)
            else:
                invoc[0] += 'p'
        if GLOBAL_DEBUG:
            invoc.extend(['-p', '123456789'])

        dirs_to_delete = []
        if kwargs.get('delete_temps', self.delete_temps):
            dirs_to_delete.append(scratch_dir)

        rpc = lambda : read_raxml_results(scratch_dir,
                                          dirs_to_delete=dirs_to_delete,
                                          temp_fs=self.temp_fs)
        job_id = kwargs.get('context_str', '') + '_raxml'
        job = DispatchableJob(invoc, result_processor=rpc, cwd=scratch_dir, context_str=job_id)
        return job

if GLOBAL_DEBUG:
    AlignerClasses = (ProbalignAligner, Clustalw2Aligner, MafftAligner, PrankAligner, OpalAligner, PadAligner, FakeAligner, CustomAligner)
    MergerClasses = (MuscleMerger, OpalMerger, PadMerger, FakeMerger, CustomMerger)
    TreeEstimatorClasses = (FastTree, Randtree, Raxml, FakeTreeEstimator, CustomTreeEstimator)
else:
    AlignerClasses = (ProbalignAligner, Clustalw2Aligner, MafftAligner, PrankAligner, OpalAligner, CustomAligner)
    MergerClasses = (MuscleMerger, OpalMerger, CustomMerger)
    TreeEstimatorClasses = (Raxml, FastTree, CustomTreeEstimator)

def get_aligner_classes():
    classes = list(AlignerClasses)
    ret = [i for i in classes if not i.section_name.startswith('custom')]
    return ret

def get_merger_classes():
    classes = list(MergerClasses)
    ret = [i for i in classes if not i.section_name.startswith('custom')]
    return ret

def get_tree_estimator_classes():
    classes = list(TreeEstimatorClasses)
    ret = [i for i in classes if not i.section_name.startswith('custom')]
    return ret

def get_external_tool_classes():
    classes = list(AlignerClasses)
    classes.extend(list(MergerClasses))
    classes.extend(list(TreeEstimatorClasses))
    ret = [i for i in classes if not i.section_name.startswith('custom')]
    return ret

