#! /usr/bin/env python

"""Main script of SATe in command-line mode
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
import re
import sys
import signal
import time
import glob
import optparse
import sate

from sate import PROGRAM_NAME, PROGRAM_VERSION, PROGRAM_LONG_DESCRIPTION, get_logger, set_timing_log_filepath, TIMING_LOG, MESSENGER
from sate.alignment import Alignment, SequenceDataset, MultiLocusDataset
from sate.configure import get_configuration, get_input_source_directory
from sate.tree import PhylogeneticTree
from sate.tools import *
from sate.satejob import *
from sate.treeholder import read_and_encode_splits
from sate.scheduler import start_worker, jobq
from sate.utility import IndentedHelpFormatterWithNL
from sate.filemgr import open_with_intermediates
from sate import filemgr

_RunningJobs = None

_LOG = get_logger(__name__)

def killed_handler(n, frame):
    global _RunningJobs
    if _RunningJobs:
        MESSENGER.send_warning("signal killed_handler called. Killing running jobs...\n")
        j = _RunningJobs
        j.kill()
    else:
        MESSENGER.send_warning("signal killed_handler called with no jobs running. Exiting.\n")
    sys.exit()

def read_input_sequences(seq_filename_list,
        datatype,
        missing=None):
    md = MultiLocusDataset()
    md.read_files(seq_filename_list=seq_filename_list,
            datatype=datatype,
            missing=missing)
    return md

def finish_sate_execution(sate_team,
                          user_config,
                          temporaries_dir,
                          multilocus_dataset,
                          sate_products):
    global _RunningJobs
    # get the RAxML model #TODO: this should check for the tree_estimator.  Currently we only support raxml, so this works...
    model = user_config.raxml.model

    options = user_config.commandline

    user_config.save_to_filepath(os.path.join(temporaries_dir, 'last_used.cfg'))
    if options.timesfile:
        f = open_with_intermediates(options.timesfile, 'a')
        f.close()
        set_timing_log_filepath(options.timesfile)
    ############################################################################
    # We must read the incoming tree in before we call the get_sequences_for_sate
    #   function that relabels that taxa in the dataset
    ######
    tree_file = options.treefile
    if tree_file:
        if not os.path.exists(tree_file):
            raise Exception('The tree file "%s" does not exist' % tree_file)
        tree_f = open(tree_file, 'rU')
        MESSENGER.send_info('Reading starting trees from "%s"...' % tree_file)
        tree_list = read_and_encode_splits(multilocus_dataset.dataset, tree_f)
        tree_f.close()
        if len(tree_list) > 1:
            MESSENGER.send_warning('%d starting trees found in "%s". The first tree will be used.' % (len(tree_list), tree_file))
        starting_tree = tree_list[0]
        score = None

    ############################################################################
    # This will relabel the taxa if they have problematic names
    #####
    multilocus_dataset.relabel_for_sate()

    options.aligned = all( [i.is_aligned() for i in multilocus_dataset] )

    ############################################################################
    # Launch threads to do work
    #####
    sate_config = user_config.get("sate")
    start_worker(sate_config.num_cpus)

    ############################################################################
    # Be prepared to kill any long running jobs
    #####
    prev_signals = []
    for sig in [signal.SIGTERM, signal.SIGABRT, signal.SIGINT]: # signal.SIGABRT, signal.SIGBUS, signal.SIGINT, signal.SIGKILL, signal.SIGSTOP]:
        prev_handler = signal.signal(sig, killed_handler)
        prev_signals.append((sig, prev_handler))

    try:
        if tree_file:
            # getting the newick string here will allow us to get a string that is in terms of the correct taxon labels
            starting_tree_str = starting_tree.compose_newick()
        else:
            MESSENGER.send_info("Performing initial tree search to get starting tree...")
            if not options.aligned:
                MESSENGER.send_info("Performing initial alignment of the entire data matrix...")
                init_aln_dir = os.path.join(temporaries_dir, 'init_aln')
                init_aln_dir = sate_team.temp_fs.create_subdir(init_aln_dir)
                delete_aln_temps = not (options.keeptemp and options.keepalignmenttemps)
                new_alignment_list= []
                for unaligned_seqs in multilocus_dataset:
                    job = sate_team.aligner.create_job(unaligned_seqs,
                                                       tmp_dir_par=init_aln_dir,
                                                       context_str="initalign",
                                                       delete_temps=delete_aln_temps)
                    _RunningJobs = job
                    jobq.put(job)
                    new_alignment = job.get_results()
                    _RunningJobs = None
                    new_alignment_list.append(new_alignment)
                for locus_index, new_alignment in enumerate(new_alignment_list):
                    multilocus_dataset[locus_index] = new_alignment
                if delete_aln_temps:
                    sate_team.temp_fs.remove_dir(init_aln_dir)
            else:
                MESSENGER.send_info("Input sequences assumed to be aligned (based on sequence lengths).")

            MESSENGER.send_info("Performing initial tree search to get starting tree...")
            init_tree_dir = os.path.join(temporaries_dir, 'init_tree')
            init_tree_dir = sate_team.temp_fs.create_subdir(init_tree_dir)
            delete_tree_temps = not options.keeptemp
            job = sate_team.tree_estimator.create_job(multilocus_dataset,
                                                    tmp_dir_par=init_tree_dir,
                                                    num_cpus=sate_config.num_cpus,
                                                    context_str="inittree",
                                                    delete_temps=delete_tree_temps)
            _RunningJobs = job
            jobq.put(job)
            score, starting_tree_str = job.get_results()
            _RunningJobs = None
            if delete_tree_temps:
                sate_team.temp_fs.remove_dir(init_tree_dir)
        _LOG.debug('We have the tree and whole_alignment, partitions...')

        sate_config_dict = sate_config.dict()

        if options.keeptemp:
            sate_config_dict['keep_iteration_temporaries'] = True
            if options.keepalignmenttemps:
                sate_config_dict['keep_realignment_temporaries'] = True

        job = SateJob(multilocus_dataset=multilocus_dataset,
                        sate_team=sate_team,
                        name=options.job,
                        status_messages=MESSENGER.send_info,
                        **sate_config_dict)
        job.tree_str = starting_tree_str
        if score is not None:
            job.store_optimum_results(new_multilocus_dataset=multilocus_dataset,
                    new_tree_str=starting_tree_str,
                    new_score=score,
                    curr_timestamp=time.time())

        _RunningJobs = job
        MESSENGER.send_info("Starting SATe algorithm on initial tree...")
        job.run(tmp_dir_par=temporaries_dir)
        _RunningJobs = None
        job.multilocus_dataset.restore_taxon_names()
        assert len(sate_products.alignment_streams) == len(job.multilocus_dataset)
        for i, alignment in enumerate(job.multilocus_dataset):
            alignment_stream = sate_products.alignment_streams[i]
            MESSENGER.send_info("Writing final alignment to %s" % alignment_stream.name)
            alignment.write(alignment_stream, file_format="FASTA")
            alignment_stream.close()


        MESSENGER.send_info("Writing final tree to %s" % sate_products.tree_stream.name)
        tree_str = job.tree.compose_newick()
        sate_products.tree_stream.write("%s;\n" % tree_str)


        #outtree_fn = options.result
        #if outtree_fn is None:
        #    if options.multilocus:
        #        outtree_fn = os.path.join(seqdir, "combined_%s.tre" % options.job)
        #    else:
        #        outtree_fn = aln_filename + ".tre"
        #MESSENGER.send_info("Writing final tree to %s" % outtree_fn)
        #tree_str = job.tree.compose_newick()
        #sate_products.tree_stream.write("%s;\n" % tree_str)


        MESSENGER.send_info("Writing final likelihood score to %s" % sate_products.score_stream.name)
        sate_products.score_stream.write("%s\n" % job.score)
    finally:
        for el in prev_signals:
            sig, prev_handler = el
            if prev_handler is None:
                signal.signal(sig, signal.SIG_DFL)
            else:
                signal.signal(sig, prev_handler)

def run_sate_from_config(user_config, sate_products):
    """
    Returns (None, None) if no temporary directory is left over from the run
    or returns (dir, temp_fs) where `dir` is the path to the temporary
    directory created for the scratch files and `temp_fs` is the TempFS
    instance used to create `dir`
    """

    multilocus_dataset = read_input_sequences(user_config.input_seq_filepaths,
            datatype=user_config.commandline.datatype,
            missing=user_config.commandline.missing)
    cmdline_options = user_config.commandline

    ############################################################################
    # Create the safe directory for temporaries
    # The general form of the directory is
    #   ${options.temporaries}/${options.job}/temp${RANDOM}
    ######
    par_dir = cmdline_options.temporaries
    if par_dir is None:
        par_dir = os.path.join(os.path.expanduser('~'), '.sate')
    cmdline_options.job = coerce_string_to_nice_outfilename(cmdline_options.job, "Job", "satejob")
    subdir = cmdline_options.job
    par_dir = os.path.abspath(os.path.join(par_dir, subdir))
    if not os.path.exists(par_dir):
        os.makedirs(par_dir) # this parent directory will not be deleted, so we don't store it in the sate_team.temp_fs

    sate_team = SateTeam(config=user_config)

    delete_dir = not cmdline_options.keeptemp

    temporaries_dir = sate_team.temp_fs.create_top_level_temp(parent=par_dir, prefix='temp')
    assert(os.path.exists(temporaries_dir))
    try:
        MESSENGER.send_info("Directory for temporary files created at %s" % temporaries_dir)
        finish_sate_execution(sate_team=sate_team,
                              user_config=user_config,
                              temporaries_dir=temporaries_dir,
                              multilocus_dataset=multilocus_dataset,
                              sate_products=sate_products)
    finally:
        if delete_dir:
            sate_team.temp_fs.remove_dir(temporaries_dir)
    if delete_dir:
        return None, None
    else:
        return temporaries_dir, sate_team.temp_fs

def coerce_string_to_nice_outfilename(p, reason, default):
    illegal_filename_pattern = re.compile(r'[^-_a-zA-Z0-9.]')
    j = "".join(illegal_filename_pattern.split(p))
    if not j:
        j = default
    if j != p:
        MESSENGER.send_warning('%s name changed from "%s" to "%s" (a safer name for filepath)' % (reason, p, j))
    return j

def sate_main(argv=sys.argv):
    '''Returns (True, dir, temp_fs) on successful execution or raises an exception.

    Where `dir` is either None or the undeleted directory of temporary files.
    and `temp_fs` is is the TempFS object used to create `dir` (if `dir` is
    not None)

    Note that if `argv` is sys.argv then the first element will be skipped, but
        if it is not the sys.argv list then the first element will be interpretted
        as an argument (and will NOT be skipped).
    '''

    _START_TIME = time.time()
    usage = """usage: %prog [options] <settings_file1> <settings_file2> ..."""
    parser = optparse.OptionParser(usage=usage,
                                    description=PROGRAM_LONG_DESCRIPTION,
                                    formatter=IndentedHelpFormatterWithNL(),
                                    version="%s v%s" % (PROGRAM_NAME, PROGRAM_VERSION))

    user_config = get_configuration()
    command_line_group = user_config.get('commandline')
    command_line_group.add_to_optparser(parser)
    sate_group = user_config.get('sate')
    sate_group.add_to_optparser(parser)
    if argv == sys.argv:
        (options, args) = parser.parse_args(argv[1:])
    else:
        (options, args) = parser.parse_args(argv)
    #if options.multilocus:
    #    sys.exit("SATe: Multilocus mode is disabled in this release.")
    config_filenames = list(args)
    for fn in config_filenames:
        if fn[0] == '"' and fn[-1] == '"':
            fn = fn[1:-1]
        if not os.path.exists(fn):
            raise Exception('The configuration (settings) file "%s" does not exist' % fn)
        try:
            user_config.read_config_filepath(fn)
        except:
            raise Exception('The file "%s" does not appear to be a valid configuration file format. It lacks section headers.' % fn)
    user_config.set_values_from_dict(options.__dict__)
    command_line_group.job = coerce_string_to_nice_outfilename(command_line_group.job, 'Job', 'satejob')

    exportconfig = command_line_group.exportconfig
    if exportconfig:
        command_line_group.exportconfig = None
        user_config.save_to_filepath(exportconfig)

        ### TODO: wrap up in messaging system
        sys.stdout.write('Configuration written to "%s". Exiting successfully.' % exportconfig )

        return True, None, None

    if user_config.commandline.input is None:
        sys.exit("ERROR: Input file(s) not specified.")

    # note: need to read sequence files first to allow SateProducts to
    # correctly self-configure
    user_config.read_seq_filepaths(src=user_config.commandline.input,
            multilocus=user_config.commandline.multilocus)
    sate_products = filemgr.SateProducts(user_config)

    MESSENGER.run_log_streams.append(sate_products.run_log_stream)
    MESSENGER.err_log_streams.append(sate_products.err_log_stream)
    temp_dir, temp_fs = run_sate_from_config(user_config, sate_products)
    _TIME_SPENT = time.time() - _START_TIME
    MESSENGER.send_info("Total time spent: %ss" % _TIME_SPENT)
    return True, temp_dir, temp_fs
