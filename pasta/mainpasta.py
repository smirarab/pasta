#! /usr/bin/env python

"""Main script of PASTA in command-line mode
"""

# This file is part of PASTA, and is forked from SATe

# PASTA like SATe is free software: you can redistribute it and/or modify
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


from math import ceil
import optparse
from random import sample
import re
import signal

from pasta import PROGRAM_NAME, PROGRAM_VERSION, PROGRAM_LONG_DESCRIPTION, set_timing_log_filepath
from pasta import filemgr
from pasta.alignment import MultiLocusDataset, compact
from pasta.configure import get_configuration
from pasta.pastajob import *
from pasta.scheduler import  stop_worker
from pasta.tools import *
from pasta.treeholder import read_and_encode_splits,\
    generate_tree_with_splits_from_tree
from pasta.utility import IndentedHelpFormatterWithNL


_RunningJobs = None

_LOG = get_logger(__name__)


def fasttree_to_raxml_model_str(datatype, model_str):
    dtu = datatype.upper()
    msu = model_str.upper()
    if dtu == "PROTEIN":
        if "-WAG" in msu:
            if "-GAMMA" in msu:
                return "PROTGAMMAWAGF"
            return "PROTCATWAGF"
        if "-GAMMA" in msu:
            return "PROTGAMMAJTTF"
        return "PROTCATJTTF"
    if "-GAMMA" in msu:
        return "GTRGAMMA"
    return "GTRCAT"
    

def get_auto_defaults_from_summary_stats(datatype, ntax_nchar_tuple_list, total_num_tax):
    """
    Returns nested dictionaries with the following keys set:

    "commandline" : ["multilocus", "datatype"],
    "pasta" : ["tree_estimator",  "aligner", "merger", "break_strategy",
              "move_to_blind_on_worse_score",  "start_tree_search_from_current",
              "after_blind_iter_without_imp_limit", "max_subproblem_size", 
              "max_subproblem_frac", "num_cpus", 
              "time_limit", "after_blind_time_without_imp_limit"],
    "fasttree" : ["model", "GUI_model']


    DO NOT delete keys from this dictionary without making sure that the GUI
        code can deal with your change! This code is used to keep the --auto
        command line option and the data-set dependent defaults of the GUI 
        in sync!
    """
    new_defaults = {}
    new_pasta_defaults = {
        'tree_estimator' : 'fasttree',
        'aligner' : 'mafft',
        'merger' : 'opal',
        'break_strategy' : 'mincluster',
        'move_to_blind_on_worse_score' : True,
        'start_tree_search_from_current' : True,
        'after_blind_iter_without_imp_limit' : -1,
        'time_limit' : -1,
        'blind_after_total_iter': 0,
        'iter_limit' : 3,
        'after_blind_time_without_imp_limit' : -1,
        'mask_gappy_sites' : total_num_tax / 1000 ,
        'build_MST' : False,
        'treeshrink_filter': False
        }
    if total_num_tax > 400:
        new_pasta_defaults['max_subproblem_size'] = 200
        new_pasta_defaults['max_subproblem_frac'] = 0
    else:
        new_pasta_defaults['max_subproblem_size'] = int(math.ceil(total_num_tax/2.0))
        new_pasta_defaults['max_subproblem_frac'] = 0.5
    if datatype.lower() == 'protein':
        new_defaults['fasttree'] = {
            'model' : '-wag -gamma -fastest',
            'GUI_model' : 'WAG+G20'
            }
    else:
        new_defaults['fasttree'] = {
            'model' : '-gtr -gamma -fastest',
            'GUI_model' : 'GTR+G20'
            }
     
    num_cpu = 1
    try:
        import multiprocessing
        num_cpu = multiprocessing.cpu_count()
    except:
        pass
    new_pasta_defaults['num_cpus'] = num_cpu
        
    new_defaults['sate'] = new_pasta_defaults
    new_commandline_defaults = {
        'datatype' : datatype.lower()
        }
    new_commandline_defaults['multilocus'] = False #bool(len(ntax_nchar_tuple_list) > 1)
    new_defaults['commandline'] = new_commandline_defaults
    #_LOG.debug('Auto defaults dictionary: %s' % str(new_defaults))
    return new_defaults

def killed_handler(n, frame):
    global _RunningJobs
    if _RunningJobs:
        MESSENGER.send_warning("signal killed_handler called. Killing running jobs...\n")
        if isinstance(_RunningJobs, list):
            for j in _RunningJobs:
                j.kill()
                MESSENGER.send_warning("kill called...\n")
        else:
            j = _RunningJobs
            j.kill()
            MESSENGER.send_warning("kill called...\n")
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

def finish_pasta_execution(pasta_team,
                          user_config,
                          temporaries_dir,
                          pasta_products,
			  multilocus_dataset):
    global _RunningJobs

    options = user_config.commandline

    user_config.save_to_filepath(os.path.join(temporaries_dir, 'last_used.cfg'))
    if options.timesfile:
        f = open_with_intermediates(options.timesfile, 'a')
        f.close()
        set_timing_log_filepath(options.timesfile)

    ############################################################################
    # Launch threads to do work
    #####
    pasta_config = user_config.get("sate")
    start_worker(pasta_config.num_cpus)
    
    
    #_LOG.debug("start reading the input alignment")
    #multilocus_dataset = read_input_sequences(user_config.input_seq_filepaths,
    #        datatype=user_config.commandline.datatype,
    #        missing=user_config.commandline.missing)
        
    ############################################################################
    # We must read the incoming tree in before we call the get_sequences_for_pasta
    #   function that relabels that taxa in the dataset
    ######
    alignment_as_tmp_filename_to_report = None
    tree_as_tmp_filename_to_report = None
    starting_tree = None
        
    tree_file = options.treefile
    if tree_file:
        if not os.path.exists(tree_file):
            raise Exception('The tree file "%s" does not exist' % tree_file)
        tree_f = open(tree_file, 'rU')
        MESSENGER.send_info('Reading starting trees from "%s"...' % tree_file)
        try:
            tree_list = read_and_encode_splits(multilocus_dataset.dataset, tree_f,
                    starting_tree=True)
        except KeyError:
            MESSENGER.send_error("Error in reading the treefile, probably due to a name in the tree that does not match the names in the input sequence files.\n")
            raise
        except:
            MESSENGER.send_error("Error in reading the treefile.\n")
            raise
        tree_f.close()
        if len(tree_list) > 1:
            MESSENGER.send_warning('%d starting trees found in "%s". The first tree will be used.' % (len(tree_list), tree_file))
        starting_tree = tree_list[0]
        score = None
        tree_as_tmp_filename_to_report = tree_file

    ############################################################################
    # This will relabel the taxa if they have problematic names
    #####
    multilocus_dataset.relabel_for_pasta()

    ############################################################################
    # This ensures all nucleotide data is DNA internally
    #####
    restore_to_rna = False
    if user_config.commandline.datatype.upper() == 'RNA':
        multilocus_dataset.convert_rna_to_dna()
        user_config.commandline.datatype = 'DNA'
        restore_to_rna = True

    export_names = True
    if export_names:
        try:
            name_filename = pasta_products.get_abs_path_for_tag('name_translation.txt')
            name_output = open(name_filename, 'w')
            safe2real = multilocus_dataset.safe_to_real_names
            safe_list = list(safe2real.keys())
            safe_list.sort()
            for safe in safe_list:
                orig = safe2real[safe][0]
                name_output.write("%s\n%s\n\n" % (safe, orig))
            name_output.close()
            MESSENGER.send_info("Name translation information saved to %s as safe name, original name, blank line format." % name_filename)
        except:
            MESSENGER.send_info("Error exporting saving name translation to %s" % name_filename)
            
    
    if options.aligned:
        options.aligned = all( [i.is_aligned() for i in multilocus_dataset] )

    ############################################################################
    # Be prepared to kill any long running jobs
    #####
    prev_signals = []
    for sig in [signal.SIGTERM, signal.SIGABRT, signal.SIGINT]: # signal.SIGABRT, signal.SIGBUS, signal.SIGINT, signal.SIGKILL, signal.SIGSTOP]:
        prev_handler = signal.signal(sig, killed_handler)
        prev_signals.append((sig, prev_handler))

    try:
        pasta_config_dict = pasta_config.dict()
        
        if (not options.two_phase) and tree_file:
            # getting the newick string here will allow us to get a string that is in terms of the correct taxon labels
            starting_tree_str = str(starting_tree)
        else:
            if not options.two_phase:
                MESSENGER.send_info("Creating a starting tree for the PASTA algorithm...")
            if (options.two_phase) or (not options.aligned):
                MESSENGER.send_info("Performing initial alignment of the entire data matrix...")
                init_aln_dir = os.path.join(temporaries_dir, 'init_aln')
                init_aln_dir = pasta_team.temp_fs.create_subdir(init_aln_dir)
                delete_aln_temps = not (options.keeptemp and options.keepalignmenttemps)
                aln_job_list = []
                query_fns = []
                for unaligned_seqs in multilocus_dataset:
                    #backbone = sorted(unaligned_seqs.keys())[0:100]
                    backbone = sample(list(unaligned_seqs.keys()), min(100,len(unaligned_seqs)))   
                    backbone_seqs = unaligned_seqs.sub_alignment(backbone)
                    
                    query_seq=list(set(unaligned_seqs.keys()) - set(backbone))
                    qn = len(query_seq)
                    chunks = min(int(4*pasta_config.num_cpus),int(ceil(qn/50.0)))
                    _LOG.debug("Will align the remaining %d sequences in %d chunks" %(qn,chunks))
                    for ch in range(0,chunks):
                        query_fn = os.path.join(init_aln_dir, "query-%d.fasta"%ch)
                        qa = unaligned_seqs.sub_alignment(query_seq[ch:qn:chunks])
                        _LOG.debug("Chunk with %d sequences built" %len(qa))
                        qa.write_filepath(query_fn)
                        query_fns.append(query_fn)
                    
                    
                    job = pasta_team.aligner.create_job(backbone_seqs,
                                                       tmp_dir_par=init_aln_dir,
                                                       context_str="initalign",
                                                       delete_temps=delete_aln_temps,
						       num_cpus=pasta_config.num_cpus)
                    aln_job_list.append(job)
                _RunningJobs = aln_job_list
                for job in aln_job_list:
                    jobq.put(job)
                
                new_alignment = compact(job.get_results())
                
                add_job_list = []
                for query_fn in query_fns:
                    job = pasta_team.hmmeralign.create_job(new_alignment, query_fn,
                                                        tmp_dir_par=init_aln_dir,
                                                        context_str="initalign",
                                                        delete_temps=delete_aln_temps)
                    add_job_list.append(job)
                _RunningJobs = None
                for job in add_job_list:
                    jobq.put(job)
                for job in add_job_list:
                    new_alignment.merge_in(compact(job.get_results()))
                    #new_alignment_list.apend(new_alignment)
                #for locus_index, new_alignment in enumerate(new_alignment_list):
                multilocus_dataset[0] = new_alignment
                
                if delete_aln_temps:
                    pasta_team.temp_fs.remove_dir(init_aln_dir)
            else:
                MESSENGER.send_info("Input sequences assumed to be aligned (based on sequence lengths).")

            MESSENGER.send_info("Performing initial tree search to get starting tree...")
            init_tree_dir = os.path.join(temporaries_dir, 'init_tree')
            init_tree_dir = pasta_team.temp_fs.create_subdir(init_tree_dir)
            delete_tree_temps = not options.keeptemp
            job = pasta_team.tree_estimator.create_job(multilocus_dataset,
                                                    tmp_dir_par=init_tree_dir,
                                                    num_cpus=pasta_config.num_cpus,
                                                    context_str="inittree",
                                                    delete_temps=delete_tree_temps,
                                                    pasta_products=pasta_products,
                                                    step_num='initialsearch',
                                                    mask_gappy_sites = pasta_config_dict['mask_gappy_sites'])
            _RunningJobs = job
            jobq.put(job)
            score, starting_tree_str = job.get_results()
            _RunningJobs = None
            alignment_as_tmp_filename_to_report = pasta_products.get_abs_path_for_iter_output("initialsearch", TEMP_SEQ_ALIGNMENT_TAG, allow_existing=True)
            tree_as_tmp_filename_to_report = pasta_products.get_abs_path_for_iter_output("initialsearch", TEMP_TREE_TAG, allow_existing=True)
            if delete_tree_temps:
                pasta_team.temp_fs.remove_dir(init_tree_dir)
        _LOG.debug('We have the tree and whole_alignment, partitions...')


        if options.keeptemp:
            pasta_config_dict['keep_iteration_temporaries'] = True
            if options.keepalignmenttemps:
                pasta_config_dict['keep_realignment_temporaries'] = True

        job = PastaJob(multilocus_dataset=multilocus_dataset,
                        pasta_team=pasta_team,
                        name=options.job,
                        status_messages=MESSENGER.send_info,
                        score=score,
                        **pasta_config_dict)
        if starting_tree is not None:            
            job.tree = generate_tree_with_splits_from_tree(starting_tree, force_fully_resolved = True)
        else:
            job.tree_str = starting_tree_str
        job.curr_iter_align_tmp_filename = alignment_as_tmp_filename_to_report
        job.curr_iter_tree_tmp_filename = tree_as_tmp_filename_to_report
        if score is not None:
            job.store_optimum_results(new_multilocus_dataset=multilocus_dataset,
                    new_tree_str=starting_tree_str,
                    new_score=score,
                    curr_timestamp=time.time())

        if options.two_phase:
            MESSENGER.send_info("Exiting with the initial tree because the PASTA algorithm is avoided when the --two-phase option is used.")
        else:
            _RunningJobs = job
            MESSENGER.send_info("Starting PASTA algorithm on initial tree...")
            job.run(tmp_dir_par=temporaries_dir, pasta_products=pasta_products)
            _RunningJobs = None

            if job.return_final_tree_and_alignment:
                alignment_as_tmp_filename_to_report = job.curr_iter_align_tmp_filename
            else:
                alignment_as_tmp_filename_to_report = job.best_alignment_tmp_filename
            
            if user_config.commandline.raxml_search_after:
                raxml_model = user_config.raxml.model.strip()
                if not raxml_model:
                    dt = user_config.commandline.datatype
                    mf = pasta_team.tree_estimator.model
                    ms =  fasttree_to_raxml_model_str(dt, mf)
                    pasta_team.raxml_tree_estimator.model = ms
                rte = pasta_team.raxml_tree_estimator
                MESSENGER.send_info("Performing post-processing tree search in RAxML...")
                post_tree_dir = os.path.join(temporaries_dir, 'post_tree')
                post_tree_dir = pasta_team.temp_fs.create_subdir(post_tree_dir)
                delete_tree_temps = not options.keeptemp
                starting_tree = None
                if user_config.sate.start_tree_search_from_current:
                    starting_tree = job.tree
                post_job = rte.create_job(job.multilocus_dataset,
                                    starting_tree=starting_tree,
                                    num_cpus=pasta_config.num_cpus,
                                    context_str="postraxtree",
                                    tmp_dir_par=post_tree_dir,
                                    delete_temps=delete_tree_temps,
                                    pasta_products=pasta_products,
                                    step_num="postraxtree",
                                    mask_gappy_sites = pasta_config_dict['mask_gappy_sites'])
                _RunningJobs = post_job
                jobq.put(post_job)
                post_score, post_tree = post_job.get_results()
                _RunningJobs = None
                tree_as_tmp_filename_to_report = pasta_products.get_abs_path_for_iter_output("postraxtree", TEMP_TREE_TAG, allow_existing=True)
                if delete_tree_temps:
                    pasta_team.temp_fs.remove_dir(post_tree_dir)
                job.tree_str = post_tree
                job.score = post_score
                if post_score > job.best_score:
                    job.best_tree_str = post_tree
                    job.best_score = post_score
            else:
                if job.return_final_tree_and_alignment:
                    tree_as_tmp_filename_to_report = job.curr_iter_tree_tmp_filename
                else:
                    tree_as_tmp_filename_to_report = job.best_tree_tmp_filename


        #######################################################################
        # Restore original taxon names and RNA characters
        #####
        job.multilocus_dataset.restore_taxon_names()
        if restore_to_rna:
            job.multilocus_dataset.convert_dna_to_rna()
            user_config.commandline.datatype = 'RNA'

        assert len(pasta_products.alignment_streams) == len(job.multilocus_dataset)
        for i, alignment in enumerate(job.multilocus_dataset):
            alignment_stream = pasta_products.alignment_streams[i]
            MESSENGER.send_info("Writing resulting alignment to %s" % alignment_stream.name)
            alignment.write(alignment_stream, file_format="FASTA")
            alignment_stream.close()


        MESSENGER.send_info("Writing resulting tree to %s" % pasta_products.tree_stream.name)
        tree_str = job.tree.compose_newick()
        pasta_products.tree_stream.write("%s;\n" % tree_str)
        pasta_products.tree_stream.close()


        #outtree_fn = options.result
        #if outtree_fn is None:
        #        outtree_fn = os.path.join(seqdir, "combined_%s.tre" % options.job)
        #    else:
        #        outtree_fn = aln_filename + ".tre"
        #MESSENGER.send_info("Writing resulting tree to %s" % outtree_fn)
        #tree_str = str(job.tree)
        #pasta_products.tree_stream.write("%s;\n" % tree_str)


        MESSENGER.send_info("Writing resulting likelihood score to %s" % pasta_products.score_stream.name)
        pasta_products.score_stream.write("%s\n" % job.score)
        pasta_products.score_stream.close()
        
        if alignment_as_tmp_filename_to_report is not None:
            MESSENGER.send_info('The resulting alignment (with the names in a "safe" form) was first written as the file "%s"' % alignment_as_tmp_filename_to_report)
        if tree_as_tmp_filename_to_report is not None:
            MESSENGER.send_info('The resulting tree (with the names in a "safe" form) was first written as the file "%s"' % tree_as_tmp_filename_to_report)

    finally:      
        stop_worker()  
        for el in prev_signals:
            sig, prev_handler = el
            if prev_handler is None:
                signal.signal(sig, signal.SIG_DFL)
            else:
                signal.signal(sig, prev_handler)

def run_pasta_from_config(user_config, pasta_products, multilocus_dataset):
    """
    Returns (None, None) if no temporary directory is left over from the run
    or returns (dir, temp_fs) where `dir` is the path to the temporary
    directory created for the scratch files and `temp_fs` is the TempFS
    instance used to create `dir`
    """

    cmdline_options = user_config.commandline

    ############################################################################
    # Create the safe directory for temporaries
    # The general form of the directory is
    #   ${options.temporaries}/${options.job}/temp${RANDOM}
    ######
    par_dir = cmdline_options.temporaries
    if par_dir is None:
        par_dir = os.path.join(os.path.expanduser('~'), '.pasta')
    cmdline_options.job = coerce_string_to_nice_outfilename(cmdline_options.job, "Job", "pastajob")
    subdir = cmdline_options.job
    par_dir = os.path.abspath(os.path.join(par_dir, subdir))
    if not os.path.exists(par_dir):
        os.makedirs(par_dir) # this parent directory will not be deleted, so we don't store it in the pasta_team.temp_fs

    pasta_team = PastaTeam(config=user_config)

    delete_dir = not cmdline_options.keeptemp

    temporaries_dir = pasta_team.temp_fs.create_top_level_temp(parent=par_dir, prefix='temp')
    assert(os.path.exists(temporaries_dir))
    try:
        MESSENGER.send_info("Directory for temporary files created at %s" % temporaries_dir)
        finish_pasta_execution(pasta_team=pasta_team,
                              user_config=user_config,
                              temporaries_dir=temporaries_dir,
                              pasta_products=pasta_products,
                              multilocus_dataset=multilocus_dataset)
#    except:
#        stop_worker()
#        raise
    finally:
        if delete_dir:
            pasta_team.temp_fs.remove_dir(temporaries_dir)
    if delete_dir:
        return None, None
    else:
        return temporaries_dir, pasta_team.temp_fs

def coerce_string_to_nice_outfilename(p, reason, default):
    illegal_filename_pattern = re.compile(r'[^-_a-zA-Z0-9.]')
    j = "".join(illegal_filename_pattern.split(p))
    if not j:
        j = default
    if j != p:
        MESSENGER.send_warning('%s name changed from "%s" to "%s" (a safer name for filepath)' % (reason, p, j))
    return j


def populate_auto_options(user_config, md, force=False):
    if user_config.commandline.input is None:
        sys.exit("ERROR: Input file(s) not specified.")
    from pasta.usersettingclasses import get_list_of_seq_filepaths_from_dir
    from pasta.alignment import summary_stats_from_parse
    try:
        if user_config.commandline.multilocus:
            fn_list = get_list_of_seq_filepaths_from_dir(user_config.commandline.input)
        else:
            fn_list = [user_config.commandline.input]
        datatype_list = [user_config.commandline.datatype.upper()]
        careful_parse = user_config.commandline.untrusted
        summary_stats = summary_stats_from_parse(fn_list, datatype_list, md, careful_parse=careful_parse)
    except:
        if user_config.commandline.auto:
            MESSENGER.send_error("Error reading input while setting options for the --auto mode\n")
        else:
            MESSENGER.send_error("Error reading input\n")
        raise
    if force or user_config.commandline.auto:
        user_config.commandline.auto = False
        auto_opts = get_auto_defaults_from_summary_stats(summary_stats[0], summary_stats[1], summary_stats[2])
        user_config.get('sate').set_values_from_dict(auto_opts['sate'])
        user_config.get('commandline').set_values_from_dict(auto_opts['commandline'])
        user_config.get('fasttree').set_values_from_dict(auto_opts['fasttree'])
    return summary_stats


def parse_user_options(argv, parser, user_config, command_line_group):
    if argv == sys.argv:
        options, args = parser.parse_args(argv[1:])
    else:
        options, args = parser.parse_args(argv)

    if options.multilocus:
        MESSENGER.send_error(''' Multilocus mode is not supported by PASTA. 
It's a legacy option inherited from SATe.''')
        sys.exit(1)
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
    command_line_group.job = coerce_string_to_nice_outfilename(command_line_group.job, 'Job', 'pastajob')

def check_user_options(user_config):
    if user_config.sate.max_subproblem_size == 1:
        MESSENGER.send_error(''' You have specified a max subproblem size of 1.
PASTA requires a max subproblem size of at least 2.  ''')
        sys.exit(1)

    

def pasta_main(argv=sys.argv):
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
    
    # This is just to read the configurations so that auto value could be set
    parse_user_options(argv, parser, user_config, command_line_group)
    
    # Read the input file, this is needed for auto values
    user_config.read_seq_filepaths(src=user_config.commandline.input,
            multilocus=user_config.commandline.multilocus)
    multilocus_dataset = read_input_sequences(user_config.input_seq_filepaths,
            datatype=user_config.commandline.datatype,
            missing=user_config.commandline.missing)

    # This is to automatically set the auto default options
    summary_stats = populate_auto_options(user_config, multilocus_dataset, force = True)

    # This is to actually read the config files and commandline args and overwrite auto value
    parse_user_options(argv, parser, user_config, command_line_group)
        
    # This is now to make sure --auto overwrites user options
    if user_config.commandline.auto or (user_config.commandline.untrusted):
        summary_stats = populate_auto_options(user_config, multilocus_dataset)
            
    check_user_options(user_config)

    if user_config.sate.mask_gappy_sites < 1:
        user_config.sate.mask_gappy_sites = int(summary_stats[2] * user_config.sate.mask_gappy_sites)
    MESSENGER.send_info("Masking alignment sites with less than %d sites before running the tree step" %user_config.sate.mask_gappy_sites)

    if user_config.commandline.raxml_search_after:
        if user_config.sate.tree_estimator.upper() != 'FASTTREE':
            sys.exit("ERROR: the 'raxml_search_after' option is only supported when the tree_estimator is FastTree")

    exportconfig = command_line_group.exportconfig
    if exportconfig:
        command_line_group.exportconfig = None
        user_config.save_to_filepath(exportconfig)

        ### TODO: wrap up in messaging system
        sys.stdout.write('Configuration written to "%s". Exiting successfully.\n' % exportconfig )

        return True, None, None

    if user_config.commandline.input is None:
        sys.exit("ERROR: Input file(s) not specified.")

    # note: need to read sequence files first to allow PastaProducts to
    # correctly self-configure
    pasta_products = filemgr.PastaProducts(user_config)
    
    export_config_as_temp = True
    if export_config_as_temp:
        name_cfg = pasta_products.get_abs_path_for_tag('pasta_config.txt')
        command_line_group.exportconfig = None
        user_config.save_to_filepath(name_cfg)
        MESSENGER.send_info('Configuration written to "%s".\n' % name_cfg )
         

    MESSENGER.run_log_streams.append(pasta_products.run_log_stream)
    MESSENGER.err_log_streams.append(pasta_products.err_log_stream)
    temp_dir, temp_fs = run_pasta_from_config(user_config, pasta_products, multilocus_dataset)
    _TIME_SPENT = time.time() - _START_TIME
    MESSENGER.send_info("Total time spent: %ss" % _TIME_SPENT)
    return True, temp_dir, temp_fs
