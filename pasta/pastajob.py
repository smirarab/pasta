#!/usr/bin/env python

# This file is part of PASTA and is forked from SATe

# PASTA, like SATe is free software: you can redistribute it and/or modify
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


import math
import os
import random
import time
import sys
import copy
from threading import Lock
from pasta import get_logger, TEMP_SEQ_UNMASKED_ALIGNMENT_TAG
from dendropy.datamodel.taxonmodel import Taxon
from dendropy import Tree
from pasta.tree import PhylogeneticTree
from io import StringIO
_LOG = get_logger(__name__)

from pasta.treeholder import TreeHolder, resolve_polytomies,\
    read_newick_with_translate
from pasta.pastaalignerjob import PASTAAlignerJob, PASTAMergerJob
from pasta import get_logger
from pasta.utility import record_timestamp
from pasta.scheduler import jobq
from pasta.filemgr import  TempFS
from pasta import TEMP_SEQ_ALIGNMENT_TAG, TEMP_TREE_TAG, MESSENGER, TEMP_SHRUNK_ALIGNMENT_TAG, TEMP_SHRUNK_TREE_TAG

# uym2 added: for minimum subsets tree
from pasta.Kruskal_MST import build_groups_MST


class PastaTeam (object):
    '''A blob for holding the appropriate merger, alignment, and tree_estimator tools
    as well as the TempFS object that keeps track of the directories that have
    been created by this process.

    '''
    def __init__(self, config):
        """Uses a configuration object `cfg` to get reference for the tools the user
        has chosen.
        """
        try:
            max_mem_mb = config.sate.max_mem_mb
            self._temp_fs = TempFS()
            self.aligner = config.create_aligner(temp_fs=self._temp_fs)
            self.aligner.max_mem_mb = max_mem_mb
            self.hmmeralign = config.create_aligner(temp_fs=self._temp_fs, name = "hmmeralign")
            self.merger = config.create_merger(temp_fs=self._temp_fs)
            self.merger.max_mem_mb = max_mem_mb
            self.tree_estimator = config.create_tree_estimator(temp_fs=self._temp_fs)
            self.raxml_tree_estimator = config.create_tree_estimator(name='Raxml', temp_fs=self._temp_fs)
            self.subsets = {} # needed for pastamerger
            self.alignmentjobs = [] # needed for pastamerger
            self.treeshrink_wrapper = config.create_treeshrink_wrapper(temp_fs=self._temp_fs)
        except AttributeError:
            raise
            raise ValueError("config cannot be None unless all of the tools are passed in.")
    def get_temp_fs(self):
        return self._temp_fs
    temp_fs = property(get_temp_fs)

class AcceptMode:
    BLIND_MODE, NONBLIND_MODE = list(range(2))

class PastaJob (TreeHolder):
    """The top-level PASTA algorithm class.  The run_pasta method iteratively
    decomposes the tree and does tree searching until the termination criterion
    is reached.
    """
    BEHAVIOUR_DEFAULTS = {  'time_limit' : 24*60*60.0,
                            'iter_limit' : -1,
                            'time_without_imp_limit' : -1,
                            'iter_without_imp_limit' : -1,
                            'blind_after_total_time' : -1,
                            'blind_after_total_iter' : -1,
                            'blind_after_time_without_imp' : -1,
                            'blind_after_iter_without_imp' : -1,
                            'after_blind_time_term_limit' : 24*60*60.0,
                            'after_blind_iter_term_limit' : -1,
                            'after_blind_time_without_imp_limit': -1,
                            'after_blind_iter_without_imp_limit': 1,
                            'move_to_blind_on_worse_score' : False,
                            'blind_mode_is_final' : True,
                            'max_subproblem_size' : 200,
			    'max_subtree_diameter': 1.3,	
                            'max_subproblem_frac' : 0.2,
                            'start_tree_search_from_current' : False,
                            'keep_realignment_temporaries' : False,
                            'keep_iteration_temporaries' : False,
                            'return_final_tree_and_alignment' : False,
                            'mask_gappy_sites' : 1,
                            'build_MST' : False,
                            'treeshrink_filter': False
                        }
    def configuration(self):
        d = {}
        for k in PASTAAlignerJob.BEHAVIOUR_DEFAULTS.keys():
            d[k] = getattr(self, k)
        for k in PastaJob.BEHAVIOUR_DEFAULTS.keys():
            d[k] = getattr(self, k)
        return d

    def __init__(self, multilocus_dataset, pasta_team, tree=None, name=None, **kwargs):
        """
        Note that the `tree` will be resolved using DendroPy's resolve_polytomies
            function.
        """
        TreeHolder.__init__(self,
                            multilocus_dataset.dataset, 
                            force_fully_resolved=True)
        if tree is not None:
            resolve_polytomies(tree, update_splits=True)

        self.blind_mode_is_final = True
        self.is_stuck_in_blind = False
        behavior = copy.copy(PastaJob.BEHAVIOUR_DEFAULTS)
        behavior.update(PASTAAlignerJob.BEHAVIOUR_DEFAULTS)
        behavior.update(kwargs)
        self.__dict__.update(behavior)

        self._job_lock = Lock()
        self.multilocus_dataset = multilocus_dataset
        self.pasta_team = pasta_team
        self.tree = tree
        self.score = kwargs.get('score', None)
        self.best_score = None

        self._tree_build_job = None
        self._pasta_decomp_job = None
        self._reset_jobs()
        self.pastamerge = True

        self._status_message_func = kwargs.get('status_messages')

        # right now max_subproblem_size can be an integer or proportion.
        # These two run modes should be separate parameters!

        # if self.max_subproblem_size <= 1.0:
        #   self.max_subproblem_size *= alignment.get_num_taxa()
        # self.max_subproblem_size = max(1, int(round(self.max_subproblem_size)))

        if isinstance(self.break_strategy, str):
            self.break_strategy = [self.break_strategy]

        self._reset_current_run_settings()
        self.killed = False

    def _reset_current_run_settings(self):
        self.start_time = None
        self.current_iteration = 0
        self.last_improvement_time = None
        self.num_iter_since_imp = 0
        self.is_stuck_in_blind = False
        self.switch_to_blind_iter = None
        self.switch_to_blind_timestamp = None
        self._termination_trigger = None
        self._blindmode_trigger = None
        self._pasta_alignment_job = None
        self._tree_build_job = None
        self.curr_iter_align_tmp_filename = None
        self.curr_iter_tree_tmp_filename = None
        self.best_tree_tmp_filename = None
        self.best_alignment_tmp_filename = None


    def _reset_jobs(self):
        self.tree_build_job = None
        self.pasta_decomp_job = None

    # We lock the job objects because kill might be called from another thread
    def get_pasta_alignment_job(self):
        return self._pasta_alignment_job

    def set_pasta_alignment_job(self, val):
        self._job_lock.acquire()
        self._pasta_alignment_job = val
        self._job_lock.release()
    pasta_aligner_job = property(get_pasta_alignment_job, set_pasta_alignment_job)

    def get_tree_build_job(self):
        return self._tree_build_job

    def set_tree_build_job(self, val):
        self._job_lock.acquire()
        self._tree_build_job = val
        self._job_lock.release()
    tree_build_job = property(get_tree_build_job, set_tree_build_job)

    def _curr_running_times(self):
        "Duration of the current run (in seconds)"
        if self.start_time is None:
            return 0.0
        curr_time = time.time()
        return curr_time - self.start_time

    def _curr_time_since_imp(self):
        "Duration of the current run (in seconds)"
        if self.start_time is None:
            return 0.0
        curr_time = time.time()
        return curr_time - self.last_improvement_time

    def _keep_iterating(self):

        # PASTA has run for at least 'iter_limit' iterations
        if self.iter_limit >= 0 and self.current_iteration >= self.iter_limit:
            self._termination_trigger = 'iter_limit'
            return False

        # PASTA has run for 'iter_without_imp_limit' iterations and the likelihood score does not improve
        if self.iter_without_imp_limit >= 0 and self.num_iter_since_imp >= self.iter_without_imp_limit:
            self._termination_trigger = 'iter_without_imp_limit'
            return False

        # PASTA is in BlIND mode and it has run 'after_blind_iter_term_limit' iterations since then
        if (self.switch_to_blind_iter is not None) and self.after_blind_iter_term_limit >= 0:
            if (self.current_iteration - self.switch_to_blind_iter) >= self.after_blind_iter_term_limit:
                self._termination_trigger = 'after_blind_iter_term_limit'
                return False

        # PASTA has run 'time_limit' seconds
        if self.time_limit >= 0.0:
            running_time = self._curr_running_times()
            if running_time > self.time_limit:
                self._termination_trigger = 'time_limit'
                return False

        # PASTA has run 'time_without_imp_limit' seconds and the likelihood score does not improve
        if self.time_without_imp_limit >= 0.0:
            time_since_last_imp = self._curr_time_since_imp()
            if time_since_last_imp > self.time_without_imp_limit:
                self._termination_trigger = 'time_without_imp_limit'
                return False

        # PASTA is in BlIND mode and it has run 'after_blind_time_term_limit' seconds since then
        if (self.switch_to_blind_timestamp is not None) and self.after_blind_time_term_limit >= 0.0:
            time_since_switching_to_blind = time.time() - self.switch_to_blind_timestamp
            if time_since_switching_to_blind > self.after_blind_time_term_limit:
                self._termination_trigger = 'after_blind_time_term_limit'
                return False

        # PASTA is in BlIND mode and it has run 'self.after_blind_time_without_imp_limit' seconds since then
        # without improvements in likelihood score
        if (self.switch_to_blind_timestamp is not None) and self.after_blind_time_without_imp_limit >= 0.0:
            time_since_last_imp = self._curr_time_since_imp()
            time_since_switching_to_blind = time.time() - self.switch_to_blind_timestamp

            if time_since_switching_to_blind > self.after_blind_time_without_imp_limit:
                self._termination_trigger = 'after_blind_time_without_imp_limit'
                return False

        # PASTA is in BlIND mode and it has run 'self.after_blind_iter_without_imp_limit' iterations since then
        # without improvements in likelihood score
        if (self.switch_to_blind_timestamp is not None) and self.after_blind_iter_without_imp_limit >= 0:
            if self.num_iter_since_imp >= self.after_blind_iter_without_imp_limit:
                self._termination_trigger = 'after_blind_iter_without_imp_limit'
                return False

        return True

    def _get_break_strategy(self, break_strategy_index):
        try:
            return self.break_strategy[break_strategy_index]
        except:
            return None

    def _get_accept_mode(self, new_score, break_strategy_index=0):
        # If we have several break_strategies,  we do not go to "BLIND_MODE" until we have tried them all out
        if break_strategy_index < (len(self.break_strategy) - 1):
            return AcceptMode.NONBLIND_MODE
        if self.is_stuck_in_blind:
            return AcceptMode.BLIND_MODE
        #if self.move_to_blind_on_worse_score and ( (new_score <= self.score) or self.score is None ):
        if self.move_to_blind_on_worse_score and (new_score <= self.score):
            self._blindmode_trigger = 'move_to_blind_on_worse_score'
            return AcceptMode.BLIND_MODE
        if (self.blind_after_total_iter >= 0) and (self.current_iteration >= self.blind_after_total_iter):
            self._blindmode_trigger = 'blind_after_total_iter'
            return AcceptMode.BLIND_MODE
        if (self.blind_after_iter_without_imp >= 0) and (self.num_iter_since_imp >= self.blind_after_iter_without_imp):
            self._blindmode_trigger = 'blind_after_iter_without_imp'
            return AcceptMode.BLIND_MODE
        if self.blind_after_total_time >= 0.0:
            running_time = self._curr_running_times()
            if running_time > self.blind_after_total_time:
                self._blindmode_trigger = 'blind_after_total_time'
                return AcceptMode.BLIND_MODE
        if self.blind_after_time_without_imp >= 0.0:
            time_since_last_imp = self._curr_time_since_imp()
            if time_since_last_imp > self.blind_after_time_without_imp:
                self._blindmode_trigger = 'blind_after_time_without_imp'
                _LOG.debug("blind_after_time_without_imp reached")

        return AcceptMode.NONBLIND_MODE

    def store_optimum_results(self, new_multilocus_dataset, new_tree_str, new_score, curr_timestamp):
        self.best_multilocus_dataset = new_multilocus_dataset.new_with_shared_meta()
        for locus_alignment in new_multilocus_dataset:
            self.best_multilocus_dataset.append(copy.copy(locus_alignment))
        self.best_tree_str = new_tree_str
        self.best_score = new_score
        self.last_improvement_time = curr_timestamp
        self.num_iter_since_imp = 0
        self.best_tree_tmp_filename = self.curr_iter_tree_tmp_filename
        self.best_alignment_tmp_filename = self.curr_iter_align_tmp_filename


    #def build_subsets_tree(self, curr_tmp_dir_par):
    def build_subsets_tree(self, curr_tmp_dir_par,build_min_tree=True):
    # uym2 added: add option for MST
        if build_min_tree:
            _LOG.debug("START building Minimum Spanning Tree")
            grouping = {}
            groupName2jobName = {}
            
            for node in self.tree._tree.leaf_node_iter():
                groupName = self.pasta_team.subsets[node.taxon.label].tmp_dir_par[len(curr_tmp_dir_par)+1:]
                grouping[node.taxon.label] = groupName.replace("/","")
                groupName2jobName[groupName] = self.pasta_team.subsets[node.taxon.label]
            
            subsets_tree = build_groups_MST(self.tree._tree,grouping)
 
            for node in subsets_tree.postorder_node_iter():
               if node.is_leaf():
                   node.taxon.label = node.taxon.label.replace("d","/d")
               node.label = node.label.replace("d","/d") 

            self.pasta_team.subsets = groupName2jobName
            MST = PhylogeneticTree(subsets_tree) 
            _LOG.debug("Spanning tree is:\n %s" %MST)
            return MST
    ###################################


        _LOG.debug("START building heuristic spanning tree")

        translate={}
        t2 = {}
        for node in self.tree._tree.leaf_node_iter():
            nalsj = self.pasta_team.subsets[node.taxon.label]            
            newname = nalsj.tmp_dir_par[len(curr_tmp_dir_par)+1:]
            translate[node.taxon.label] = newname
            t2[newname] = set([nalsj])            
        subsets_tree = PhylogeneticTree(Tree.get(data=self.tree_str,schema='newick'))
        for node in subsets_tree._tree.leaf_node_iter():
            node.alignment_subset_job = t2[translate[node.taxon.label]]
            #node.alignment_subset_job = t2[node.taxon]
        del t2
        del translate
        _LOG.debug("leafs labeled")        
        #subsets_tree._tree.infer_taxa()
        #_LOG.debug("fake taxa inferred")                   
        #Then make sure the tree is rooted at a branch (not at a node). 
        if len(subsets_tree._tree.seed_node.child_nodes()) > 2:
            for c in subsets_tree._tree.seed_node.child_nodes():
                if c.edge.is_internal():
                    break
            subsets_tree._tree.is_rooted = True
            subsets_tree._tree.reroot_at_edge(c.edge,length1=c.edge.length/2., 
                                              length2=c.edge.length/2., suppress_unifurcations=False)                        
        _LOG.debug("Subset Labeling (start):\n%s" %str(subsets_tree.compose_newick(suppress_rooting=False))[0:5000])
        #_LOG.debug("Subset Labeling (start):\n%s" %str(len(subsets_tree._tree.seed_node.child_nodes())))
        # Then label internal branches based on their children, and collapse redundant edges. 
        for node in subsets_tree._tree.postorder_internal_node_iter():
            # my label is the intersection of my children, 
            # unless the intersection is empty, in which case it is the union
            if not hasattr(node, "alignment_subset_job") or node.alignment_subset_job is None:
                node.alignment_subset_job = set.intersection(*[c.alignment_subset_job for c in node.child_nodes()])
                if not node.alignment_subset_job:
                    node.alignment_subset_job = set.union(*[c.alignment_subset_job for c in node.child_nodes()])
            # Now go ahead and prune any child whose label encompasses my label. 
            # Use indexing instead of iteration, because with each collapse, 
            # new children can be added, and we want to process them as well.                         
            i = 0;
            while i < len(node.child_nodes()):                                
                c = node.child_nodes()[i]
                if node.alignment_subset_job.issubset(c.alignment_subset_job):
                    # Dendropy does not collapsing and edge that leads to a tip. Remove instead
                    if c.child_nodes():
                        c.edge.collapse()                                    
                    else:
                        node.remove_child(c)                      
                else:
                    i += 1
            
            node.label = "+".join(nj.tmp_dir_par[len(curr_tmp_dir_par)+1:] for nj in node.alignment_subset_job)
            if node.is_leaf():
                node.taxon = subsets_tree._tree.taxon_namespace.new_taxon(label=node.label)
            
        _LOG.debug("Before final round, the tree is:\n %s" %str(subsets_tree.compose_newick(suppress_rooting=False))[0:5000])
        # Now, the remaining edges have multiple labels. These need to
        # be further resolved. Do it by minimum length
        #   First find all candidate edges that we might want to contract
        candidate_edges = set()
        for e in subsets_tree._tree.postorder_edge_iter():
            if e.tail_node and e.head_node.alignment_subset_job.intersection(e.tail_node.alignment_subset_job):
                candidate_edges.add( (e.length,e) )
        #   Then sort the edges, and start removing them one by one
        #   only if an edge is still having intersecting labels at the two ends                                                    
        candidate_edges = sorted(candidate_edges, key=lambda x:  x[0] if x[0] else -1)       
        for (el, edge) in candidate_edges:
            I = edge.tail_node.alignment_subset_job.intersection(edge.head_node.alignment_subset_job)
            if I:
                edge.tail_node.alignment_subset_job = I 
                if edge.head_node.child_nodes():
                    #edge.collapse(adjust_collapsed_head_children_edge_lengths=True)
                    edge.collapse()
                else:
                    edge.tail_node.remove_child(edge.head_node)
        # Make sure the tree is correct, remove the actual jobs
        # from nodes (can cause deep-copy problems), assign a label to each
        # node, and keep a mapping between the labels and actual alignment job objects
        self.pasta_team.subsets = {} # Let this now map from subset labels to the actual alignment jobs
        for node in subsets_tree._tree.postorder_node_iter():
            assert len(node.alignment_subset_job) == 1
            nalsj = node.alignment_subset_job.pop()
            node.alignment_subset_job = None 
            node.label = nalsj.tmp_dir_par[len(curr_tmp_dir_par)+1:]#only find last part of the name
            self.pasta_team.subsets[node.label] = nalsj
            if node.is_leaf():
                # Add a dummy taxon, or else dendropy can get confused
                node.taxon = subsets_tree._tree.taxon_namespace.new_taxon(label=node.label)
        #subsets_tree._tree.infer_taxa()
        _LOG.debug("Spanning tree is:\n %s" %subsets_tree)
        labels = [nd.label for nd in subsets_tree._tree.postorder_node_iter()]
        if len(set(labels)) != len(labels):
            import collections
            raise Exception("Duplicate names found %s" %"\n".join
                   (item for item, count in 
                    collections.Counter(labels).items() if count > 1))
           
        return subsets_tree
        
    def run(self, tmp_dir_par, pasta_products=None):
        assert(os.path.exists(tmp_dir_par))

        self._reset_current_run_settings()
        self._reset_jobs()

        self.start_time = time.time()
        self.last_improvement_time = self.start_time

        num_non_update_iter = 0

        configuration = self.configuration()
        # Here we check if the max_subproblem_frac is more stringent than max_subproblem_size
        frac_max = int(math.ceil(self.max_subproblem_frac*self.tree.n_leaves))
        if frac_max > self.max_subproblem_size:
            configuration['max_subproblem_size'] = frac_max
        MESSENGER.send_info('Max subproblem set to {0}'.format(
                configuration['max_subproblem_size']))
        if configuration['max_subproblem_size'] >= self.tree.n_leaves:
            MESSENGER.send_warning('''\n
WARNING: you have specified a max subproblem ({0}) that is equal to or greater
    than the number of taxa ({0}). Thus, the PASTA algorithm will not be invoked
    under the current configuration (i.e., no tree decomposition will occur).
    If you did not intend for this behavior (which you probably did not since
    you are using PASTA) please adjust your settings for the max subproblem and
    try running PASTA again. If you intended to use PASTA to align your data with
    the specified aligner tool *without* any decomposition, you can ignore this
    message.\n'''.format(configuration['max_subproblem_size'],
                       self.tree.n_leaves))
        if configuration['max_subproblem_size'] == 1:
             MESSENGER.send_error(''' You have specified a max subproblem size of 1. PASTA requires a max subproblem size of at least 2.  ''')
             sys.exit(1)

        delete_iteration_temps = not self.keep_iteration_temporaries
        delete_realignment_temps = delete_iteration_temps or (not self.keep_realignment_temporaries)
        configuration['delete_temps'] = delete_realignment_temps

        while self._keep_iterating():
            record_timestamp(os.path.join(tmp_dir_par, 'start_pastaiter_timestamp.txt'))

            # create a subdirectory for this iteration
            curr_iter_tmp_dir_par = os.path.join(tmp_dir_par, 'step' + str(self.current_iteration))
            curr_iter_tmp_dir_par = self.pasta_team.temp_fs.create_subdir(curr_iter_tmp_dir_par)
            _LOG.debug('directory %s created' % curr_iter_tmp_dir_par)
            break_strategy_index = 0
            this_iter_score_improved = False

            while True:
                break_strategy =  self._get_break_strategy(break_strategy_index)
                if not bool(break_strategy):
                    break
                context_str = "iter%d-%s" % (self.current_iteration, break_strategy)
                # create a subdirectory for this iteration/break_strategy
                curr_tmp_dir_par = os.path.join(curr_iter_tmp_dir_par, break_strategy)
                curr_tmp_dir_par = self.pasta_team.temp_fs.create_subdir(curr_tmp_dir_par)
                record_timestamp(os.path.join(curr_tmp_dir_par, 'start_align_timestamp.txt'))
                # Align (with decomposition...)
                self.status('Step %d. Realigning with decomposition strategy set to %s' % (self.current_iteration, break_strategy))
                if self.killed:
                    raise RuntimeError("PASTA Job killed")
                tree_for_aligner = self.get_tree_copy()
                aligner = PASTAAlignerJob(multilocus_dataset=self.multilocus_dataset,
                                         pasta_team=self.pasta_team,
                                         tree=tree_for_aligner,
                                         tmp_base_dir=curr_tmp_dir_par,
                                         reset_recursion_index=True,
                                         skip_merge=self.pastamerge,
                                         **configuration)
                self.pasta_aligner_job = aligner
                aligner.launch_alignment(break_strategy=break_strategy,
                                         context_str=context_str)                
                if self.pastamerge:
                    _LOG.debug("Build PASTA merge jobs")
                    subsets_tree = self.build_subsets_tree(curr_tmp_dir_par,self.build_MST)
                    if len(self.pasta_team.subsets) == 1:
                        # can happen if there are no decompositions
                        for job in self.pasta_team.alignmentjobs:
                            jobq.put(job)
                        new_multilocus_dataset = list(self.pasta_team.subsets.values())[0].get_results()
                    else:
                        pariwise_tmp_dir_par = os.path.join(curr_tmp_dir_par, "pw")
                        pariwise_tmp_dir_par = self.pasta_team.temp_fs.create_subdir(pariwise_tmp_dir_par)    
                        pmj = PASTAMergerJob(multilocus_dataset=self.multilocus_dataset,
                                             pasta_team=self.pasta_team,
                                             tree=subsets_tree,
                                             tmp_base_dir=pariwise_tmp_dir_par,
                                             reset_recursion_index=True,   
                                             #delete_temps2=False,                                      
                                             **configuration)
                                                
                        pmj.launch_alignment(context_str=context_str)
                        
                        # Start alignment jobs
                        for job in self.pasta_team.alignmentjobs:
                            jobq.put(job)
                            
                            
                        new_multilocus_dataset = pmj.get_results()
                        del pmj  
                    
                    self.pasta_team.alignmentjobs = []
                    self.pasta_team.subsets = {}                                                                  
                else:          
                    new_multilocus_dataset = aligner.get_results()
                
                _LOG.debug("Alignment obtained. Preparing for tree.")
                self.pasta_aligner_job = None
                del aligner

                record_timestamp(os.path.join(curr_tmp_dir_par, 'start_treeinference_timestamp.txt'))
                
                # Tree inference
                if self.start_tree_search_from_current:
                    start_from = self.tree
                else:
                    start_from = None
                self.status('Step %d. Alignment obtained. Tree inference beginning...' % (self.current_iteration))
                if self.killed:
                    raise RuntimeError("PASTA Job killed")                             
           
            
                tbj = self.pasta_team.tree_estimator.create_job(new_multilocus_dataset,
                                                               starting_tree=start_from,
                                                               num_cpus=self.num_cpus,
                                                               context_str=context_str + " tree",
                                                               tmp_dir_par=curr_tmp_dir_par,
                                                               delete_temps=delete_iteration_temps,
                                                               pasta_products=pasta_products,
                                                               step_num=self.current_iteration,
                                                               mask_gappy_sites = self.mask_gappy_sites)
                prev_curr_align = self.curr_iter_align_tmp_filename
                prev_curr_tree = self.curr_iter_tree_tmp_filename
                self.curr_iter_align_tmp_filename = pasta_products.get_abs_path_for_iter_output(self.current_iteration, TEMP_SEQ_ALIGNMENT_TAG, allow_existing=True)
                self.curr_iter_tree_tmp_filename = pasta_products.get_abs_path_for_iter_output(self.current_iteration, TEMP_TREE_TAG, allow_existing=True)

                self.tree_build_job = tbj
                jobq.put(tbj)
                new_score, new_tree_str = tbj.get_results()
                self.tree_build_job = None
                del tbj
                if self.killed:
                    raise RuntimeError("PASTA Job killed")

                record_timestamp(os.path.join(curr_tmp_dir_par, 'end_treeinference_timestamp.txt'))
                curr_timestamp = time.time()

                                
                _LOG.debug("Tree obtained. Checking for acceptance.")

                this_iter_score_improved = ( self.best_score is None ) or ( new_score > self.best_score )

                accept_iteration =  ( this_iter_score_improved or 
                                      self._get_accept_mode(new_score=new_score, break_strategy_index=break_strategy_index) == AcceptMode.BLIND_MODE )

                if self.score is None:
                    self.score = new_score


                if accept_iteration:
                    if this_iter_score_improved:
                        self.status('realignment accepted and score improved.')
                    else:
                        self.status('realignment accepted despite the score not improving.')

                    self.status('current score: %s, best score: %s' % (self.score, self.best_score) )
                     
                    
                    _LOG.debug("Realignment and tree are accepted. Start treeshrink filtering.")
                   
                    if self.treeshrink_filter:
                        self.status("Running TreeShrink on the estimated tree")
                        aln_fn = self.curr_iter_align_tmp_filename
                        tree_fn = self.curr_iter_tree_tmp_filename
                        tsj = self.pasta_team.treeshrink_wrapper.create_job(aln_fn,
                                                                       new_multilocus_dataset[0].datatype,
                                                                       tree_fn,
                                                                       context_str=context_str + " treeshrink",
                                                                       tmp_dir_par=curr_tmp_dir_par,
                                                                       delete_temps=delete_iteration_temps,
                                                                       pasta_products=pasta_products,
                                                                       step_num=self.current_iteration)
                        jobq.put(tsj) 
                        shrunk_aln,shrunk_tree_str = tsj.get_results()
                        new_multilocus_dataset = shrunk_aln
                        new_tree_str = shrunk_tree_str
                        self.curr_iter_align_tmp_filename = pasta_products.get_abs_path_for_iter_output(self.current_iteration, TEMP_SHRUNK_ALIGNMENT_TAG, allow_existing=True)
                        self.curr_iter_tree_tmp_filename = pasta_products.get_abs_path_for_iter_output(self.current_iteration, TEMP_SHRUNK_TREE_TAG, allow_existing=True)
                    else:
                        self.status("TreeShrink option has been turned off!")
                    
                    self.score = new_score
                    self.multilocus_dataset = new_multilocus_dataset
                    self.tree_str = new_tree_str
                    
                    if this_iter_score_improved:
                        self.store_optimum_results(new_multilocus_dataset,
                                new_tree_str,
                                new_score,
                                curr_timestamp)

                    if self._get_accept_mode(new_score=new_score, break_strategy_index=break_strategy_index) == AcceptMode.BLIND_MODE:
                        if self.blind_mode_is_final:
                            self.is_stuck_in_blind = True
                            if self.switch_to_blind_timestamp is None:
                                if self._blindmode_trigger:
                                    _LOG.debug("Blind runmode trigger = %s" % self._blindmode_trigger)
                                self.switch_to_blind_iter = self.current_iteration
                                self.switch_to_blind_timestamp = curr_timestamp
        
                    # we do not want to continue to try different breaking strategies for this iteration so we break
                    break
                else:
                    self.status('realignment NOT accepted.')
                    self.curr_iter_align_tmp_filename = prev_curr_align
                    self.curr_iter_tree_tmp_filename = prev_curr_tree 

                break_strategy_index += 1

                # self.status('current score: %s, best score: %s' % (self.score, self.best_score) )
                
            if not this_iter_score_improved:
                self.num_iter_since_imp += 1
        
                
            self.current_iteration += 1

        if self._termination_trigger:
            _LOG.debug("Termination trigger = %s" % self._termination_trigger)
        record_timestamp(os.path.join(tmp_dir_par, 'end_pastaiter_timestamp.txt'))

        ### TODO: if configuration is 'return_final_iter_TreeAndAlignpair', then skip the following three lines
        if not self.return_final_tree_and_alignment:
            self.multilocus_dataset = self.best_multilocus_dataset.new_with_shared_meta()
            for locus_alignment in self.best_multilocus_dataset:
                self.multilocus_dataset.append(copy.copy(locus_alignment))
            self.tree_str = self.best_tree_str
            self.score = self.best_score
        else:
            assert self.multilocus_dataset is not None
            assert self.tree_str is not None
            assert self.score is not None

    def status(self, message):
        if self._status_message_func:
            self._status_message_func(message)

    def kill(self):
        self.killed = True
        j = self.pasta_aligner_job
        if j:
            self.pasta_aligner_job = None
            j.kill()
        j = self.tree_build_job
        if j:
            self.tree_build_job = None
            j.kill()
