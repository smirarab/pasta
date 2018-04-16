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


import os
import copy
from threading import Lock
from pasta import get_logger
from pasta.tree import PhylogeneticTree
from dendropy.datamodel.treemodel import Tree
from pasta.alignment import  CompactAlignment
_LOG = get_logger(__name__)

from pasta.treeholder import TreeHolder
from pasta.scheduler import jobq, new_merge_event
from pasta.scheduler import TickableJob

from pasta.new_decomposition import midpoint_bisect, min_cluster_size_bisect

# uym2 modified: add min_size option
def bisect_tree(tree, breaking_edge_style='mincluster',min_size=0,max_size=None,max_diam=None):
    """Partition 'tree' into two parts
    """
    _LOG.debug("bisecting tree...")
    # uym2: midpoint decomposition (not in used for now)
    if breaking_edge_style == 'midpoint':
        _LOG.debug("breaking by midpoint")
        t1,t2 = midpoint_bisect(tree._tree,min_size=min_size)	
        if t1 is None or t2 is None:
            return None,None
        tree1 = PhylogeneticTree(t1)
        tree2 = PhylogeneticTree(t2)   
        return tree1,tree2
    ###############

    # uym2: min_cluster decomposition
    if breaking_edge_style == 'mincluster':
        _LOG.debug("breaking using min-cluster strategy")
        t1,t2 = min_cluster_size_bisect(tree._tree,max_size)
        tree1 = PhylogeneticTree(t1) if t1 else None
        tree2 = PhylogeneticTree(t2) if t2 else None
        return tree1,tree2
    ###############

    _LOG.debug("breaking by centroid")
    e = tree.get_breaking_edge(breaking_edge_style)
    _LOG.debug("breaking_edge length = %s, %s" % (e.length, breaking_edge_style) )
    snl = tree.n_leaves
    tree1, tree2 = tree.bipartition_by_edge(e)
    _LOG.debug("Tree 1 has %s nodes, tree 2 has %s nodes" % (tree1.n_leaves, tree2.n_leaves) )
    assert snl == tree1.n_leaves + tree2.n_leaves
    return tree1, tree2


class PASTAAlignerJob(TreeHolder, TickableJob):
    """A class that performs one alignment in the PASTA algorithm.
    
    If the size of the tree is <= than the threshold problem size then this
        will  be accomplished by calling the aligner (e.g. MAFFT).

    If the size of the tree is larger than the threshold problem size then the
        alignment is accomplished by decomposing:
            - the tree into two parts,
            - Creating a PASTAAlignerJob to align each part,
            - Merging the resulting alignments with the merger tool (e.g. Opal).
            
    """
    BEHAVIOUR_DEFAULTS = {  'break_strategy' : tuple(['mincluster']) ,
                            'max_subproblem_size' : 50,
            			    'max_subtree_diameter': 2.5,
                            'min_subproblem_size': 0,
                            'delete_temps' : True}
    RECURSION_INDEX = 0
    def __init__(self, 
                multilocus_dataset, 
                pasta_team, 
                tree,
                tmp_base_dir,
                tmp_dir_par=None,
                reset_recursion_index=False,
                skip_merge = False,
                **kwargs):
        self._job_lock = Lock()
        self._merge_queued_event = new_merge_event()
        TickableJob.__init__(self)
        TreeHolder.__init__(self, multilocus_dataset.dataset)
        behavior = copy.copy(PASTAAlignerJob.BEHAVIOUR_DEFAULTS)
        behavior.update(kwargs)
        for k in list(PASTAAlignerJob.BEHAVIOUR_DEFAULTS.keys()):
            setattr(self, k, behavior[k])
        self.multilocus_dataset = multilocus_dataset
        self.pasta_team = pasta_team
        self.tree = tree
        _LOG.debug("Number of children of the root at the start: %d" % len(tree._tree.seed_node.child_nodes()))
        self._subjob1 = None
        self._subjob2 = None
        self._align_job_list = None
        self._merge_job_list = None
        self.tmp_base_dir = tmp_base_dir
        self.context_str = ''
        self.killed = False
        self._dirs_to_cleanup = []
        self.tmp_dir_par = tmp_dir_par
        self.skip_merge = skip_merge
        if self.tmp_dir_par == None:
            self.tmp_dir_par = self.tmp_base_dir
        if reset_recursion_index:
            self.__class__.RECURSION_INDEX = 0
        self.finished = False

    def configuration(self):
        d = {}
        for k in list(PASTAAlignerJob.BEHAVIOUR_DEFAULTS.keys()):
            d[k] = getattr(self, k)
        return d

    def _reset_jobs(self):
        self.subjob1 = None
        self.subjob2 = None
        self.align_job_list = None
        self.merge_job_list = None

    def postprocess(self):
        self.tick_praents()

    def on_dependency_ready(self):
        ''' This is called when the "child" jobs are finished. 
        If self is a "alignment" subproblem, we just need to tick the parent job(s).
        If self is a "merge" subproblem, the first time all dependencies are met
        is when the children jobs finish. At this point, the actual merge job
        is started and queued, and that new merge job is added as a child of self.
        Now, when this new "merge" job finishes, self is completed, and we need
        to tick the parent(s).'''
        #_LOG.debug("Dependency ready!")
        if self.killed:
            raise RuntimeError("PastaAligner Job killed")
        if not self.align_job_list:
            if not self.merge_job_list:
                assert self.subjob1 and self.subjob2
                if not self.skip_merge:
                    self._start_merger()
                    return
        self.postprocess()
        

    def _start_merger(self):
        '''Blocks until the two "subjobs" are done 
        (with new implementation, they will be done, before this is even called)
            creates the merger job and puts it in the jobs queue, 
            cleans up the alignment subdirectories,
            signals an event that signifies the fact that the merge job is on queue,
            and then returns.
        
        Called by wait()
        '''
        if self.killed:
            raise RuntimeError("PastaAligner Job killed")
        assert(self.subjob1 is not None)
        result1 = self.subjob1.get_results()
        if self.killed:
            raise RuntimeError("PastaAligner Job killed")
        self.subjob1 = None
        assert(self.subjob2 is not None)
        result2 = self.subjob2.get_results()
        self.subjob2 = None
        if self.killed:
            raise RuntimeError("PastaAligner Job killed")
        assert(result1.get_num_loci() == result2.get_num_loci())
        
        mj_list = []
        for n, r1 in enumerate(result1):
            r2 = result2[n]
            cs = self.context_str + " merger" + str(n)
            mj = self.pasta_team.merger.create_job(r1,
                                                  r2,
                                                  tmp_dir_par=self.tmp_dir_par,
                                                  delete_temps=self.delete_temps,
                                                  context_str=cs)
            mj.add_parent_tickable_job(self)
            self.add_child(mj)
                        
            if self.killed:
                raise RuntimeError("PastaAligner Job killed")
            mj_list.append(mj)

        self.merge_job_list = mj_list
        for mj in mj_list:
            jobq.put(mj)

        if self.delete_temps:
            for d in self._dirs_to_cleanup:
                self.pasta_team.temp_fs.remove_dir(d)

        self._merge_queued_event.set()
    
    def _get_subjob_dir(self, num):
        '''Creates a numbered directory d1, d2, etc within tmp_dir_par.
        
        Called in bipartition_by_tree, and the directories are cleaned up
        at the end of _start_merger.
        '''
        assert(self.tmp_base_dir)
        rn = "r%d" % PASTAAlignerJob.RECURSION_INDEX
        dn = "d%d" % num
        r_dir = os.path.join(self.tmp_base_dir, rn)
        sd = os.path.join(r_dir, dn)
        if not os.path.exists(r_dir):
            self.pasta_team.temp_fs.create_subdir(r_dir)
        full_path_to_new_dir = self.pasta_team.temp_fs.create_subdir(sd)
        self._dirs_to_cleanup.append(full_path_to_new_dir)
        return full_path_to_new_dir


    def launch_alignment(self, tree=None, break_strategy=None, context_str=None):
        '''Puts a alignment job(s) in the queue and then return None
        
        get_results() must be called to get the alignment. Note that this call 
        may not be trivial in terms of time (the tree will be decomposed, lots
        of temporary files may be written...), but the call does not block until
        completion of the alignments.
        Rather it queues the alignment jobs so that multiple processors can be 
        exploited if they are available.
        '''
        if self.killed:
            raise RuntimeError("PastaAligner Job killed")

        if break_strategy is not None:
            self.break_strategy = break_strategy
        break_strategy = self.break_strategy
        if tree is not None:
            self.tree = tree
        self.expected_number_of_taxa = self.multilocus_dataset.get_num_taxa() # for debugging purposes
        self._reset_jobs()
        prefix = "self.multilocus_dataset.get_num_taxa = %d" % self.expected_number_of_taxa
        self.context_str = context_str
        if self.context_str is None:
            self.context_str = ''
        _LOG.debug("Comparing expected_number_of_taxa=%d and max_subproblem_size=%d\n" % (self.expected_number_of_taxa,  self.max_subproblem_size))

        if self.expected_number_of_taxa <= self.max_subproblem_size:
            _LOG.debug("%s...Calling Aligner" % prefix)
            aj_list = []
            for index, single_locus_sd in enumerate(self.multilocus_dataset):
                aj = self.pasta_team.aligner.create_job(single_locus_sd,
                                                       tmp_dir_par=self.tmp_dir_par,
                                                       delete_temps=self.delete_temps,
                                                       context_str=self.context_str + " align" + str(index))
                aj.add_parent_tickable_job(self)
                self.add_child(aj)
                
                aj_list.append(aj)
                if self.killed:
                    raise RuntimeError("PastaAligner Job killed")
                
                self.pasta_team.alignmentjobs.append(aj)
            
            self.align_job_list = aj_list
            
            if self.skip_merge:
                for taxa in self.tree.leaf_node_names():
                    self.pasta_team.subsets[taxa]=self
            else:
                for aj in aj_list:
                    jobq.put(aj)
        else:
            # added by uym2 on August 1st 2017
            subjob1, subjob2 = self.bipartition_by_tree(break_strategy)
            if subjob1 is None or subjob2 is None:
                return
            _LOG.debug("%s...Recursing" % prefix)
            # create the subjobs
            # the next line was modified by uym2 (August 1st 2017)
            self.subjob1 = subjob1
            self.subjob2 = subjob2
            # store this dir so we can use it in the merger
            if self.killed:
                raise RuntimeError("PastaAligner Job killed")

            self.subjob1.add_parent(self)
            self.subjob2.add_parent(self)
            self.add_child(self.subjob1)
            self.add_child(self.subjob2)

            self.subjob1.launch_alignment(break_strategy=break_strategy)
            if self.killed:
                raise RuntimeError("PastaAligner Job killed")
            self.subjob2.launch_alignment(break_strategy=break_strategy)
            if self.killed:
                raise RuntimeError("PastaAligner Job killed")

    # We lock the job objects because kill might be called from another thread
    def get_subjob1(self):
        return self._subjob1
    def set_subjob1(self, val):
        self._job_lock.acquire()
        self._subjob1 = val
        self._job_lock.release()
    subjob1 = property(get_subjob1, set_subjob1)

    def get_subjob2(self):
        return self._subjob2
    def set_subjob2(self, val):
        self._job_lock.acquire()
        self._subjob2 = val
        self._job_lock.release()
    subjob2 = property(get_subjob2, set_subjob2)

    def get_merge_job_list(self):
        return self._merge_job_list
    def set_merge_job_list(self, val):
        self._job_lock.acquire()
        self._merge_job_list = val
        self._job_lock.release()
    merge_job_list = property(get_merge_job_list, set_merge_job_list)

    def get_align_job_list(self):
        return self._align_job_list
    def set_align_job_list(self, val):
        self._job_lock.acquire()
        self._align_job_list = val
        self._job_lock.release()
    align_job_list = property(get_align_job_list, set_align_job_list)

    def get_allow_launch(self):
        return self._allow_launch
    def set_allow_launch(self, val):
        self._job_lock.acquire()
        self._allow_launch = val
        self._job_lock.release()
    allow_launch = property(get_allow_launch, set_allow_launch)

    def kill(self):
        self.killed = True
        _LOG.debug("killing pasta aligner jobs")
        j = self.subjob2
        if j:
            _LOG.debug("Killing subjob2")
            j.kill()
        j = self.subjob1
        if j:
            _LOG.debug("Killing subjob1")
            j.kill()            
        j_list = self.merge_job_list
        if j_list:
            for j in j_list:
                if j:
                    _LOG.debug("Killing merge job")
                    j.kill()
        j_list = self.align_job_list
        if j_list:
            for j in j_list:
                if j:
                    _LOG.debug("Killing align job")
                    j.kill()
        

    def wait(self):
        if self.killed:
            raise RuntimeError("PastaAligner Job killed")
        if self.finished:
            return
        try:
            j_list = self.align_job_list
            if j_list:
                for j in j_list:
                    j.wait()
            else:
                if not self.merge_job_list:
                    _LOG.debug("waiting for merge events to be queued")
                    self._merge_queued_event.wait()
                    _LOG.debug("merge events to be queued. signal received.")
                if self.merge_job_list:
                    for j in self.merge_job_list:
                        return j.wait()
        except KeyboardInterrupt:
            self.kill()

    def bipartition_by_tree(self, option):
        _LOG.debug("tree before bipartition by %s = %s ..." % (option, self.tree.compose_newick()[0:200]))

# uym2 modified: add min_size option
        tree1, tree2 = bisect_tree(self.tree, breaking_edge_style=option,min_size=self.min_subproblem_size,max_size=self.max_subproblem_size)
        if tree1 is None or tree2 is None:
            return [None,None]
        assert tree1.n_leaves > 0
        assert tree2.n_leaves > 0
        assert tree1.n_leaves + tree2.n_leaves == self.tree.n_leaves

        _LOG.debug("tree1 = %s ..." % tree1.compose_newick() [0:200])
        _LOG.debug("tree2 = %s ..." % tree2.compose_newick() [0:200])

        multilocus_dataset1 = self.multilocus_dataset.sub_alignment(tree1.leaf_node_names())
        multilocus_dataset2 = self.multilocus_dataset.sub_alignment(tree2.leaf_node_names())
        sd1 = self._get_subjob_dir(1)
        sd2 = self._get_subjob_dir(2)
        PASTAAlignerJob.RECURSION_INDEX += 1
        configuration = self.configuration()
        return [PASTAAlignerJob(multilocus_dataset=multilocus_dataset1,
                                pasta_team=self.pasta_team,
                                tree=tree1,
                                tmp_base_dir=self.tmp_base_dir,
                                tmp_dir_par=sd1,
                                skip_merge=self.skip_merge,
                                **configuration),
                PASTAAlignerJob(multilocus_dataset=multilocus_dataset2,
                                pasta_team=self.pasta_team,
                                tree=tree2,
                                tmp_base_dir=self.tmp_base_dir,
                                tmp_dir_par=sd2,
                                skip_merge=self.skip_merge,
                                **configuration)]

    def get_results(self):
        self.wait()
        if self.killed:
            raise RuntimeError("PastaAligner Job killed")
        j_list = self.align_job_list
        if j_list:
            r = self.multilocus_dataset.new_with_shared_meta()
            for j in j_list:
                r.append(j.get_results())
            #self.align_job_list = None
            self.finished = True
        else:
            j_list = self.merge_job_list
            if j_list:
                r = self.multilocus_dataset.new_with_shared_meta()
                for j in j_list:
                    r.append(j.get_results())
                #self.merge_job_list = None
                self.finished = True
            else:
                return None # this can happen if jobs are killed
        return r

class PASTAMergerJob(PASTAAlignerJob):
    
    def __init__(self, 
                multilocus_dataset,
                pasta_team, 
                tree,
                tmp_base_dir,
                tmp_dir_par=None,
                reset_recursion_index=False,
                delete_temps2 = None,
                **kwargs):
        PASTAAlignerJob.__init__(self, 
                                 multilocus_dataset,
                                 pasta_team, 
                                 tree, 
                                 tmp_base_dir, 
                                 tmp_dir_par, 
                                 reset_recursion_index,
                                 **kwargs                                 
                                 )
        if delete_temps2 is not None:
            self.delete_temps = delete_temps2
        self.skip_merge = None
        self.result_alignment = None
        self.results_lock = Lock()

    def launch_alignment(self, context_str=None):
        '''
        '''
        if self.killed:
            raise RuntimeError("PastaAligner Job killed")

        self._reset_jobs()
        self.context_str = context_str
        if self.context_str is None:
            self.context_str = ''
        node_count = self.tree.count_nodes()
        _LOG.debug("Recursive merge on a branch with %d subsets" % (node_count))
        prefix = "subsets tree: %s" %self.tree.compose_newick()[0:200]
        if node_count == 2:
            nodes = self.tree._tree.nodes()
            _LOG.debug("%s ... pairwise merge " % prefix)
            self.skip_merge = False
            self.subjob1 = self.pasta_team.subsets[nodes[0].label]           
            self.subjob2 = self.pasta_team.subsets[nodes[1].label]
            
            self.subjob1.add_parent(self)
            self.add_child(self.subjob1)

            self.subjob2.add_parent(self)
            self.add_child(self.subjob2)                                        
        else:
            _LOG.debug("%s ... recursing further " % prefix)
            self.skip_merge = True
            
            # Reroot near centroid edge
            ce = self.tree.get_centroid_edge(spanning=True)
            nr = ce.head_node if not ce.head_node.is_leaf() else ce.tail_node
            self.tree._tree.reroot_at_node(nr,suppress_unifurcations=False)            
            _LOG.debug("rerooted to: %s ..." % self.tree.compose_newick()[0:200])   
            # For each path from root to its children, create a new merge job         
            merge_job_list = []
            nr = self.tree._tree.seed_node
            children = nr.child_nodes()
            for keepchild in children:                
                remchilds = []                
                for remchild in children:
                    if remchild != keepchild:
                        remchilds.append(nr.reversible_remove_child(remchild, suppress_unifurcations=False))
                t1 = PhylogeneticTree(Tree(self.tree._tree))
                remchilds.reverse()
                for child in remchilds:
                    nr.reinsert_nodes(child)
                _LOG.debug("child = %s ..." % t1.compose_newick()[0:200])
                multilocus_dataset1 = self.multilocus_dataset.new_with_shared_meta()
                
                if t1.count_nodes() == 2:            
                    ns = t1._tree.nodes()
                    tmp_dir_par = self.get_pairwise_temp_dir(ns[0].label, ns[1].label)
                else:
                    tmp_dir_par = self.tmp_base_dir                    
                configuration = self.configuration()
                cj = PASTAMergerJob(multilocus_dataset=multilocus_dataset1,
                                    pasta_team=self.pasta_team,
                                    tree=t1,
                                    tmp_base_dir=self.tmp_base_dir,
                                    tmp_dir_par= tmp_dir_par,
                                    delete_temps2=False,
                                    **configuration)
                cj.add_parent(self)
                self.add_child(cj)                                
                merge_job_list.append(cj);
                        
            self.merge_job_list = merge_job_list
            
            # now launch these new merge jobs
            for merge_job in self.merge_job_list:
                if self.killed:
                    raise RuntimeError("PastaAligner Job killed")
                merge_job.launch_alignment()

            self._merge_queued_event.set()
            
            if self.killed:
                raise RuntimeError("PastaAligner Job killed")
        return

    def get_results(self):
        
        self.wait()
        
        self.results_lock.acquire()        
        if self.result_alignment is not None:
            self.results_lock.release()
            return self.result_alignment
        if self.killed:
            raise RuntimeError("PastaAligner Job killed")
        if self.align_job_list:
            raise RuntimeError("hmm, this should be empty")
        else:
            j_list = self.merge_job_list
            if j_list:
                # These are merge jobs that need transitivity merging             
                if self.skip_merge:
                    r = self.multilocus_dataset.new_with_shared_meta()
                    r.append(j_list[0].get_results()[0]) #TODO: this should be changed to be multi-locus
                    for j in j_list[1:]:
                        r[0].merge_in(j.get_results()[0]) #TODO: this should be changed to be multi-locus
                        j.clear_results_object()                        
                    #assert all(x.is_aligned() for x in r)
                else: # These are pairwise merges
                    r = self.multilocus_dataset.new_with_shared_meta()
                    for j in j_list:
                        a = CompactAlignment()
                        a.update_from_alignment(j.get_results())
                        r.append(a)
                self.result_alignment = r
                self.finished = True
            else:
                r = None # this can happen if jobs are killed
        self.results_lock.release()
        return r
        
    def clear_results_object(self):
        self.results_lock.acquire()
        for alg in self.result_alignment: 
            alg.clear()
        self.result_alignment = None
        self.results_lock.release()                                        
        
    def get_pairwise_temp_dir(self, label1, label2):
        ''' Get a temp file for a pairwise merge
        '''
        assert(self.tmp_base_dir)
        label = "%s_%s" %(label1,label2)
        sd = os.path.join(self.tmp_base_dir, label.replace("/","").replace("\\", ""))
        full_path_to_new_dir = self.pasta_team.temp_fs.create_subdir(sd)
        self._dirs_to_cleanup.append(full_path_to_new_dir)
        return full_path_to_new_dir

    def postprocess(self):
        self.get_results()
        self.tick_praents()
