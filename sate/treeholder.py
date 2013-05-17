#!/usr/bin/env python

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

import dendropy

from cStringIO import StringIO
from sate.tree import PhylogeneticTree
from sate.errors import TaxaLabelsMismatchError
from sate import get_logger, MESSENGER
_LOG = get_logger(__name__)

# Provide a random number generator
import random
POLYTOMY_RNG = random.Random()

from dendropy.treesplit import delete_outdegree_one

def resolve_polytomies(tree, update_splits=False, rng=None):
    """
    Copied from more recent DendroPy than the version that we bundle...
    
    Arbitrarily resolve polytomies using 0-length splits.

    If `rng` is an object with a sample() method then the polytomy will be
        resolved by sequentially adding (generating all tree topologies
        equiprobably
        rng.sample() should behave like random.sample()
    If `rng` is not passed in, then polytomy is broken deterministically by
        repeatedly joining pairs of children.
    """
    from dendropy import Node
    polytomies = []
    if rng is None:
        rng = POLYTOMY_RNG
    for node in tree.postorder_node_iter():
        if len(node.child_nodes()) > 2:
            polytomies.append(node)
    for node in polytomies:
        children = node.child_nodes()
        nc = len(children)
        if nc > 2:
            #if nc == 3 and node.parent_node is None:
            #    continue
            to_attach = children[2:]
            for child in to_attach:
                node.remove_child(child)
            attachment_points = children[:2] + [node]
            while len(to_attach) > 0:
                next_child = to_attach.pop()
                next_sib = rng.sample(attachment_points, 1)[0]
                next_attachment = Node()
                next_attachment.edge.length = 0.0
                p = next_sib.parent_node
                if p is None:
                    c_list = list(next_sib.child_nodes())
                    next_sib.add_child(next_attachment)
                    next_sib.add_child(next_child)
                    for child in c_list:
                        next_sib.remove_child(child)
                        next_attachment.add_child(child)
                else:
                
                    p.add_child(next_attachment)
                    p.remove_child(next_sib)
                    next_attachment.add_child(next_sib)
                    next_attachment.add_child(next_child)
                attachment_points.append(next_attachment)
    if update_splits:
        tree.update_splits()

def check_taxon_labels(taxon_set, dataset):
    ts = set([i for tset in dataset.taxon_sets for i in tset.labels()])
    ds = set(taxon_set.labels())
    extra = ds - ts
    missing = ts - ds
    return extra, missing

def read_trees_into_dataset(dataset, tree_stream, starting_tree=False):
    if starting_tree:        
        try:
            dataset.read_from_stream(tree_stream,
                schema='NEWICK', taxon_set=dataset.taxon_sets[0])
        except KeyError as e:
            m = str(e)
            m = m[1:m.find("TaxonSet")] + "sequences but present in tree"            
            raise TaxaLabelsMismatchError(                
                 'There are taxon label mismatches between the starting tree '
                 'and sequences...\n'
                 '%s\n' %m)
        st = dataset.tree_lists[-1][0]        
        if len(st.leaf_nodes()) != len(dataset.taxon_sets[0]):
            missing = [t.label for t in set(dataset.taxon_sets[0]) - set((n.taxon for n in st.leaf_nodes()))]
            raise TaxaLabelsMismatchError(                
                 'There are taxon label mismatches between the starting tree '
                 'and sequences...\n'
                 'In sequences, not tree: {0}\n'.format(','.join(missing)) )
        _LOG.debug("reading tree finished")
    elif dataset.taxon_sets:
        dataset.read_from_stream(tree_stream, schema='NEWICK', taxon_set=dataset.taxon_sets[0])
    else:
        dataset.read_from_stream(tree_stream, schema='NEWICK')
    return  dataset.tree_lists[-1]

def read_and_encode_splits(dataset, tree_stream, starting_tree=False):
    """Reads the file-like object `tree_stream` as a source of trees for the
    the taxa found in dataset. and then encodes the splits of the nodes of the trees.
    This is a convenience function that bridges between dendropy 2 and 3 API's
    """
    _LOG.debug("NOT covered in tests")
    tree_list = read_trees_into_dataset(dataset, tree_stream,
            starting_tree=starting_tree)
    assert len(tree_list) == 1
    delete_outdegree_one(tree_list[0])
    return tree_list

def generate_tree_with_splits_from_str(tree_str, dataset, force_fully_resolved=False):
    '''Uses `tree_str` and `dataset` to create a PhylogeneticTree object
    and calls `calc_splits` on the object before returning it.
    '''
    _LOG.debug("start generating tree from string")
    tree_stream = StringIO(tree_str)
    tree_list = read_and_encode_splits(dataset, tree_stream)
    t = tree_list[0]
    if force_fully_resolved:
        _LOG.debug("start resolving polytomies")
        resolve_polytomies(t, update_splits=True)
        _LOG.debug("end of resolving polytomies")
    t = PhylogeneticTree(t)
    _LOG.debug("calculating splits")
    t.calc_splits()
    _LOG.debug("end generating tree from string")
    return t

class TreeHolder(object):
    '''Uses the tree attribute to provide a `tree_str` property, but also
        enables setting of the `tree_str` to update the tree.
    '''
    def __init__(self, dataset, force_fully_resolved=False):
        self.dataset = dataset
        self.tree = None
        self._force_fully_resolved=force_fully_resolved

    def get_tree_str(self):
        return self.tree.compose_newick() if self.tree else None

    def set_tree_str(self, tree_str):
        self.tree = generate_tree_with_splits_from_str(tree_str, 
                                                       self.dataset,
                                                       self._force_fully_resolved)

    tree_str = property(get_tree_str, set_tree_str)

    def get_tree_copy(self):
        '''Returns a deep copy of the tree instance.'''
        return generate_tree_with_splits_from_str(self.tree_str, self.dataset)

