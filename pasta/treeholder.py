#!/usr/bin/env python

# This file is part of PASTA and is forked from SATe

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

import dendropy

try:
    from cStringIO import StringIO
except:
    from io import StringIO
from pasta.tree import PhylogeneticTree
from pasta.errors import TaxaLabelsMismatchError
from pasta import get_logger
from dendropy import Tree, TreeList
_LOG = get_logger(__name__)
from dendropy.datamodel.treemodel import _convert_node_to_root_polytomy as convert_node_to_root_polytomy

# Provide a random number generator
import random

POLYTOMY_RNG = random.Random()

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
    _LOG.debug("start resolving polytomies")
    from dendropy import Node
    polytomies = []
    if rng is None:
        rng = POLYTOMY_RNG
    for node in tree.postorder_node_iter():
        if len(node.child_nodes()) > 2:
            polytomies.append(node)
            
    _LOG.debug("Found %d polytomies" %len(polytomies))
    for node in polytomies:
        children = node.child_nodes()
        nc = len(children)
        if nc > 2:
            if nc == 3 and node.parent_node is None:
                continue
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
    _LOG.debug("polytomies resolution - updating splits")
    if update_splits:
        tree.update_splits()
    _LOG.debug("polytomies resolved.")
    
def check_taxon_labels(taxon_namespace, dataset):
    ts = set([i for tset in dataset.taxon_namespaces for i in tset.labels()])
    ds = set(taxon_namespace.labels())
    extra = ds - ts
    missing = ts - ds
    return extra, missing


def read_newick_with_translate(stream,taxon_namespace):
    """
    Instantiates and returns a `DataSet` object based on the
    NEWICK-formatted contents read from the file-like object source
    `stream`.
    """
    ts = [t for t in tree_source_iter(stream=stream,
            translate_dict=translate_dict,
            allow_repeated_use=True)]
    return ts[0]



def read_trees_into_dataset(dataset, tree_stream, starting_tree=False, preserve_underscores=True):
    if starting_tree:        
        try:
            dataset.read_from_stream(tree_stream,
                schema='NEWICK', taxon_namespace=dataset.taxon_namespaces[0],preserve_underscores=preserve_underscores)
        except KeyError as e:
            m = str(e)
            m = m[1:m.find("TaxonSet")] + "sequences but present in tree"            
            raise TaxaLabelsMismatchError(                
                 'There are taxon label mismatches between the starting tree '
                 'and sequences...\n'
                 '%s\n' %m)
        st = dataset.tree_lists[-1][0]        
        if len(st.leaf_nodes()) != len(dataset.taxon_namespaces[0]):
            missing = [t.label for t in set(dataset.taxon_namespaces[0]) - set((n.taxon for n in st.leaf_nodes()))]
            raise TaxaLabelsMismatchError(                
                 'There are taxon label mismatches between the starting tree '
                 'and sequences...\n'
                 'In sequences, not tree: {0}\n'.format(','.join(missing)) )
        _LOG.debug("reading tree finished")
    elif dataset.taxon_namespaces:
        dataset.read_from_stream(tree_stream, schema='NEWICK', taxon_namespace=dataset.taxon_namespaces[0])
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
    #from dendropy.legacy.treesplit import delete_outdegree_one
    #delete_outdegree_one(tree_list[0])
    tree_list[0].suppress_unifurcations()
    convert_node_to_root_polytomy(tree_list[0].seed_node)
    return tree_list

def generate_tree_with_splits_from_str(tree_str, dataset, force_fully_resolved=False):
    '''Uses `tree_str` and `dataset` to create a PhylogeneticTree object
    and calls `calc_splits` on the object before returning it.
    '''
    _LOG.debug("start generating tree from string %s" %tree_str[0:200])
    tree_stream = StringIO(tree_str)
    tree_list = read_and_encode_splits(dataset, tree_stream)
    t = tree_list[0]
    _LOG.debug("tree  generated from string %s" %str(t)[0:200])
    #_LOG.debug("tree rooting %s" %str(t.is_rooted))
    return generate_tree_with_splits_from_tree(t, force_fully_resolved)
    
def generate_tree_with_splits_from_tree(t, force_fully_resolved=False):    
    if force_fully_resolved:        
        resolve_polytomies(t, update_splits=False)
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

