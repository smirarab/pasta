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


from cStringIO import StringIO
from sate.tree import PhylogeneticTree
from sate import get_logger
_LOG = get_logger(__name__)

from dendropy.treesplit import delete_outdegree_one

def read_trees_into_dataset(dataset, tree_stream):
    if dataset.taxon_sets:
        dataset.read_from_stream(tree_stream, schema='NEWICK', taxon_set=dataset.taxon_sets[0])
    else:
        dataset.read_from_stream(tree_stream, schema='NEWICK')
    return  dataset.tree_lists[-1]

def read_and_encode_splits(dataset, tree_stream):
    """Reads the file-like object `tree_stream` as a source of trees for the
    the taxa found in dataset. and then encodes the splits of the nodes of the trees.
    This is a convenience function that bridges between dendropy 2 and 3 API's
    """
    _LOG.debug("NOT covered in tests")
    tree_list = read_trees_into_dataset(dataset, tree_stream)
    assert len(tree_list) == 1
    delete_outdegree_one(tree_list[0])
    return tree_list

def generate_tree_with_splits_from_str(tree_str, dataset):
    '''Uses `tree_str` and `dataset` to create a PhylogeneticTree object
    and calls `calc_splits` on the object before returning it.
    '''
    tree_stream = StringIO(tree_str)
    tree_list = read_and_encode_splits(dataset, tree_stream)
    t = PhylogeneticTree(tree_list[0])
    t.calc_splits()
    return t

class TreeHolder(object):
    '''Uses the tree attribute to provide a `tree_str` property, but also
        enables setting of the `tree_str` to update the tree.
    '''
    def __init__(self, dataset):
        self.dataset = dataset
        self.tree = None

    def get_tree_str(self):
        return self.tree.compose_newick() if self.tree else None

    def set_tree_str(self, tree_str):
        self.tree = generate_tree_with_splits_from_str(tree_str, self.dataset)

    tree_str = property(get_tree_str, set_tree_str)

    def get_tree_copy(self):
        '''Returns a deep copy of the tree instance.'''
        return generate_tree_with_splits_from_str(self.tree_str, self.dataset)

