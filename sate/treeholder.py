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
from sate import get_logger
from dendropy.dataobject.taxon import Taxon
from dendropy.dataio.newick import tree_source_iter
from dendropy.dataobject.tree import Tree
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
    _LOG.debug("polytomies resolution - updating splits")
    if update_splits:
        tree.update_splits()
    _LOG.debug("polytomies resolved.")
    
def check_taxon_labels(taxon_set, dataset):
    ts = set([i for tset in dataset.taxon_sets for i in tset.labels()])
    ds = set(taxon_set.labels())
    extra = ds - ts
    missing = ts - ds
    return extra, missing

from dendropy.dataio.nexustokenizer import RootingInterpreter, StrToTaxon
from dendropy import dataobject
from dendropy.utility import containers

def stt_require_taxon(stt,label):
        t = stt.get_taxon(label)
        if t is not None:
            return t
        t = dataobject.Taxon(label=label)
        stt.taxon_set.append(t)
        stt.label_taxon[label] = t
        return stt._returning(t, label)

def read_newick_with_translate(stream,translate_dict):
    """
    Instantiates and returns a `DataSet` object based on the
    NEWICK-formatted contents read from the file-like object source
    `stream`.
    """
    ts = [t for t in tree_source_iter(stream=stream,
            translate_dict=translate_dict,
            allow_repeated_use=True)]
    return ts[0]


def tree_from_token_stream(stream_tokenizer, **kwargs):
    """
    Processes a (SINGLE) TREE statement. Assumes that the input stream is
    located at the beginning of the statement (i.e., the first non-comment
    token should be the opening parenthesis of the tree definition).

    str_to_taxon kwarg (if used) must supply the StrToTaxon interface).
    """
    translate_dict = kwargs.get("translate_dict", None)
    encode_splits = kwargs.get("encode_splits", False)
    rooting_interpreter = kwargs.get("rooting_interpreter", RootingInterpreter(**kwargs))
    finish_node_func = kwargs.get("finish_node_func", None)
    edge_len_type = kwargs.get("edge_len_type", float)
    taxon_set = kwargs.get("taxon_set", None)
    suppress_internal_node_taxa = kwargs.get("suppress_internal_node_taxa", False)
    store_tree_weights = kwargs.get("store_tree_weights", False)
    extract_comment_metadata = kwargs.get('extract_comment_metadata', False)
    case_sensitive_taxon_labels = kwargs.get('case_sensitive_taxon_labels', False)
    allow_repeated_use = kwargs.get('allow_repeated_use', False)
    stream_tokenizer_extract_comment_metadata_setting = stream_tokenizer.extract_comment_metadata
    stream_tokenizer.extract_comment_metadata = extract_comment_metadata
    if taxon_set is None:
        taxon_set = dataobject.TaxonSet()
    tree = dataobject.Tree(taxon_set=taxon_set)

    stream_tokenizer.tree_rooting_comment = None # clear previous comment
    stream_tokenizer.clear_comment_metadata()
    token = stream_tokenizer.read_next_token()
    if not token:
        return None
    tree.is_rooted = rooting_interpreter.interpret_as_rooted(stream_tokenizer.tree_rooting_comment)
#    if stream_tokenizer.tree_rooting_comment is not None:
#        tree.is_rooted = rooting_interpreter.interpret_as_rooted(stream_tokenizer.tree_rooting_comment)
#    elif rooting_interpreter.interpret_as_rooted(stream_tokenizer.tree_rooting_comment):
#        tree_is_rooted = True

    if store_tree_weights and stream_tokenizer.tree_weight_comment is not None:
        try:
            weight_expression = stream_tokenizer.tree_weight_comment.split(' ')[1]
            tree.weight = eval("/".join(["float(%s)" % cv for cv in weight_expression.split('/')]))
        except IndexError:
            pass
        except ValueError:
            pass
        stream_tokenizer.tree_weight_comment = None

    if encode_splits:
        if len(taxon_set) == 0:
            raise Exception("When encoding splits on a tree as it is being parsed, a "
                + "fully pre-populated TaxonSet object must be specified using the 'taxon_set' keyword " \
                + "to avoid taxon/split bitmask values changing as new Taxon objects are created " \
                + "and added to the TaxonSet.")
        if tree.is_rooted:
            tree.split_edges = {}
        else:
            atb = taxon_set.all_taxa_bitmask()
            d = containers.NormalizedBitmaskDict(mask=atb)
            tree.split_edges = d
        split_map = tree.split_edges

    stt = kwargs.get('str_to_taxon')
    if stt is None:
        stt = StrToTaxon(taxon_set,
                translate_dict,
                allow_repeated_use=allow_repeated_use,
                case_sensitive=case_sensitive_taxon_labels)

    tree.seed_node = dataobject.Node()
    curr_node = tree.seed_node
    if encode_splits:
        curr_node.edge.split_bitmask = 0L

    ### NHX format support ###
    def store_node_comments(active_node):
        if stream_tokenizer.comments:
            active_node.comments.extend(stream_tokenizer.comments)

    def store_comment_metadata(target):
        if extract_comment_metadata:
            if stream_tokenizer.has_comment_metadata():
                comment_metadata = stream_tokenizer.comment_metadata
                try:
                    target.comment_metadata.update(comment_metadata)
                except AttributeError:
                    target.comment_metadata = comment_metadata
                stream_tokenizer.clear_comment_metadata()
            elif not hasattr(target, "comment_metadata"):
                target.comment_metadata = {}

    # store and clear comments
    tree.comments = stream_tokenizer.comments
    stream_tokenizer.clear_comments()
    store_comment_metadata(tree)

    while True:
        if not token or token == ';':
            if curr_node is not tree.seed_node:
                raise stream_tokenizer.data_format_error("Unbalanced parentheses -- not enough ')' characters found in tree description")
            if encode_splits:
                split_map[curr_node.edge.split_bitmask] = curr_node.edge
            break
        if token == '(':
            if not curr_node.parent_node:
                if curr_node.child_nodes():
                    raise stream_tokenizer.data_format_error("Unexpected '(' after the tree description.  Expecting a label for the root or a ;")
            tmp_node = dataobject.Node()
            if encode_splits:
                tmp_node.edge.split_bitmask = 0L
            curr_node.add_child(tmp_node)
            curr_node = tmp_node
            token = stream_tokenizer.read_next_token()
            store_node_comments(curr_node)
            store_comment_metadata(curr_node)
        elif token == ',':
            tmp_node = dataobject.Node()
            if curr_node.is_leaf() and not curr_node.taxon:
#                 curr_node.taxon = taxon_set.Taxon(oid="UNAMED_" + str(id(curr_node)), label='')
#                 taxon_set.add(curr_node.taxon)
                raise stream_tokenizer.data_format_error("Missing taxon specifier in a tree -- found either a '(,' or ',,' construct.")
            p = curr_node.parent_node
            if not p:
                raise stream_tokenizer.data_format_error("Comma found one the 'outside' of a newick tree description")
            if encode_splits:
                tmp_node.edge.split_bitmask = 0L
                e = curr_node.edge
                u = e.split_bitmask
                split_map[u] = e
                p.edge.split_bitmask |= u
            if finish_node_func is not None:
                finish_node_func(curr_node, tree)
            p.add_child(tmp_node)
            curr_node = tmp_node
            token = stream_tokenizer.read_next_token()
            store_node_comments(curr_node)
            store_comment_metadata(curr_node)
        else:
            if token == ')':
                if curr_node.is_leaf() and not curr_node.taxon:
                    raise stream_tokenizer.data_format_error("Missing taxon specifier in a tree -- found either a '(,' or ',,' construct.")
                p = curr_node.parent_node
                if not p:
                    raise stream_tokenizer.data_format_error("Unbalanced parentheses -- too many ')' characters found in tree description")
                if encode_splits:
                    e = curr_node.edge
                    u = e.split_bitmask
                    p.edge.split_bitmask |= u
                    split_map[u] = curr_node.edge
                if finish_node_func is not None:
                    finish_node_func(curr_node, tree)
                curr_node = p
            else:
                is_leaf = curr_node.is_leaf()
                if is_leaf:
                    if curr_node.taxon:
                        raise stream_tokenizer.data_format_error("Multiple labels found for the same leaf (taxon '%s' and label '%s')" % (str(curr_node.taxon), token))
                    try:
                        t = stt_require_taxon(stt,label=token)
                    except StrToTaxon.MultipleTaxonUseError, e:
                        raise stream_tokenizer.data_format_error(e.msg)
                else:
                    if curr_node.label:
                        raise stream_tokenizer.data_format_error("Multiple labels found for the same leaf (taxon '%s' and label '%s')" % (curr_node.label, token))
                    if suppress_internal_node_taxa:
                        t = None
                    else:
                        try:
                            t = stt.get_taxon(label=token)
                        except StrToTaxon.MultipleTaxonUseError, e:
                            raise stream_tokenizer.data_format_error(e.msg)
                if t is None:
                    curr_node.label = token
                else:
                    curr_node.taxon = t
                    if encode_splits:
                        try:
                            cm = t.split_bitmask
                        except:
                            cm = 1 << (stt.index(t))
                        e = curr_node.edge
                        e.split_bitmask = cm
                        split_map[cm] = e

            token = stream_tokenizer.read_next_token()
            store_node_comments(curr_node)
            store_comment_metadata(curr_node)
            if token == ':':
                edge_length_str = stream_tokenizer.read_next_token(ignore_punctuation='-+.')
                store_node_comments(curr_node)
                store_comment_metadata(curr_node)
                if not edge_length_str:
                    raise stream_tokenizer.data_format_error("Expecting a branch length after : but encountered the end of the tree description" )
                try:
                    curr_node.edge.length = edge_len_type(edge_length_str)
                except:
                    curr_node.edge.length = edge_length_str
                token = stream_tokenizer.read_next_token()
                store_node_comments(curr_node)
                store_comment_metadata(curr_node)
    stream_tokenizer.extract_comment_metadata = stream_tokenizer_extract_comment_metadata_setting
    return tree

dendropy.dataio.nexustokenizer.tree_from_token_stream = tree_from_token_stream
 
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

