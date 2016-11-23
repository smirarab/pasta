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

from cStringIO import StringIO
from pasta.tree import PhylogeneticTree
from pasta.errors import TaxaLabelsMismatchError
from pasta import get_logger
from dendropy import Tree, TreeList
_LOG = get_logger(__name__)

# Provide a random number generator
import random

# StrToTaxon class from Dendropy 3
from dendropy.utility.error import DataParseError
class StrToTaxon(object):

    class MultipleTaxonUseError(DataParseError):
        def __init__(self, *args, **kwargs):
            DataParseError.__init__(self, *args, **kwargs)

    def __init__(self,
            taxon_namespace,
            translate_dict=None,
            allow_repeated_use=False,
            case_sensitive=True):
        """
        __init__ creates a StrToTaxon object with the requested policy of taxon
        repitition.
        If `allow_repeated_use` is True, then get_taxon and require_taxon
        can be called multiple times to get the same taxon.  If it is false then
        calling the functions with the same label will generate a DataParseError
        indicating that the taxon has been used multiple times."""
        self.taxon_namespace = taxon_namespace
        self.case_sensitive = case_sensitive
        if translate_dict is not None:
            self.translate = translate_dict
        else:
            if self.case_sensitive:
                self.translate = containers.OrderedCaselessDict()
            else:
                self.translate = {}
        if self.case_sensitive:
            self.label_taxon = {}
        else:
            self.label_taxon = containers.OrderedCaselessDict()
        for t in self.taxon_namespace:
            self.label_taxon[t.label] = t
        if allow_repeated_use:
            self.returned_taxon_namespace = None
        else:
            self.returned_taxon_namespace = set()

    def _returning(self, t, label):
        #_LOG.debug("Checking if we can return %s" % str(t))
        if (self.returned_taxon_namespace is not None) and (t is not None):
            if t in self.returned_taxon_namespace:
                raise StrToTaxon.MultipleTaxonUseError(message="Taxon %s used twice (it appears as %s the second time)" % (str(t), label))
            self.returned_taxon_namespace.add(t)
        return t

    def get_taxon(self, label):
        t = self.translate.get(label)
        if t is None:
            t = self.label_taxon.get(label)
        if t is None:
            for taxon in self.taxon_namespace:
                if taxon.label == label:
                    t = taxon
                    break
        if t is not None:
            self.label_taxon[label] = t
            return self._returning(t, label)
        else:
            return None

    def require_taxon(self, label):
        t = self.get_taxon(label)
        if t is not None:
            return t
        t = Taxon(label=label)
        self.taxon_namespace.add(t)
        self.label_taxon[label] = t
        return self._returning(t, label)

    def index(self, t):
        return self.taxon_namespace.index(t)

# RootingInterpreter class from Dendropy 3
class RootingInterpreter(object):

    def evaluate_as_rooted_kwargs(kwdict, default):
        if "as_rooted" in kwdict and "as_unrooted" in kwdict \
                and (kwdict["as_rooted"] != (not kwdict["as_unrooted"])):
            raise TypeError("Conflict rooting specification: 'as_rooted'=%s and 'as_unrooted'=%s" \
                    % (kwdict["as_rooted"], kwdict["as_unrooted"]))
        if "as_rooted" in kwdict:
            if kwdict["as_rooted"] is None:
                return default
            if kwdict["as_rooted"] is True or kwdict["as_rooted"] is False:
                return kwdict["as_rooted"]
            else:
                raise ValueError("Invalid value for 'as_rooted' (expecting True/False, but received '%s')" \
                        % kwdict["as_rooted"])
        elif "as_unrooted" in kwdict:
            if kwdict["as_unrooted"] is None:
                return default
            if kwdict["as_unrooted"] is True or kwdict["as_unrooted"] is False:
                return not kwdict["as_unrooted"]
            else:
                raise ValueError("Invalid value for 'as_unrooted' (expecting True/False, but received '%s')" \
                        % kwdict["as_rooted"])
        else:
            return default

    evaluate_as_rooted_kwargs = staticmethod(evaluate_as_rooted_kwargs)

    def evaluate_default_as_rooted_kwargs(kwdict, default=None):
        if "default_as_rooted" in kwdict and "default_as_unrooted" in kwdict \
                and (kwdict["default_as_rooted"] != (not kwdict["default_as_unrooted"])):
            raise TypeError("Conflict rooting specification: 'default_as_rooted'=%s and 'default_as_unrooted'=%s" \
                    % (kwdict["default_as_rooted"], kwdict["default_as_unrooted"]))
        if "default_as_rooted" in kwdict:
            if kwdict["default_as_rooted"] is True or kwdict["default_as_rooted"] is False:
                return kwdict["default_as_rooted"]
            else:
                raise ValueError("Invalid value for 'default_as_rooted' (expecting True/False, but received '%s')" \
                        % kwdict["default_as_rooted"])
        elif "default_as_unrooted" in kwdict:
            if kwdict["default_as_unrooted"] is True or kwdict["default_as_unrooted"] is False:
                return not kwdict["default_as_unrooted"]
            else:
                raise ValueError("Invalid value for 'default_as_unrooted' (expecting True/False, but received '%s')" \
                        % kwdict["default_as_rooted"])
        else:
            return default

    evaluate_default_as_rooted_kwargs = staticmethod(evaluate_default_as_rooted_kwargs)

    def __init__(self, **kwargs):
        self._as_rooted = RootingInterpreter.evaluate_as_rooted_kwargs(kwargs, None)
        self._default_as_rooted = RootingInterpreter.evaluate_default_as_rooted_kwargs(kwargs, False)

    def update(self, **kwargs):
        self._as_rooted = RootingInterpreter.evaluate_as_rooted_kwargs(kwargs, self._as_rooted)
        self._default_as_rooted = RootingInterpreter.evaluate_default_as_rooted_kwargs(kwargs, self._default_as_rooted)

    def _get_as_rooted(self):
        return self._as_rooted

    def _set_as_rooted(self, val):
        self._as_rooted = val

    as_rooted = property(_get_as_rooted, _set_as_rooted)

    def _get_as_unrooted(self):
        return not self._as_rooted

    def _set_as_unrooted(self, val):
        self._as_rooted = not val

    as_unrooted = property(_get_as_unrooted, _set_as_unrooted)

    def _get_default_as_rooted(self):
        return self._default_as_rooted

    def _set_default_as_rooted(self, val):
        self._default_as_rooted = val

    default_as_rooted = property(_get_default_as_rooted, _set_default_as_rooted)

    def _get_default_as_unrooted(self):
        return not self._default_as_rooted

    def _set_default_as_unrooted(self, val):
        self._default_as_rooted = not val

    default_as_unrooted = property(_get_default_as_unrooted, _set_default_as_unrooted)

    def interpret_as_rooted(self, tree_rooting_comment=None, **kwargs):
        if self.as_rooted is not None:
            return self.as_rooted
        elif tree_rooting_comment is not None:
            return tree_rooting_comment.upper() == "&R"
        else:
            return self.default_as_rooted

    def interpret_as_unrooted(self, **kwargs):
        return not self.interpret_as_rooted(**kwargs)

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
    
def check_taxon_labels(taxon_namespace, dataset):
    ts = set([i for tset in dataset.taxon_namespaces for i in tset.labels()])
    ds = set(taxon_namespace.labels())
    extra = ds - ts
    missing = ts - ds
    return extra, missing

from dendropy.utility import container

def stt_require_taxon(stt,label):
        t = stt.get_taxon(label)
        if t is not None:
            return t
        print "Adding taxon", label
        t = Taxon(label=label)
        stt.taxon_namespace.append(t)
        stt.label_taxon[label] = t
        return stt._returning(t, label)

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
    taxon_namespace = kwargs.get("taxon_namespace", None)
    suppress_internal_node_taxa = kwargs.get("suppress_internal_node_taxa", False)
    store_tree_weights = kwargs.get("store_tree_weights", False)
    extract_comment_metadata = kwargs.get('extract_comment_metadata', False)
    case_sensitive_taxon_labels = kwargs.get('case_sensitive_taxon_labels', False)
    allow_repeated_use = kwargs.get('allow_repeated_use', False)
    stream_tokenizer_extract_comment_metadata_setting = stream_tokenizer.extract_comment_metadata
    stream_tokenizer.extract_comment_metadata = extract_comment_metadata
    if taxon_namespace is None:
        taxon_namespace = TaxonSet()
    tree = Tree(taxon_namespace=taxon_namespace)

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
        if len(taxon_namespace) == 0:
            raise Exception("When encoding splits on a tree as it is being parsed, a "
                + "fully pre-populated TaxonSet object must be specified using the 'taxon_namespace' keyword " \
                + "to avoid taxon/split bitmask values changing as new Taxon objects are created " \
                + "and added to the TaxonSet.")
        if tree.is_rooted:
            tree.split_edges = {}
        else:
            atb = taxon_namespace.all_taxa_bitmask()
            d = containers.NormalizedBitmaskDict(mask=atb)
            tree.split_edges = d
        split_map = tree.split_edges

    stt = kwargs.get('str_to_taxon')
    if stt is None:
        stt = StrToTaxon(taxon_namespace,
                translate_dict,
                allow_repeated_use=allow_repeated_use,
                case_sensitive=case_sensitive_taxon_labels)

    tree.seed_node = Node()
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
            tmp_node = Node()
            if encode_splits:
                tmp_node.edge.split_bitmask = 0L
            curr_node.add_child(tmp_node)
            curr_node = tmp_node
            token = stream_tokenizer.read_next_token()
            store_node_comments(curr_node)
            store_comment_metadata(curr_node)
        elif token == ',':
            tmp_node = Node()
            if curr_node.is_leaf() and not curr_node.taxon:
#                 curr_node.taxon = taxon_namespace.Taxon(oid="UNAMED_" + str(id(curr_node)), label='')
#                 taxon_namespace.add(curr_node.taxon)
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
    from dendropy.legacy.treesplit import delete_outdegree_one
    #delete_outdegree_one(tree_list[0])
    tree_list[0].suppress_unifurcations()
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

