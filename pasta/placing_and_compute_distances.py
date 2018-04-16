#! /usr/bin/env python

from dendropy import Tree
from decompose_lib import decompose_by_diameter, compute_group_distance_matrix,place_group_onto_tree 
import sys
import os
from pasta import get_logger
_LOG = get_logger(__name__)


intree_file = sys.argv[1]
grouping_file = sys.argv[2]
nleaf_file = sys.argv[3]
distance_file = sys.argv[4]

t = Tree.get_from_path(intree_file,'newick')

grouping = {}
with open(grouping_file,'r') as f:
    for line in f:
        name, taxon = line.split()
        grouping[taxon] = name
_LOG.info('computing treeMap ... ')
treeMap = place_group_onto_tree(t,grouping)

D = compute_group_distance_matrix(t,treeMap)


with open(distance_file,'w') as f:
    for A,B in D:
        f.write(A + " " + B + " " + str(D[(A,B)]) + "\n")

with open(nleaf_file,'w') as f:
    for name in sorted(treeMap):
        node = treeMap[name]
        f.write(name +  " " + str(node.nleaf) + "\n")
