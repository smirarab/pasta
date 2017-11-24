#! /usr/bin/env python

from dendropy import Tree
#from tree import PhylogeneticTree
from decompose_lib import decompose_by_diameter, compute_group_distance_matrix,place_group_onto_tree 
import sys
import os

intree_file = sys.argv[1]
grouping_file = sys.argv[2]
nleaf_file = sys.argv[3]
distance_file = sys.argv[4]
#min_size = int(sys.argv[4])
#max_diam = float(sys.argv[5])

t = Tree.get_from_path(intree_file,'newick')
#T = PhylogeneticTree(t)

grouping = {}
with open(grouping_file,'r') as f:
    for line in f:
        name, taxon = line.split()
        grouping[taxon] = name
#print(grouping)		
print('computing treeMap ... ')
treeMap = place_group_onto_tree(t,grouping)
#treeMap = T.decompose_tree(max_size,'centroid',decomp_strategy="brlen")#,minSize=min_size,maxDiam=max_diam)
#treeMap = decompose_by_diameter(t,'centroid',max_size=max_size)
#treeNum = len(treeMap.keys())
#print(treeMap)
#print('writing tree ... ')


D = compute_group_distance_matrix(t,treeMap)

#print(D)

with open(distance_file,'w') as f:
    for A,B in D:
        f.write(A + " " + B + " " + str(D[(A,B)]) + "\n")

with open(nleaf_file,'w') as f:
    for name in sorted(treeMap):
        node = treeMap[name]
        f.write(name +  " " + str(node.nleaf) + "\n")
	#print(name + " " + str(node.nleaf) + "\n")
