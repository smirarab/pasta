# uym2 added
# June 2017
# utils for tree decomposition


from dendropy import Tree, Node
try:
    from queue import Queue # python 3
except:
    from Queue import Queue # python 2
#from tree import PhylogeneticTree
from sepp import get_logger

_LOG = get_logger(__name__)

def decompose_by_diameter(a_tree,strategy,max_size=None,min_size=None,max_diam=None):
    def __ini_record__():
        for node in a_tree.postorder_node_iter():
               __updateNode__(node)
    
    def __find_midpoint_edge__(t):
        u = t.seed_node.bestLCA.anchor
        d = 0
        while (d+u.edge_length < t.seed_node.diameter/2):
            d += u.edge_length
            u = u.parent_node
        return u.edge
    
    def __find_centroid_edge__(t):
        u = t.seed_node
        product = 0
        acc_nleaf = 0

        while not u.is_leaf():
            max_child = None
            max_child_nleaf = 0
            for ch in u.child_node_iter():
                if ch.nleaf > max_child_nleaf:
                    max_child_nleaf = ch.nleaf
                    max_child = ch
            acc_nleaf += (u.nleaf-max_child.nleaf)
            new_product = max_child.nleaf * acc_nleaf
            if new_product <= product:
                break
            product = new_product
            u = max_child

        return u.edge

    def __bisect__(t,e):
#        e = __find_centroid_edge__(t)
        
        u = e.tail_node
        v = e.head_node

        u.remove_child(v)
        t1 = Tree(seed_node = v)

        if u.num_child_nodes() == 1:
            p = u.parent_node
            v = u.child_nodes()[0]
            l_v = v.edge_length
            u.remove_child(v)
            if p is None: # u is the seed_node; this means the tree runs out of all but one side
                t.seed_node = v
                return t,t1
            l_u = u.edge_length
            p.remove_child(u)
            p.add_child(v)
            v.edge_length = l_u+l_v
            u = p

        while u is not None:
            __updateNode__(u)
            u = u.parent_node

        return t,t1

    def __clean_up__(t):
        for node in t.postorder_node_iter():
            delattr(node,"nleaf")
            delattr(node,"anchor")
#            delattr(node,"maxheight")
            delattr(node,"maxdepth")
            delattr(node,"diameter")
#            delattr(node,"topo_diam")
            delattr(node,"bestLCA")

    def __updateNode__(node):
        if node.is_leaf():
            node.anchor = node
#            node.maxheight = 0
            node.maxdepth = 0
            node.diameter = 0
#            node.topo_diam = 0
            node.bestLCA = node
            node.nleaf = 1
            return

#        n1 = -1
#        n2 = -1
        d1 = -1
        d2 = -1
        anchor1 = None
        anchor2 = None
        node.diameter = 0
#        node.topo_diam = 0
        node.bestLCA = None
        node.nleaf = 0

        for ch in node.child_node_iter():
               node.nleaf += ch.nleaf
#               n = ch.maxheight + 1
               d = ch.maxdepth + ch.edge_length
#               if n > n1:
#                   n2 = n1
#                   n1 = n
#                   anchor2 = anchor1
#                   anchor1 = ch.anchor
#               elif n > n2:
#                   n2 = n
#                   anchor2 = ch.anchor
               if d > d1:
                   d2 = d1
                   d1 = d
                   anchor2 = anchor1
                   anchor1 = ch.anchor
               elif d > d2:
                   d2 = d
                   anchor2 = ch.anchor
               if ch.diameter > node.diameter:
                   node.diameter = ch.diameter
                   node.bestLCA = ch.bestLCA
#               node.diameter = max(ch.diameter,node.diameter)

#        node.diameter = max(d1+d2, node.diameter)
        node.maxdepth = d1
#        node.maxheight = n1
        node.anchor = anchor1
        if d1+d2 > node.diameter:
            node.diameter = d1+d2
            node.bestLCA = node

    def __get_breaking_edge__(t,edge_type):
        if t.seed_node.nleaf <= max_size and t.seed_node.diameter <= max_diam:
            return None
        if edge_type == 'midpoint':
            e = __find_midpoint_edge__(t)
        elif edge_type == 'centroid':
            e = __find_centroid_edge__(t)
        else:
            _LOG.warning("Invalid decomposition type! Please use either 'midpoint' or 'centroid'")
            return None

        n = e.head_node.nleaf
        if (n < min_size) or (t.seed_node.nleaf - n) < min_size:
            return None
        return e

    def __check_stop__(t):
        return ( (t.seed_node.nleaf <= max_size and t.seed_node.diameter <= max_diam) or
                 (t.seed_node.nleaf//2 < min_size) )     

    def __break_by_MP_centroid__(t):
        e = __get_breaking_edge__(t,'midpoint')
        if e is None:
#            print("Midpoint failed. Trying centroid decomposition...")
            e = __get_breaking_edge__(t,'centroid')
#        else:
#            print("Successfully splitted by midpoint")
        return e

    def __break(t):
        if strategy == "centroid":
            return __get_breaking_edge__(t,'centroid')
        elif strategy == "midpoint":
            return __break_by_MP_centroid__(t)
        else:
            raise Exception("strategy not valid: %s" %strategy)

    _LOG.debug("Starting brlen decomposition ...")
    tqueue = Queue()
    __ini_record__()
    min_size = min_size if min_size else 0
    max_size = max_size if max_size else a_tree.seed_node.nleaf
    max_diam = max_diam if max_diam else a_tree.seed_node.diameter

    # try using midpoint
    e = __break(a_tree)

    if e is None:
        __clean_up__(a_tree)
        return [a_tree]
        
    treeMap = [] 
    tqueue.put((a_tree,e))
    while not tqueue.empty():
        t,e = tqueue.get()
        t1,t2 = __bisect__(t,e)
        e1 = __break(t1)
        if e1 is None:
             __clean_up__(t1)            
             treeMap.append(t1)
        else:
            tqueue.put((t1,e1))
        e2 = __break(t2)
        if e2 is None:
             __clean_up__(t2)
             treeMap.append(t2)
        else:
            tqueue.put((t2,e2))

    return treeMap

def midpoint_bisect(tree):
######## working ##########
