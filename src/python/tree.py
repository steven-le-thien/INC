# import ete3
import dendropy

def query(x, node, D):
    """
        A generic query function. Given a leaf X to insert while querying internal node NODE,
        return the neighbor node X should belong in.
    """
    [(n1,l1),(n2,l2),(n3,l3)] = node.leaf_samples.items()
    v1 = D(x,l1) + D(l2, l3)
    v2 = D(x,l2) + D(l1, l3)
    v3 = D(x,l3) + D(l1, l2)
    best_v, best_n = min([(v1,n1),(v2,n2),(v3,n3)], key=lambda x:x[0])
    return best_n


class PNode(object):
    """ A generic container class with pointers to other PNodes """
    
    """ Plan: instead have self.edges_nodes, mapping PEdge objects to PNode objects
        PEdge objects will store weight, both nodes, and have a bipartitions_dict
            the bipartitions_dict maps a unique ID for each constraint tree
            to the bipartition this edge currently represents for leaves in that constraint tree
        
    """
    
    def __init__(self, val=None):
        # maps neighbor Nodes to float weight values
        self.edge_weights = dict()
        # maps neighbor nodes to leaf nodes
        self.leaf_samples = dict()
        # a value for this particular node
        self.val = val
    def __str__(self):
        if self.is_leaf():
            return "{0}".format(self.get_val())
        else:
            return "Internal {0}".format(repr(self))
    def get_val(self):
        return self.val
    
    def set_val(self, val):
        self.val = val
    
    def get_neighbors(self):
        """ Returns a mutable copy of neighbors """
        return list(self.edge_weights.keys())
    
    def get_neighbors_weights(self):
        """ Return (copy of) tuples of (neighbor, weight) tuples"""
        return list(self.edge_weights.items())
    
    def remove_neighbor(self, neighbor):
        # single direction only
        del self.edge_weights[neighbor] 
    
    def add_neighbor(self, neighbor, weight):
        # single direction only
        self.edge_weights[neighbor] = weight
    
    def is_leaf(self):
        return len(self.get_neighbors()) <= 1

    def neighbors_not(self, neighbor):
        """ Return the two neighbors that are not NEIGHBOR """
        if neighbor not in self.leaf_samples:
            raise ValueError("neighbor is not a neighbor of this node")
        r = []
        for n in self.leaf_samples:
            if n is not neighbor:
                r.append(n)
        return r
    def leaf_samples_not_for(self, neighbor):
        """ Return the two leaf samples that are not from neighbor """
        n1,n2 = self.neighbors_not(neighbor)
        return [self.leaf_samples[n1], self.leaf_samples[n2]]
    
    def add_leaf(self, neighbor, leaf, close_wt, far_wt, leaf_wt,D):
        """
            Break the edge between self and neighbor, 
                adding an internal node and a leaf hanging from the internal node
            Use close_wt, far_wt,leaf_wt values as the edge weights
            close_wt - (self, internal)
            far_wt   - (internal, neighbor)
            leaf_wt  - (internal, leaf)
            
            Also add the leaf samples for the new internal node. 
                If we break (u,v) use leaf samples u->v and v->u, along 
                with the new leaf as the samples for the new node.
                If either of u is a leaf, can be any arbitrary leaf
                on the other side of v (same for u by symmetry)
                
        """
        ## create new internal node
        internal = PNode()
        
        # remove old edges
        self.remove_neighbor(neighbor)
        neighbor.remove_neighbor(self)
        
        # add new edges
        self.add_neighbor(internal, close_wt)
        internal.add_neighbor(self, close_wt)
        neighbor.add_neighbor(internal, far_wt)
        internal.add_neighbor(neighbor, far_wt)
        leaf.add_neighbor(internal, leaf_wt)
        internal.add_neighbor(leaf, leaf_wt)
        
        
        # add leaf samples
        internal.leaf_samples[leaf] = leaf
        if self.is_leaf():
            internal.leaf_samples[self] = self
            if neighbor.is_leaf():
                internal.leaf_samples[neighbor] = neighbor
            else:
                # self is leaf, neighbor is internal.
                # use the closest one of neighbor's leaf samples
                l1,l2 = neighbor.leaf_samples_not_for(self)
                internal.leaf_samples[neighbor] = min(l1,l2, key=lambda x:D(x,leaf))
        else:
            if neighbor.is_leaf():
                internal.leaf_samples[neighbor] = neighbor
                l1,l2 = self.leaf_samples_not_for(neighbor)
                internal.leaf_samples[self] = min(l1,l2,key=lambda x:D(x,leaf))
            else:
                # both internal nodes, just use their samples
                internal.leaf_samples[self] = neighbor.leaf_samples[self]
                internal.leaf_samples[neighbor] = self.leaf_samples[neighbor]
                    
        # modify leaf_sample dicts for self and neighbor
        if not self.is_leaf():
            self.leaf_samples[internal] = self.leaf_samples[neighbor]
            del self.leaf_samples[neighbor]
            if len(self.leaf_samples) != 3:
                print(self.leaf_samples)
                raise RuntimeError()
        if not neighbor.is_leaf():
            neighbor.leaf_samples[internal] = neighbor.leaf_samples[self]
            del neighbor.leaf_samples[self]
            if len(neighbor.leaf_samples) != 3:
                print(neighbor.leaf_samples)
                raise RuntimeError()
        
        if len(internal.leaf_samples) != 3:
            print(internal.leaf_samples)
            raise RuntimeError()
            
    
    def make_newick_string(self):
        """ Return the newick string associated with this tree. 
            
            Root at an internal node.
            If only 1 leaf, root is leaf, 
            If only 2 leaves, arbitrary center in the middle"""
        newick_str = ""
        if self.is_leaf():
            # exactly one neighbor, but keep with the iterative interface anyways
            for neighbor in self.get_neighbors():
                
                if len(neighbor.get_neighbors()) > 1:
                    #root from internal node
                    return neighbor.make_newick_string()
                else:
                    # two-leaf tree
                    newick_str = "({0}:{2},{1}:{2});".format(str(self.get_val()),\
                                                   neighbor.get_val(), \
                                                   self.edge_weights[neighbor]/2)
        
            if newick_str == "":
                # single leaf tree
                newick_str = "({0}:0);".format(self.get_val())
            
        # define function that adds nodes through recursive DFS in one direction
        else:
            def make_ete_DFS(node, seen):
                final_strs = []
                for neighbor, dist in node.get_neighbors_weights():
                    # do not recursively call DFS on parent node
                    if neighbor not in seen:
                        seen.add(neighbor)
                        child_tree = make_ete_DFS(neighbor, seen)
                        final_strs.append("{0}:{1}".format(child_tree, dist))
                       
                # if final_str is still blank, node must be a leaf
                if not final_strs:
                    return str(node.get_val())
                
                return "({0})".format(','.join(final_strs))

            newick_str = '{0};'.format(make_ete_DFS(self, set([self])))
        return newick_str
    # def make_ete_tree(self):
    #     """ Return ete3 tree, rooted at this node if internal. For visualization. """
    #     return ete3.Tree(self.make_newick_string())
    
    def make_dendro_tree(self):
        """ Return dendropy tree, rooted at this node, if internal, or neighbor node if leaf. """
        return dendropy.Tree.get(data=self.make_newick_string(),schema='newick')
    
    def encode_bipartitions(self):
        """ Return a list of biparitions associated with the edges of this tree
            Functionally, I just have DendroPy do all the work
            """
        return self.make_dendro_tree().encode_bipartitions()
    
class PTree(object):
    """ A binary tree - where all internal nodes have degree 3 
        All edges are double sided.
    """
    def __init__(self, search_fn):
        self.leaves = set()
        self.root = None # an arbitrary root node. for now just the first leaf inserted
        self.search_fn = search_fn # function(start_leaf,insert_leaf, D) that returns edge to break (v1,v2)
    def num_leaves(self):
        return len(self.leaves)
    def get_leaves(self):
        return self.leaves
    def make_newick_string(self):
        return self.root.make_newick_string()
    
    # def make_ete_tree(self):
    #     """ Return ete3 tree, rooted at root node. """
    #     return self.root.make_ete_tree()
    
    def make_dendro_tree(self):
        return self.root.make_dendro_tree()
    
    def encode_bipartitions(self):
        return self.root.get_bipartitions()
    
    def get_root(self):
        return self.root
    
    def set_root(self, new_root):
        self.root = new_root
    
    def insert_leaf(self, leaf, D,start_leaf =None):
        """
            insert LEAF into this tree, using distance function D(i,j)
            use self.SEARCH_FN(start,insert,D) to find the edge to break
            then break that edge, computing the new edges' distances
        """
        # add leaf to list of leaves
        self.leaves.add(leaf)
        
        # case: first leaf added to tree
        if self.num_leaves() == 1:
            self.set_root(leaf)
        
        # case: second leaf added to tree - just have them be neighbors
        elif self.num_leaves() == 2:
            self.get_root().add_neighbor(leaf, D(leaf,self.get_root()))
            leaf.add_neighbor(self.get_root(), D(leaf,self.get_root()))
        
        # case: third leaf - have a single internal node
        elif self.num_leaves() == 3:
            root = self.get_root()
            neighbor = root.get_neighbors()[0]
            # A,B,C are root(leaf), neighbor(leaf), new_leaf
            AB = D(root, neighbor)
            BC = D(neighbor, leaf)
            AC = D(root, leaf)
            close_wt = (AC + AB - BC)/2
            far_wt = (AB + BC - AC)/2
            leaf_wt = (AC + BC - AB)/2
            root.add_leaf(neighbor,leaf,close_wt,far_wt,leaf_wt,D)
       
        # case: 4+th leaf - figure out which edge to break from search_fn and break it
        else:
            # find the new edge to break
            if not start_leaf:
                start_leaf = self.get_root()
            [v1,v2] = self.search_fn(start_leaf,leaf,D)
            
            # case 1: edge to break is internal
            if not v1.is_leaf() and not v2.is_leaf():
                # (A1,B1),(C2,D2) are leaf samples for v1, v2 respectively
                v1_other_neighbors = set(v1.get_neighbors()) - set([v2])
                v2_other_neighbors = set(v2.get_neighbors()) - set([v1])

                [A1,B1] = v1.leaf_samples_not_for(v2)
                [C2,D2] = v2.leaf_samples_not_for(v1) 

                #p is close weight, q is far weight, r is leaf weight
                r = (D(A1,leaf) + D(C2,leaf) - D(A1,C2))/2
                p = (D(A1,leaf) + D(B1,leaf) - D(A1,B1))/2 - r
                q = (D(leaf,C2) + D(leaf,D2) - D(C2,D2))/2 - r
                v1.add_leaf(v2, leaf, p, q, r,D)
            
            # case 2 - edge to break is not internal
            else:
                # a leaf hangs from edge we wish to break
                # let A be this leaf
                A = [v2,v1][v1.is_leaf()]
                other = [v1,v2][v1.is_leaf()]
                [B,C] = other.leaf_samples_not_for(A)
                #p is close weight, q is far weight, r is leaf weight
                # variables i_j represent the algebraic quantity "i+j"
                r_q = (D(B,leaf) + D(C,leaf) - D(B,C))/2
                p_r = D(A,leaf)
                p_q = (D(A,B) + D(A,C) - D(B,C))/2
                
                p = (p_r - r_q + p_q)/2
                r = p_r - p
                q = p_q - p
                A.add_leaf(other, leaf, p,q,r,D)
                

