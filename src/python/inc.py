from tree import PNode, PTree
import queries
from mst import prims_mst
from dist import mat_to_func, parse_distance_matrix
from nj import NJ
import numpy as np
import random
import dendropy
from accuracy import satisfies_constraint
import argparse

def voting_dfs(insert_leaf, D, query_func, q, weighted, 
               start_leaf=None, start_edge=None, end_edge=None, edge_candidates=None,
               warnings=False):
    """ 
        start_leaf - a pointer to any arbitrary leaf in the tree
        D - a function mapping pairs of leaves to their D-distance
        q - if not None, use only valid queries (where diameter <= q)
            else use all queries
        weighted - if True, weight queries as 1/diam**2. else give only one vote
        
        return - a list of edges with maximal vote weight
    """
    msg = "start_leaf={0}, start_edge={1}, end_edge={2} bad arguments".format(
            start_leaf, start_edge, end_edge)
    assert start_leaf or (start_edge and end_edge), msg


    insert_taxon = insert_leaf.get_val()
    # maps edges (u,v) to (relative) scores
    score_of = dict() 
    # use DFS search, starting from the edge subtending START_LEAF 
    score = 0
    # INVARIANT: all edges are ordered to be (n1,n2) 
    # where n1 is closer to START_LEAF than n2
    if start_edge and end_edge:
        first_edge = start_edge
    else:
        first_edge = (start_leaf, start_leaf.get_neighbors()[0])
    stack = [(score,first_edge)]
    
    exists_valid_query = True
    if q and not weighted:
        # in round 1, check for condition where no valid query exists
        exists_valid_query = False
    i = 0
    while stack:
        score, edge = stack.pop()
        i += 1
        score_of[edge] = score
        u,v = edge
        if v.is_leaf() or edge == end_edge:
            continue
            
        # two children vertices of the deeper node v
        c1, c2 = v.neighbors_not(u)
        
        # determine the diameter
        l1,l2,l3 = v.leaf_samples.values()
        l1,l2,l3 = l1.get_val(), l2.get_val(),l3.get_val()
        diameter, l_a, l_b = max((D(l1,l2),l1,l2),
                                 (D(l1,l3),l1,l3),
                                 (D(l1,insert_taxon),l1,insert_taxon),
                                 (D(l2,l3),l2,l3),
                                 (D(l2,insert_taxon),l2,insert_taxon), 
                                 (D(l3,insert_taxon),l3,insert_taxon)
                                 )
        # if q is None or nonpositive, any query is valid
        vote_wt = (not q) or (q <= 0) or (diameter <= q)
        # reweight if weighted, else keep as-is
        vote_wt *= diameter**(-16.0 * weighted)
        
        if not vote_wt:
            # don't query this node as we will ignore it.
            stack.append((score, (v, c1)))
            stack.append((score, (v, c2)))
            continue
        elif q and not weighted:
            # in round 1, check for condition where no valid query exists
            exists_valid_query = True
        
        best_neighbor = query_func(insert_leaf, v, D)
        if best_neighbor is u:
            stack.append((score - vote_wt, (v, c1)))
            stack.append((score - vote_wt, (v, c2)))
        elif best_neighbor is c1:
            stack.append((score + vote_wt, (v, c1)))
            stack.append((score, (v, c2)))
        elif best_neighbor is c2:
            stack.append((score, (v, c1)))
            stack.append((score + vote_wt, (v, c2)))
        else:
            raise RuntimeError("query_func did not return a neighbor internal_node")
    
    # find all candidate edges with the maximum score
    # at this point, directedness is not important, need to treat edge like an unordered pair
    max_score = float("-inf")
    max_edges = set()
    for edge,score in score_of.items():
        edge = frozenset(edge)
        if edge_candidates and edge not in edge_candidates:
            continue
        if score > max_score:
            max_edges = set([edge])
            max_score = score
        elif score == max_score:
            max_edges.add(edge)
        else: # do nothing
            pass
    
    if warnings and q and not exists_valid_query:
        import warnings
        msg = "No valid query exists.\nTree:\n{0}.\nLeaf:\n{1}\nMatrix:\n{2}\n"
        msg = msg.format(first_edge[0].make_newick_string(), 
                         insert_leaf.get_val(), 
                         D(get_mat=True))
        # print(msg)
        warnings.warn(msg)

    return max_edges

def random_dfs(start_leaf=None, start_edge=None, end_edge=None):
    """
        Given a start_leaf into the tree, pick one edge an random
    """
    msg = "start_leaf={0}, start_edge={1}, end_edge={2} bad arguments".format(
            start_leaf, start_edge, end_edge)
    assert start_leaf or (start_edge and end_edge), msg

    best_score, best_edge = float("-inf"), None
    first_score = random.random()
    if start_edge and end_edge:
        first_edge = start_edge
    else:
        first_edge = (start_leaf, start_leaf.get_neighbors()[0])
    stack = [(first_score,first_edge)]
    while stack:
        score, edge = stack.pop()
        best_score, best_edge = max((best_score, best_edge), (score, edge))

        u,v = edge
        if v.is_leaf() or edge == end_edge:
            continue
            
        # two children vertices of the deeper node v
        c1, c2 = v.neighbors_not(u)
        stack.append((random.random(), (v, c1)))
        stack.append((random.random(), (v, c2)))
    return best_edge


def voting_with(query_func, q, r1=True, r2=True, r3=True, warnings=False):
    """ 
    Return a voting search function that uses QUERY_FUNC to query internal nodes,
    returning, the edge that wins the most votes (or by some tie breaking scheme)
    interface:
        QUERY_FUNC(insert_leaf, internal_node, D)
    """
    
    # The voting goes in at most three rounds, and fewer if no ties in a previous round
    # In each round, we query all nodes and ask for their votes using DFS over the edges
    # If one round ends in a tie, the remaining edges go to the next round, 
    #     while others are permanently eliminated
    # Round 1: all valid queries get 1 vote. invalid get 0
    # Round 2: all valid queries get 1/diam**2 vote weight. invalid get 0
    # Round 3: all (valid and invalid) queries get 1/diam**2 vote weight
    # Round 4: pick uniformly randomly among the remaining edges.
    # note: diam is maximum pairwise D-distance of (l1,l2,l3,insert_leaf)
        
    
    def voting_search(insert_leaf, D, start_leaf=None, start_edge=None, end_edge=None):
        """ 
            args:
                start_leaf - a pointer to any arbitrary leaf in the tree
                insert_leaf - the leaf to be inserted
                D - a function mapping pairs of leaves to their D-distance
            kwargs:
                start_edge -    for constrained_inc
                                directed edge specifying direction to start voting.
                                candidates do include START_EDGE
                                queries should start with START_EDGE[1]

                end_edge -      for constrained_inc
                                last valid edge to vote on. no candidates are considered
                                past end_edge.
            
            return - the edge to break to insert insert_leaf, 
                        determined by 4 rounds of voting with tie-breaking
        """
        msg = "start_leaf={0}, start_edge={1}, end_edge={2} bad arguments".format(
                start_leaf, start_edge, end_edge)
        assert start_leaf or (start_edge and end_edge), msg

        # round 1
        candidates = []
        if r1:
            candidates = voting_dfs(insert_leaf, D, query_func, q, False, 
                                start_leaf, start_edge, end_edge, candidates, warnings=warnings)
        # round 2
        if r2 and len(candidates) != 1:
            candidates = voting_dfs(insert_leaf, D, query_func, q, True,
                                    start_leaf, start_edge, end_edge, candidates, warnings=warnings)
        # round 3
        if r3 and len(candidates) != 1:
            candidates = voting_dfs(insert_leaf, D, query_func, None, True,
                                    start_leaf, start_edge, end_edge, candidates, warnings=warnings)
        
        # if r1 and r2 and r3 are not selected, candidates is still empty - pick a winner at random
        if not candidates:
            candidates = [random_dfs(start_leaf, start_edge, end_edge)]

        # round 4
        if len(candidates) > 1:
            # print("RANDOM ROUND 4 with {0} candidates {1}".format(len(candidates), candidates))
            candidates = [random.choice(list(candidates))]

        # one candidate left
        assert len(candidates) == 1, "VOTING FAILURE"
        winner = candidates.pop()
        return winner
    
    return voting_search

def clique_distance(clique, vertex, D, q):
    """
        Given a clique, an external vertex, and distance metric D
        return the maximum pairwise distance of vertex to all members in the clique
        If the distance s at least Q, return infinity
    """

    # Clique distance is the max distance to any clique member
    # dist = float("-inf")
    # for member in clique:
    #     dist = max(dist, D(member,vertex))
    #     if dist >= q:
    #         return float("inf")
    # return dist
    
    # Clique distance in the min distance to any clique member
    # dist = float("inf")
    # for member in clique:
    #     dist = min(dist, D(member,vertex))
    #     if dist >= q:
    #         return float("inf")
    # return dist

    # Clique distance is the average distance to all clique members
    total_dist = 0.0
    for member in clique:
        d = D(member,vertex)
        if d >= q:
            return float("inf")
        total_dist += d
    return total_dist / len(clique)


def partition(names, D, q):
    """
        Given NAMES of length N, distance function D, and threshold Q,
        partition NAMES into disjoint cliques with edges at most Q,
        where no partition is larger than sqrt(N)
        
        Things to play with 
        It's arbitrary greedy search

        what is the "closest" item to a clique? 
            see implementation in clique_distance
    """
    N = len(names)
    names = set(names) # copy the list of names
    
    partitions_list = [] # list of sets of names
    while names:
        curr_partition = set()
        curr_partition.add(names.pop())
        # iterate over list of remaining names and find the D-closest
        # that is, the name with the least D-distance to the current clique
        # where the D-distance is the max pairwise distance to all clique members
        
        while len(curr_partition) < np.sqrt(N):
            min_dist, min_name = float("inf"), None
            for curr_name in names:
                curr_dist = clique_distance(curr_partition, curr_name, D, q)
                if curr_dist == float("inf"):
                    # curr_name is too far away. don't bother trying to update
                    continue
                min_dist,min_name = min((min_dist,min_name), (curr_dist,curr_name))
            # at this point, min_name tells me the D-closest node that is at most q away
            if not min_name:
                # no other node is closer than q away, so just break and create another partition
                break
            else:
                curr_partition.add(min_name)
                names.remove(min_name)
            
        partitions_list.append(curr_partition)
    
    # CHECKS FOR CORRECTNESS
    # for i in range(len(partitions_list)):
    #     for j in range(i):
    #         assert not partitions_list[i].intersection(partitions_list[j])

    return partitions_list

def get_bipartition(taxa_set, leaf):
    """
        find the bipartition of TAXA induced by LEAF 
        input:
            taxa_set - a set of labels to find the bipartition of
            leaf - a dendropy node for the leaf that induces the bipartition,
                    also a pointer into the tree
        return: 
            two sets of names that bipartitions TAXA
    
        implmentation:
            do DFS post_order (recursive since it's easier to think about)
    """
    def dfs_find_bipartition(u,v, taxa_set):
        """
            Given directed edge (u,v), return the LCA of TAXA_SET using DFS postorder

            return:
                if v is the LCA, return the bipartition - the left and right sets
                if v is ABOVE the LCA, return the bipartition
                else, return the set of all TAXA at or below V
        """
        # base case: v is leaf
        if v.is_leaf():
            if v.taxon.label in taxa_set:
                return {v.taxon.label}
            return set()

        # recursive case: v is not leaf
        child_sets = []
        for adj in v.adjacent_nodes():
            if adj is u:
                continue
            ret = dfs_find_bipartition(v, adj, taxa_set)
            if type(ret) is tuple:
                assert len(ret) == 2, "BIPARTITION ERROR"
                return ret
            assert type(ret) is set, "BIPARTITION ERROR"
            child_sets.append(ret)

        assert len(child_sets) == 2, "BIPARTITION ERROR"
        size1 = len(child_sets[0])
        size2 = len(child_sets[1])
        # adding is OK since we should never have duplicates
        if size1 + size2 == len(taxa_set):
            return tuple(child_sets)

        smaller = child_sets[size1 > size2]
        larger = child_sets[size1 <= size2]

        larger.update(smaller)
        return larger
    
    adj = leaf.adjacent_nodes()[0]
    assert not adj.is_leaf()
    return dfs_find_bipartition(leaf, adj,taxa_set)

def get_start_and_end(ptree, U, V):
    """
        find the two LCA of U and V in the current PTREE
        input:
            ptree - a PTree (the tree currently being constructed)
            U - a set of labels
            V - a set of labels
        return:
            (start_edge, end_edge) where start_edge is directed left-right 
            and specfies the direction of end_edge

        implementation:
            starting from an arbitrary root find the LCA of U
            if U is the arbitrary root, we have no information, 
                try again with V (guaranteed to have non root LCA)
            from this LCA as the new root, find the LCA of the other
            these are the two LCA's to return
    """
    def find_LCA(u,v, taxa_set):
        """
            Given directed edge (u,v) find LCA of TAXA_SET

            inputs:
                u,v - PNodes describing an edge in the tree
                taxa_set - a set of string names to find the LCA of

                return:
                    if v is the LCA, the directed edge (v,u)
                    if v is ABOVE the LCA, return  (v,u)
                    else, return the set of all TAXA in TAXA_SET at or below V

        """
        sets_so_far = set()

        if v.is_leaf() and v.get_val() in taxa_set:
            if len(taxa_set) == 1:
                return (v,u)
            sets_so_far.add(v)

        for adj in v.get_neighbors():
            if adj is u:
                continue
            ret = find_LCA(v, adj, taxa_set)
            if type(ret) == tuple:
                assert len(ret) == 2, "LCA ERROR"
                return ret
            assert type(ret) == set, "LCA ERROR"
            bigger = max(sets_so_far, ret, key=len)
            other = [sets_so_far, ret][bigger == sets_so_far]
            bigger.update(other)
            sets_so_far = bigger
            assert len(sets_so_far) <= len(taxa_set), "LCA ERROR"
            if len(sets_so_far) == len(taxa_set):
                return (v,u)
        assert len(sets_so_far) < len(taxa_set), "LCA ERROR"
        return sets_so_far

    # find the LCA of U. This may not return an edge towards the LCA of V
    # search from this 'LCA' to find the LCA. This returns an edge towards LCA of U (since we started search there)
    # search one more time for the LCA of U, but from the LCA of V. This should be a good LCA as well
    root = ptree.get_root().get_neighbors()[0]
    pre_lca = find_LCA(None, root, U)
    assert isinstance(pre_lca, tuple), "LCA ERROR"
    v_lca = find_LCA(None, pre_lca[0], V)
    assert isinstance(v_lca, tuple), "LCA ERROR"
    assert v_lca[1] != None, "LCA ERROR"
    u_lca = find_LCA(v_lca[0], v_lca[1], U)
    assert isinstance(v_lca, tuple), "LCA ERROR"
    assert v_lca[1] != None, "LCA ERROR"

    # CHECKS FOR CORRECTNESS
    # assert find_LCA(v_lca[1], v_lca[0], V) == v_lca
    # assert find_LCA(v_lca[0], v_lca[1], U) == u_lca
    
    # assert find_LCA(u_lca[1], u_lca[0], U) == u_lca
    # assert find_LCA(u_lca[0], u_lca[1], V) == v_lca
    return u_lca, (v_lca[1],v_lca[0])

def constraint_inc(names, matrix, query_func, 
                    r1=True, r2=True,r3=True, 
                    q_0=None, mult=8.0, ordering=None, 
                    constraints=[], warnings=False,
                    **kwargs):
    """
        The single, all-pupose function that implementations call with different arguments

        Estimate a tree represented by MATRIX, by iteratively inserting leaves into it,
        subject to the constraint trees in constraints
        inputs:
            names - an ordered list of names of taxa
            matrix - matrix of pairwise distances, indexed by the indices of NAMES (same ordering)
            constraints - a list of newick string trees, whose leaves are DISJOINT partition
                            of the taxa in NAMES. this implementation will accept constraint
                            lists that do not fully cover NAMES as long as they are DISJOINT
    """
    # convert matrix to "function" interface
    # true_tree = kwargs.get("true_tree", None) # newick tree
    # wrong_list = kwargs.get("wrong_list", None) # list of (leaf_insert,annotated_ptree) 
                                                # annotated_ptree has nodes that are annotated as "good" or "bad"
                                                # bad edges (ones that dont exist in the true tree) cause both (u,v) labelled to be "bad"
    # plan, at each step can manually recover which edges are not correct.
    # because of how ete3 works, annotate the nodes instead
    # ete3: get node with name "C" in tree t: C= t&"C"


    D = mat_to_func(names, matrix)
    if not q_0 or not ordering:
        s_node, q_0, ordering = prims_mst(names,D)
    
    # implement constraint_inc
    #   for each item we insert
    #       identify which constraint tree it exists in (if any)
    #       if it does, for each tree, keep a running list of wich
    #       taxa we've inserted up to this point for that tree
    #       before insertion, ID the LCA for all those taxa with insert and root (BFS)
    #       determine the bipartition induced by this LCA on the taxa
    #       determine which path induces the same biparition in the growin tree
    if constraints:
        names2constraint_inds = dict()
        names2constraint_leaf = dict()
        for n in names:
            names2constraint_inds[n] = None
        for i in range(len(constraints)):
            dendro_constraint = dendropy.Tree.get(data=constraints[i],schema='newick')
            for leaf in dendro_constraint.leaf_node_iter():
                names2constraint_inds[leaf.taxon.label] = i
                # map names to their leaves so we don't repeatedly search for them later
                names2constraint_leaf[leaf.taxon.label] = leaf

        added_so_far = [set() for _ in range(len(constraints))]

    # intialize tree
    ptree = PTree(search_fn=voting_with(query_func, mult*q_0, r1=r1,r2=r2,r3=r3, warnings=warnings))
    for name in ordering:
        before_p_str = ptree.make_newick_string()
        # if unconstrained, insert arbitrarily
        if not constraints:
            ptree.insert_leaf(PNode(name), D, start_leaf=ptree.get_root())
            continue
        # if no constraint tree for this name, insert arbitrarily        
        if name not in names2constraint_inds:
            # no constraint associated with this name. insert arbitrarily
            ptree.insert_leaf(PNode(name), D, start_leaf=ptree.get_root())
            continue
        
        idx = names2constraint_inds[name]
        added = added_so_far[idx]
        # if only 0,1,2 nodes inserted, any insertion is fine
        if len(added) < 3:
            ptree.insert_leaf(PNode(name), D, start_leaf=ptree.get_root())
            added.add(name)
            continue

        # must determine the bipartition that adding name causes in the constraint tree
        U,V = get_bipartition(added, names2constraint_leaf[name])
        # find the two end points of this bipartition in the current ptree
        # make sure one of them is a directed edge specifiying the direction of the other
        start_edge, end_edge = get_start_and_end(ptree, U,V)
        ptree.insert_leaf(PNode(name), D, start_edge=start_edge, end_edge=end_edge)
        added.add(name)

        # CHECKS FOR CORRECTNESS (these are redundant, for debugging)
        # assert satisfies_constraint(ptree.make_newick_string(), constraints[names2constraint_inds[name]])
        # for i,constraint in enumerate(constraints):
        #     assert satisfies_constraint(ptree.make_newick_string(), constraint)
    
    # CHECKS FOR CORRECTNESS
    for constraint in constraints:
        assert satisfies_constraint(ptree.make_newick_string(), constraint)
    
    return ptree.make_newick_string()

def nj_constrained_inc_with(query_func, **kwargs):
    def NJ_CONSTRAINED_INC(names, matrix, warnings=False):
        names2inds = dict()
        for i in range(len(names)):
            names2inds[names[i]] = i

        # convert matrix to "function" interface
        D = mat_to_func(names, matrix)
        s_node, q_0, ordered_names = prims_mst(names,D)

        p_list = partition(names, D, 8*q_0)

        # construct constraint trees
        constraints = []
        for clique in p_list:
            # convert clique (set) to list. ordering needed for NJ
            clique = list(clique)
            if len(clique) < 4:
                # doing NJ on this partition won't do anything here
                continue
            indices = [names2inds[name]  for name in clique]
            submatrix = matrix[indices,:][:,indices]
            constraint_tree = NJ(clique, submatrix)
            constraints.append(constraint_tree)
        return constraint_inc(names, matrix, query_func, q_0=q_0, ordering=ordered_names, constraints=constraints,warnings=warnings, **kwargs)
    return NJ_CONSTRAINED_INC



# INC_R1 is INC with only round 1 of voting: valid queries get 1 vote, invalid queries get no vote
# ties are broken by picking uniformly at random.
def INC_R1(names,matrix):
    return constraint_inc(names, matrix, queries.query, r1=True, r2=False, r3=False)

# INC_R1 but with NJ constraints
def INC_NJ_R1(names, matrix, warnings=False):
    return nj_constrained_inc_with(queries.query, r1=True, r2=False, r3=False)(names, matrix, warnings)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-phylip', help='the full path of the phylip distance matrix')
    parser.add_argument('--out', help='desired output file name of the NEWICK txt file')
    parser.add_argument('--inc', default=False,action="store_true", help='perform INC')
    parser.add_argument('--inc-nj', default=False,action="store_true", help='perform INC_NJ')
    args = parser.parse_args()

    with open(args.out, 'w') as out_file:
        print("Reading file...")
        names, matrix = parse_distance_matrix(args.in_phylip)

        if args.inc:
            print("Building INC Tree for tree...")
            newick_tree = INC_R1(names, matrix)
            out_file.write("INC\n")
            out_file.write(newick_tree + "\n")

        if args.inc_nj:
            print("Building INC_NJ Tree for tree...")
            newick_tree = INC_NJ_R1(names, matrix)
            out_file.write("INC_NJ\n")
            out_file.write(newick_tree + "\n")


####################################################
### Below are different version of INC that were ###
### explored, but not used for experiments.      ###
####################################################

# INC, but with each node query, we take the nearest (at most) 5 leaves in each direction
# and use those (at most) 5 leaves to summarize the subtree in that direction,
# to compare then with the leave to be inserted.
def INC_KNJ(names, matrix):
    return constraint_inc(names, matrix, queries.nj_query_with(5))

#INC_KNJ but using constraint trees as well
def INC_NJ_KNJ(names, matrix, warnings=False):
    return nj_constrained_inc_with(queries.nj_query_with(5))(names, matrix, warnings)


# INC with several rounds of different tie-breaking
def INC(names, matrix, warnings=False):
    """
        - Estimate a tree, represented by MATRIX by iteratively inserting leaves into it    
    """
    return constraint_inc(names, matrix, queries.query, warnings=warnings)

# INC with several tie-breaking rounds, and with constraints.
def INC_NJ(names, matrix, warnings=False):
    return nj_constrained_inc_with(queries.query)(names, matrix, warnings)

# INC where votes are made randomly for each internal node
def INC_RAND(names,matrix):
    return constraint_inc(names, matrix, queries.query, r1=False, r2=False, r3=False)

# different versions of INC with different rounds
def INC_R2(names,matrix):
    return constraint_inc(names, matrix, queries.query, r1=False, r2=True, r3=False)

def INC_R3(names,matrix):
    return constraint_inc(names, matrix, queries.query, r1=False, r2=False, r3=True)

def INC_R2R3(names,matrix):
    return constraint_inc(names, matrix, queries.query, r1=False, r2=True, r3=True)