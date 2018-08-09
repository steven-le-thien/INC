from tree import query
from mst import prims_mst
import random

def voting_with(query_func, q):
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
        
    def voting_dfs(start_leaf, insert_leaf, D, q, weighted, edge_candidates=None):
        """ 
            start_leaf - a pointer to any arbitrary leaf in the tree
            D - a function mapping pairs of leaves to their D-distance
            q - if not None, use only valid queries (where diameter <= q)
                else use all queries
            weighted - if True, weight queries as 1/diam**2. else give only one vote
            
            return - a list of edges with maximal vote weight
        """
        # maps edges (u,v) to (relative) scores
        score_of = dict() 
        # use DFS search, starting from the edge subtending START_LEAF 
        score = 0
        # INVARIANT: all edges are ordered to be (n1,n2) 
        # where n1 is closer to START_LEAF than n2
        edge = (start_leaf, start_leaf.get_neighbors()[0])
        queue = [(score,edge)]
        
        exists_valid_query = True
        if q and not weighted:
            # in round 1, check for condition where no valid query exists
            exists_valid_query = False
        while queue:
            score, edge = queue.pop()
            score_of[edge] = score
            u,v = edge
            if v.is_leaf():
                continue
                
            # two children vertices of the deeper node v
            c1, c2 = v.neighbors_not(u)
            
            # determine the diameter
            l1,l2,l3 = v.leaf_samples.values()
            diameter = max(D(l1,l2),D(l1,l3),D(l1,insert_leaf),
                           D(l2,l3),D(l2,insert_leaf), D(l3,insert_leaf))
            
            # if q is None or nonpositive, any query is valid
            vote_wt = (not q) or (q <= 0) or (diameter <= q)
            # reweight if weighted, else keep as-is
            vote_wt *= diameter**(-2.0 * weighted)
            
            if not vote_wt:
                # don't query this node as we will ignore it.
                queue.append((score, (v, c1)))
                queue.append((score, (v, c2)))
                continue
            elif q and not weighted:
                # in round 1, check for condition where no valid query exists
                exists_valid_query = True
            
            best_neighbor = query_func(insert_leaf, v, D)
            if best_neighbor is u:
                queue.append((score - vote_wt, (v, c1)))
                queue.append((score - vote_wt, (v, c2)))
            elif best_neighbor is c1:
                queue.append((score + vote_wt, (v, c1)))
                queue.append((score, (v, c2)))
            elif best_neighbor is c2:
                queue.append((score, (v, c1)))
                queue.append((score + vote_wt, (v, c2)))
            else:
                raise RuntimeError("query_func did not return a neighbor internal_node")
        
        # find all candidate edges with the maximum score
        max_score = float("-inf")
        max_edges = set()
        for edge,score in score_of.items():
            if edge_candidates and edge not in edge_candidates:
                continue
            if score > max_score:
                max_edges = set([edge])
                max_score = score
            elif score == max_score:
                max_edges.add(edge)
            else: # do nothing
                pass
        
        if q and not exists_valid_query:
            import warnings
            msg = "No valid query exists.\nTree:\n{0}.\nLeaf:\n{1}\nMatrix:\n{2}\n"
            msg = msg.format(best_edge[0].make_newick_string(), insert_leaf.get_val(), D)
            print(msg)
            warnings.warn(msg)

        return max_edges
    
    def voting_search(start_leaf, insert_leaf, D):
        """ 
            start_leaf - a pointer to any arbitrary leaf in the tree
            insert_leaf - the leaf to be inserted
            D - a function mapping pairs of leaves to their D-distance
            
            return - the edge to break to insert insert_leaf, 
                        determined by 4 rounds of voting with tie-breaking
        """
        # round 1
        candidates = voting_dfs(start_leaf,insert_leaf,D, q, False)
        # round 2
        if len(candidates) > 1:
            candidates = voting_dfs(start_leaf,insert_leaf,D, q, True, candidates)
        # round 3
        if len(candidates) > 1:
            candidates = voting_dfs(start_leaf, insert_leaf,D, None, True, candidates)
        # round 4
        if len(candidates) > 1:
            candidates = [random.choice(list(candidates))]
        
        # one candidate left
        assert len(candidates) == 1, "VOTING FAILURE"
        winner = candidates.pop()
        return winner
    
    return voting_search

