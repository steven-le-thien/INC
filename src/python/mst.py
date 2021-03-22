from tree import PNode
def prims_mst(names, D, ret_attach=False):
    """
        names   - names of leaves
        D       - function mapping pairs of names to their interleaf distances 
                    (describes the edges in the FC graph)
        ret_attach - if true,return an extra list describing the order of 'attaching nodes'.
            For the ith node in the prim-ordering, the 'attaching node' is that 
            node that the ith node was connected to upon being added to the MST
            This is None for the 0th node.
        
        return 
            1. a pointer to a PNode in the MST of this graph (unnecessary?)
            2. q_0, the largest weight in this MST
            3. list of the names, in order of they were inserted into the MST
            4. (optional) list of names for 'attaching nodes' for local search
        Use 1-indexing for the nodes 
        
        """
    # Implementation details
    # Since the graph is fully connected, don't bother with priority queues
    # and instead just do a naive "update everyone's distance" version of Prims
    # that tracks the closest distance so far
    # The rt remains O(E) = O(n^2) and avoids unnecessary PQ operations
    
    # Adapted from an MST agnostic version of this algorithm 
    # from before I read the paper
    
    n = len(names)
    q_0 = float("-inf")
    # map each vertex V still on the fringe, to its distance from the growing MST
    # and the corresponding node in the MST that V is closest to
    leaf_dists = dict((PNode(name),(float("inf"), None)) for name in names)
    closest_leaf = leaf_dists.popitem()[0]
    leaves_prim_order = [closest_leaf]
    if ret_attach:
        attaching_nodes = []
    while leaf_dists:
        closest_dist = float("inf")
        next_closest = None
        # update leaf_dists dict, tracking the nearest leaf
        for l in leaf_dists:
            leaf_dists[l] = min(
                            leaf_dists[l], 
                            (D(l.get_val(), closest_leaf.get_val()), closest_leaf),
                                   key=lambda x:x[0])
            if leaf_dists[l][0] < closest_dist:
                closest_dist, adj_node = leaf_dists[l]
                next_closest = l
        # next_closest is the next node to insert into the MST. insert now
        closest_leaf = next_closest
        adj_node.add_neighbor(closest_leaf, closest_dist)
        # update q_0 with the length of this edge
        q_0 = max(q_0, closest_dist)
        # updates lists
        leaves_prim_order.append(closest_leaf)
        if ret_attach:
            attaching_nodes.append(adj_node)
        del leaf_dists[closest_leaf]
        
    # return the list of names of the nodes in the tree
    leaves = [l.get_val() for l in leaves_prim_order]
    if ret_attach:
        ret_att_nodes = [None] + [a.get_val() for a in attaching_nodes]
        
    if ret_attach:
        return closest_leaf, q_0, leaves, ret_att_nodes
    return closest_leaf, q_0, leaves