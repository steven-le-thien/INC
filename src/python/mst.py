from tree import PNode
def prims_mst(n, D, ret_attach=False):
    """
        n - number of leaves (all nodes in the fully connected graph)
        D - function mapping pairs of Pnodes to their interleaf distances 
            (describes the edges in the FC graph)
        ret_attach - if true,return an extra list describing the 'attaching nodes'
            for the ith node in the prim-ordering, the 'attaching node' is that 
            node that the ith node was connected to upon being added to the MST
            This is None for the 0th node.
        
        return 
            1. a pointer to a PNode in the MST of this graph (unnecessary?)
            2. q_0, the largest weight in this MST
            3. list of all node, in order they were inserted into the MST
            (4. list of 'attaching nodes' for local search )
        Use 1-indexing for the nodes 
        
        """
    # Implementation details
    # Since the graph is fully connected, don't bother with priority queues
    # and instead just do a naive "update everyone's distance" version of Prims
    # that tracks the closest distance so far
    # The rt remains O(E) = O(n^2) and avoids unnecessary PQ operations
    
    # Adapted from an MST agnostic version of this algorithm 
    # from before I read the paper
    
    q_0 = float("-inf")
    # map each vertex V still on the fringe, to its distance from the growing MST
    # and the corresponding node in the MST that V is closest to
    leaf_dists = dict((PNode(i),(float("inf"), None)) for i in range(1,n+1))
    closest_leaf = leaf_dists.popitem()[0]
    leaves_prim_order = [closest_leaf]
    if ret_attach:
        attaching_nodes = []
    while leaf_dists:
        closest_dist = float("inf")
        next_closest = None
        # update leaf_dists dict, tracking the nearest leaf
        for l in leaf_dists:
            leaf_dists[l] = min(leaf_dists[l], (D(l, closest_leaf), closest_leaf),
                                   key=lambda x:x[0])
            if leaf_dists[l][0] < closest_dist:
                closest_dist, adj_node = leaf_dists[l]
                next_closest = l
        # next_closest is the next node to insert into the MST
        closest_leaf = next_closest
        adj_node.add_neighbor(closest_leaf, closest_dist)
        # update q_0 with the length of this edge
        q_0 = max(q_0, closest_dist)
        # updates lists
        leaves_prim_order.append(closest_leaf)
        if ret_attach:
            attaching_nodes.append(leaf_dists[closest_leaf][1])
        del leaf_dists[closest_leaf]
        
    # make a copy of the list of leaves so that the user receives fresh, 
    # unconnected nodes
    new_nodes = [PNode(i) for i in range(1,n+1)]
    leaves = [new_nodes[l.get_val() - 1] for l in leaves_prim_order]
    if ret_attach:
        ret_att_nodes = [None] + [new_nodes[a.get_val()-1] for a in attaching_nodes]
        
    if ret_attach:
        return closest_leaf, q_0, leaves, ret_att_nodes
    return closest_leaf, q_0, leaves