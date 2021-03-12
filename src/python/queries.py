import queue
import numpy as np
from tree import PNode

def query(x, node, D):
    """
        A generic query function. Given a leaf X to insert while querying internal node NODE,
            which has three "leaf samples" for each direction,
        return the neighbor node X should belong in.

        See Tree.py for the definition of the PNode class
    """
    [(n1,l1),(n2,l2),(n3,l3)] = node.leaf_samples.items()
    val_x = x.get_val()
    val_1 = l1.get_val()
    val_2 = l2.get_val()
    val_3 = l3.get_val()
    s1 = D(val_x,val_1) + D(val_2, val_3)
    s2 = D(val_x,val_2) + D(val_1, val_3)
    s3 = D(val_x,val_3) + D(val_1, val_2)
    best_v, best_n = min([(s1,n1),(s2,n2),(s3,n3)], key=lambda x:x[0])
    return best_n



def k_closest(parent,child, k):
    """
        Return a collection of at most k PNodes that are leaves strictly below node parent in direction of child

        Return the k closest, found through BFS
    """
    k_closest = set()

    fringe = queue.Queue()
    fringe.put((parent,child))
    while not fringe.empty() and len(k_closest) < k:
        parent, child = fringe.get()
        if child.is_leaf():
            k_closest.add(child)
            continue

        for neighbor in child.get_neighbors():
            if neighbor is parent:
                continue
            fringe.put((child, neighbor))
    return k_closest


def min_q_ij(dm, leaf_list, node_to_ind):
    min_inds = (None,None)
    min_q = float("inf")
    N = len(leaf_list)
    # pre-calculate total dists
    total_dists = dict()
    for i in range(N):
        li = leaf_list[i]
        indi = node_to_ind[li]
        total_dists[indi] = 0
        for j in range(N):
            lj = leaf_list[j]
            indj = node_to_ind[lj]
            total_dists[indi] += dm[indi,indj]

    # calculate q values
    for i in range(N):
        for j in range(i):
            li, lj = leaf_list[i], leaf_list[j]
            indi, indj = node_to_ind[li],node_to_ind[lj]
            q_ij = (N - 2) * dm[indi,indj] - total_dists[indi] - total_dists[indj]

            (min_q, min_inds) = min((min_q, min_inds), (q_ij, (indi,indj)))

    return min_inds


def nj_query_with(num_leaves):
    def nj_query(x, node, D):
        """
            Query NODE from a binary tree to determine which of the three subtrees that X should belong in

            Use Neighbor joining to summarize the K closest leaves in each of the three directions (or fewer than K leaves).

            Implementation:
                We have three sets of k leaves. For each of the three sets, use NJ to merge them into a single "summary leaf"
                In each set, use the Q criterion locally within that set to merge two leaves into a cherry
                    After a merge, update the distance of this cherry from every other leaf (both in and out of the set)
        """
        assert not node.is_leaf()

        # D(label_1, label_2) is a functional interface and doesn't allow for alteration of the matrix

        neighbors = node.get_neighbors() # list of size 3
        leaves = [k_closest(node, neighbor, num_leaves) for neighbor in neighbors]# list of 3 sets, each of k leaves
        

        ## Set up the distance matrix
        node2ind = dict()
        ind2node = dict()

        node2ind[x] = 0
        ind2node[0] = x
        i = 1
        for l_set in leaves:
            for leaf in l_set:
                node2ind[leaf] = i
                ind2node[i] = leaf
                i += 1
        total_leaves = len(ind2node)
        dm = np.zeros((total_leaves,total_leaves))
        for i in range(total_leaves):
            for j in range(i):
                dm[i,j] = D(ind2node[i].get_val(), ind2node[j].get_val())
                dm[j,i] = dm[i,j]

        ## for each of the three leaf lists, reduce them to a single node, updating distances
        for l_set in leaves:
            while len(l_set) > 1:
                # find indices to reduce
                i,j = min_q_ij(dm, list(l_set), node2ind)
                # join i and j
                node_i, node_j = ind2node[i], ind2node[j]
                merged_node = PNode("({0},{1})".format(node_i.get_val(), node_j.get_val()))
                l_set.remove(node_i)
                l_set.remove(node_j)
                l_set.add(merged_node)
                # delete nodes
                del node2ind[node_j]
                del node2ind[node_i]
                # delete index j from use, but reuse i for the merged node
                del ind2node[j]
                ind2node[i] = merged_node
                node2ind[merged_node] = i
                # 1 calculate distances of new node to all others
                old_d_ij = dm[i,j]
                for k in ind2node.keys():
                    dm[i,k] = 0.5 * (dm[i,k] + dm[j,k] - old_d_ij)
                    dm[k,i] = dm[i,k]

        for l_set in leaves:
            assert len(l_set) == 1
        assert len(node2ind) == 4
        assert len(ind2node) == 4

        merged_nodes = [l_set.pop() for l_set in leaves]
        [ind0, ind1, ind2] = [node2ind[merged_n] for merged_n in merged_nodes]
        indx = node2ind[x]
        assert indx == 0


        v0 = dm[indx, ind0] + dm[ind1, ind2]
        v1 = dm[indx, ind1] + dm[ind0, ind2]
        v2 = dm[indx, ind2] + dm[ind0, ind1]
        n0, n1, n2 = neighbors
        best_v, _, best_n = min((v0, 0, n0), (v1, 1, n1), (v2, 2, n2))
        return best_n

    return nj_query