import numpy as np
import random
# for each bit, decide whether to flip a coin with probability p
# if we flip a coin, the bit is changed with probability 1/2, same otherwise.
# this is equivalent to...
# for each bit, flip it with probability p/2
def noisy_mutation(bit_string, branch_len, scale):
    """Given a bit_string, for each bit, flip coin with prob p, do nothing prob 1-p
        If we flip a coin, change with prob a 1/2, stay same with prob 1/2
        where p = 1 - exp(-branch_len/scale)
        
        But this is equivalent to
        For each bit, flip it with probability 0.5(1-exp(-branch_len)) 
        """
    # 1 if we flip, 0 if we decide not to flip
    do_flip = np.random.binomial(1,(1-np.exp(-branch_len/scale))*0.5, bit_string.shape)
    return np.bitwise_xor(bit_string, do_flip)

def frac_hamming_dist(s1,s2):
    """Find the fractional hamming distance between bitstrings s1 and s2"""
    return np.sum(np.abs(s1-s2))/s1.size

def noisy_dist(s1,s2,scale):
    """ Given the bit sequences corresponding to s1, s2, give the estimated distance
        taxa
        
        If the change probability is p = 0.5(1-exp(-branch_len/scale)), and p is also the 
        fractional hamming distance, then 
        branch_len = -scale*ln(1-2p)
        """
    ham_dist = frac_hamming_dist(s1,s2)
    return -scale * np.log(1-2*ham_dist), ham_dist

def measure_epsilon(orig, new, smallest_branch_len):
    return np.amax(np.abs(orig-new))/smallest_branch_len

def create_noisy_dist_matrix(dendro_tree,num_bits,scale=None,ret_ham_dist=False):
    """ Given a dendropy tree, sample a noisy distance matrix from it using 
        DNA-based evolution model. 
        
        DENDRO_TREE - the dendropy tree to be mutated
        NUM_BITS - specifies the number of bits in the sequence to be mutated
                   make this smaller for more noise
        SCALE    - specifies the multiplier that each edge will be divided by
                   before mutating. multiply SCALE * log(1-2h) to get a distance
                   that in expectation is the true distance
                   
                   Rule of thumb - tune this until the maximum hamming distance is ~0.4-0.45
        MAX_HAM_DIST - specifies whether to return the maximum hamming distance between leaves in the tree
        """
    # assign arbitrary node an arbitrary bit string
    # do DFS over adjacent nodes from the first node
    # over each edge, do a noisy_mutation according to the edge length
    # once all leaves have a bit string, 
    #   construct the final matrix using noisy_dist
    
    if not scale:
        # by default, use the largest inter-leaf distance
        # this can change to some other behavior if needed
        scale = float("-inf")
        tpdm = dendro_tree.phylogenetic_distance_matrix()
        for tax1, tax2 in tpdm.distinct_taxon_pair_iter():
            scale = max(scale, tpdm.distance(tax1,tax2))
        
    # first sequence
    sequence = np.random.binomial(1,0.5,int(num_bits))
    for edge in dendro_tree.preorder_edge_iter():
        if edge.tail_node == None:
            # The very first edge has None for the tail
            # assign the head_node to be initial sequence
            edge.head_node.bit_sequence = sequence
        else:
            # otherwise, mutate the bit sequence in TAIL using edge length as the parameter
            # assign this new bit sequence to the HEAD
            edge.head_node.bit_sequence = noisy_mutation(edge.tail_node.bit_sequence, edge.length, scale)

    # all leaves should now have a bit sequence associated with them. construct the 
    # noisy distance matrix based on these sequences
    leaves = dendro_tree.leaf_nodes()
    num_leaves = len(leaves)
    d_matrix = np.zeros((num_leaves,num_leaves))
    max_ham_dist = float("-inf")
    for i in range(num_leaves):
        leafi = leaves[i]
        li = int(leafi.taxon.label) - 1
        for j in range(i):
            leafj = leaves[j]
            lj = int(leafj.taxon.label) - 1
            d_matrix[li,lj], ham_dist = noisy_dist(leafi.bit_sequence,leafj.bit_sequence, scale)
            d_matrix[lj,li] = d_matrix[li,lj]
            max_ham_dist = max(max_ham_dist, ham_dist)
    if ret_ham_dist:
        return d_matrix, max_ham_dist
    return d_matrix