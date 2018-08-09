from tree import query
from mst import prims_mst

def local_search_with(query_func):
    """
        Using QUERY_FUNC to determine which direction to go locally,
        return local_search function that uses it
        
        interface:
            QUERY_FUNC(insert_leaf, internal_node, D)
    """
    def local_search(start_leaf, insert_leaf, D):
        """
            Starting from START_LEAF in the tree, do local hill climbing search
                to find the best node,
            use distance matrix D in sub calls to node queries
            
            return the edge to break in order to insert INSERT_LEAF
        """
        visited = set()
        # don't actually start at a leaf but at an adjacent internal node
        curr_node = start_leaf.get_neighbors()[0]
        curr_neighbors = curr_node.get_neighbors()
        neighbors_not_visited = set(curr_neighbors) - visited
        while len(neighbors_not_visited) > 1:
            #if curr_node.is_leaf():
            #    raise RuntimeError("Bad")
            best_neighbor = query_func(insert_leaf,curr_node, D)
            visited.update(set(curr_neighbors) - set([best_neighbor]))
            prev_node = curr_node
            curr_node = best_neighbor
            if curr_node.is_leaf():
                break
            curr_neighbors = curr_node.get_neighbors()
            neighbors_not_visited = set(curr_neighbors) - visited
        # return the edge between curr_node and prev_node
        return (curr_node, prev_node)
    return local_search

