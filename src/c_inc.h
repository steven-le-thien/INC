#ifndef DIST_H
#define DIST_H

typedef struct inc_grp{
	int 		n_taxa;
	int 		n_ctree;

	BT * 		g_tree;
	BT **		c_tree;		// n_ctree-length array of constraint trees
	int * 		ordering; 	// n_taxa-length array of insertion order of leaf taxon into the tree indexed through the master index
	int * 		visited;	// n_taxa-lenght array of visted taxa in the building tree 	
						
	float ** 	d; // n_taxa by n_taxa array of distance matrix
} INC_GRP;

typedef struct mapping{
	int 	n_taxa;

	int * 	master_to_ctree; 	// n_taxa length array mapping the master indexing scheme to the index of the binary constraint tree
	int *	master_to_cidx;  	// n_taxa length array mapping the master indexing scheme to the index within the binary constraint tree
	char **	master_to_name;		// n_taxa length array mapping the master indexing scheme to its name
} MAP_GRP;

typedef struct mst{
	int 	n_taxa;
	int * 	prim_ord;		// n_taxa length array denoting the Prim's insertion order in building the mst
	int * 	prim_par;		// n_taxa length array denoting the parent of each node in the mst (this can be viewed as the tree itself)
	float 	max_w;			// maximum edge weight in the mst
} MST_GRP;

// Special edges 
typedef struct s_edge{
	int 	c; 			// child node, usually the starting node
	int 	p;			// parent node, usually the avoiding node
} S_EDGE;

typedef struct vote{
	int 	n_taxa;	
	int * 	vote;			// at-most-4-n_taxa length array storing the vote for each edge in the growing tree. this is allocated at maximum (4 * num_sequence) but only a few entried can be valid

	S_EDGE 	valid_st; 		// the valid 'component', indexed in growing tree scheme
	S_EDGE 	c_lca; 			// lca of taxa in the induced constraint tree, this follows the constraint's indexing scheme
	S_EDGE 	st_lca;			// lca of taxa in the growing tree following the first bipartition, this follows the growing tree's indexing scheme
	S_EDGE 	nd_lca;			// lca of taxa in the growing tree following the second bipartition, this follow the growing tree's indexing scheme

	S_EDGE 	ins;			// insertion edge (edge with the most vote), this follws the growing tree's indexing scheme

} VOTE_GRP;

typedef struct ins{
	int 	c;				// insertion child
	int 	d
} INS_GRP;	

// Trees
typedef struct edge{ // this is an adjacency list edge, meaning that it is directed and that the source is know when the structure is accessed
	int 	dest; 			// destination 
	int 	master_idx;		// indexing into the vote array (only for the voting), it maybe possible to remove this field
	int 	sample; 		// leaf sample at the destination
} BT_edge;

typedef struct tree{
	int n_node;  					// number of nodes in the tree
	BT_edge * adj_list[3]; 			// n_node array of adjacent nodes
	int * degree; 					// n_node array of degee
	int * master_idx_map;			// n_node indexing from the constrained tree's scheme to the master scheme, this is not a surjection
} BT;

#endif