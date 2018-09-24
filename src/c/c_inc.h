// File in inc_ml, created by Thien Le in July 2018

#ifndef C_INC_H
#define C_INC_H

typedef struct ml_options{
	char * input_alignment;		// currently only accepts path
	char * output_prefix;		// currently only accepts path
	char * init_tree_name;		// currently only accepts path 
	char * init_d_name;			// currently only accepts path
	char * guide_tree_name; 	// currently only accepts path

	// Use distance matrix 
	int use_distance_matrix;

	// Constraint trees settings
	int use_constraint;
	int recompute_constraint_trees;

	int use_subtree_for_constraint_trees;
	int use_raxml_for_constraint_trees;
	int use_fasttree_for_constraint_trees;

	// Quartet method
	int use_four_point_method_with_distance;
	int use_four_point_method_with_tree;        
	int use_new_quartet_raxml;
	int use_ml_method;

	// Distance setting
	// DISTANCE_MODEL                          "logDet"
	char * distance_model; //0 for no-distance mode1 for JC, 2 for log-det, ;

	// PASTA decomposition setting
	int ss_threshold; 
} ml_options;


// Trees
typedef struct m{ // this is an adjacency list edge, meaning that it is directed and that the source is know when the structure is accessed
	int 	dest; 			// destination 
	int 	master_idx;		// indexing into the vote array (only for the voting), it maybe possible to remove this field
	int 	sample; 		// leaf sample at the destination
	double 	weight;
} BT_edge;

typedef struct tree{
	int n_node;  					// number of nodes in the tree
	BT_edge ** adj_list; 			// n_node array of adjacent nodes
	int * degree; 					// n_node array of degee
	int * master_idx_map;			// n_node indexing from the constrained tree's scheme to the master scheme, this is not a surjection
} BT;

// Whole graphs 
typedef struct list{
	int 		dest;
	struct list * 	next;
} ADJ_LIST;

typedef struct graph{
	int n;
	ADJ_LIST ** last_node; // O(1) insertion
	ADJ_LIST ** adjacency_list;
} GRAPH;

// Groups
typedef struct inc_grp{
	int 		n_taxa;
	int 		n_ctree;

	BT * 		gtree;
	BT **		ctree;		// n_ctree-length array of constraint tree
	int * 		visited;	// n_taxa-lenght array of visted taxa in the building tree, this is in the master indexing scheme
						
	float ** 	d; 			// n_taxa by n_taxa array of distance matrix
	float ** 	dm; 		// n_taxa by n_taxa array of distance matrix on the FastTree

	ml_options * master_ml_options; 
} INC_GRP;

typedef struct mapping{
	int 	n_taxa;

	int *	master_to_midx;  	// n_taxa length array mapping the master indexing scheme to the index within the binary constraint tree

	int * 	master_to_ctree; 	// n_taxa length array mapping the master indexing scheme to the index of the binary constraint tree
	int *	master_to_cidx;  	// n_taxa length array mapping the master indexing scheme to the index within the binary constraint tree
	int * 	master_to_gidx;
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

	// WARNING: redo initialization if the struct changes
	S_EDGE 	valid_st; 		// the valid 'component', indexed in growing tree scheme
	S_EDGE 	c_lca; 			// lca of taxa in the induced constraint tree, this follows the constraint's indexing scheme
	S_EDGE 	st_lca;			// lca of taxa in the growing tree following the first bipartition, this follows the growing tree's indexing scheme
	S_EDGE 	nd_lca;			// lca of taxa in the growing tree following the second bipartition, this follow the growing tree's indexing scheme

	S_EDGE 	ins;			// insertion edge (edge with the most vote), this follws the growing tree's indexing scheme

	int 	ctree_idx; 
} VOTE_GRP;

typedef struct options{
    int num_options;
    int num_trees;

    int input_index;
    char * input_name;

    int output_index;
    char * output_name;

    int tree_index;
    char ** tree_names ;
} option_t;

typedef struct heap_node{
    int key;
    int value;
} heap_node;

typedef struct min_heap{
    int size;
    int capacity;
    int * pos;
    heap_node ** heap;
} min_heap;

typedef struct msa {
    int     N;          // size of one sequence
    int     num_seq;    // number of sequences
    char**  msa;
    char**  name;
} msa_t;

typedef struct rmq{
    // Statistics
    int         n;
    int         num_blk;
    int         num_exhaust;
    int         blk_sz;

    // Arrays and mappings
    int         * a;
    int         * blk_idx_to_config;

    // Tables
    int         ** sparse_table;        // num_blk x log2(num_blk) array, where A[i][j] is the index of the block with the smallest value from the jth block to the j + 2^i block
    int         *** exhaustive_table;   // config_num x blk_sz x blk_sz array where A[i][j][k] is the relative index of the smallest value from the jth position to the kth position for the ith configuration
} RMQ_T;

typedef struct lca_t{
	// LCA stuff
	int 	* euler_tour;
	int 	* level;
	int 	* first_occurence;
	int 	* visited;
	BT 		* tree;
	RMQ_T 	* RMQ;

	// Smallest distance stuff
	double 	* d_from_root;
} LCA_T;

extern int constraint_inc_main(int argc, char ** argv, ml_options * master_ml_options);

#endif