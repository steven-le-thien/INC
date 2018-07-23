// File in HMMDecompositionDecision, created by Thien Le in July 2018

#ifndef TREE_H
#define TREE_H
typedef struct tree{
	int n;  			// number of leaves in the tree
	int root; 			// root of the tree, -1 if the tree is supposed to be unrooted
	int * adj_list[3]; 	// binary tree assumes that degree of nodes is no larger than 3
	int * degree; 		// degree of a node in the tree
	int * parent_map;	// identify the parent, -1 if the node is a leaf
	int * is_in_building_tree // boolean array to check if a current leaf is in the building tree
	char ** name_map;
} BT;

extern int parse_tree(BT ** constrained_trees, option_t * options);
extern void tree_destructor(BT * tree); 		// null destructor
extern BT * read_newick(char * filename); // newick constructor
extern int is_leave(int cur_node, BT * tree);

extern int init_growing_tree(BT * tree, int * ordering);

#endif