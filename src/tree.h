// File in HMMDecompositionDecision, created by Thien Le in July 2018

#ifndef TREE_H
#define TREE_H
typedef struct tree{
	int n;  			// number of nodes in the tree
	int root; 			// root of the tree, -1 if the tree is supposed to be unrooted
	int * adj_list[3]; 	// binary tree assumes that degree of nodes is no larger than 3
	int * degree; 		// degree of a node in the tree
	int * index_in_master_name_map;
} BT;

extern int parse_tree(BT ** constrained_trees, option_t * options, int * constraint_to_master_map, int master_to_ctreeindex, int num_sequence, char ** master_to_ctree_map);
extern void tree_destructor(BT * tree); 		// null destructor
extern int is_leave(int cur_node, BT * tree);

extern int init_growing_tree(BT * tree, int * ordering, int * in_building);
extern int init_in_building(int * in_building, int n);

#endif