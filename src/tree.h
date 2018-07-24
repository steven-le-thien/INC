// File in HMMDecompositionDecision, created by Thien Le in July 2018

#ifndef TREE_H
#define TREE_H

#include "c_inc.h"

// Initializations
extern int parse_tree(INC_GRP * meta, MAP_GRP * map, option_t * options);
extern int init_growing_tree(BT * tree, int * ordering, int * in_building, int num_sequence, int adj_in_mst);
extern int init_in_building(int * in_building, int n);

// Modifiers
extern int attach_leaf_to_edge(BT * growing_tree, int x, int addition_edge_parent, int additional_edge_child, int adj_in_mst);

#endif