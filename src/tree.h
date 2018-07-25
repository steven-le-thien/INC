// File in HMMDecompositionDecision, created by Thien Le in July 2018

#ifndef TREE_H
#define TREE_H

#include "c_inc.h"

// Initializations
extern int parse_tree(INC_GRP * meta, MAP_GRP * map, option_t * options);
extern int init_growing_tree(INC_GRP * meta, MST_GRP * mst);

// Modifiers
extern int attach_leaf_to_edge(INC_GRP * meta, MST_GRP * mst, VOTE_GRP * vote, int i);
#endif