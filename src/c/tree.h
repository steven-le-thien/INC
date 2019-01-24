// File in inc_ml, created by Thien Le in July 2018

#ifndef TREE_H
#define TREE_H

#include "c_inc.h"

// Initializations
extern int parse_tree(INC_GRP * meta, MAP_GRP * map, ml_options * options);
extern int init_growing_tree(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst);
extern int write_newick(BT * tree, char * filename, char ** name_map);

BT * read_newick(MAP_GRP * map, char * filename, int tree_idx);

// Modifiers
extern int attach_leaf_to_edge(
    INC_GRP * meta,  
    MAP_GRP * map,
    MST_GRP * mst, 
    VOTE_GRP * vote, 
    int i
);
#endif