// File in inc_ml, created by Thien Le in July 2018

#ifndef TREE_H
#define TREE_H

#include "c_inc.h"

// Initializations
extern int parse_tree(INC_GRP * meta, MAP_GRP * map, ml_options * options);
extern int init_growing_tree(INC_GRP * meta, MAP_GRP * map, MST_GRP * mst);
extern int write_newick(BT * tree, char * filename, char ** name_map);
extern int get_degree(BT * tree, int idx);

extern int get_adj(BT * tree, int idx, int order);
extern void set_edge_master_idx(BT * tree, int ini, int dest, int val);
extern int get_edge_master_idx(BT * tree, int ini, int dest);
extern int get_edge_sample(BT * tree, int idx, int order);

BT * read_newick(MAP_GRP * map, char * filename, int tree_idx);

// Modifiers
extern int attach_leaf_to_edge(
  INC_GRP * meta,  
  MAP_GRP * map,
  MST_GRP * mst, 
  VOTE_GRP * vote, 
  int i
);

static const char N_GTREE[] = "gtree is null\n"; 

static const char F_MK_UNW_MAT_IN_PARSE_TREE[]  
  = "make_unweighted_matrices failed in parse_tree\n";

static const char LESS_THAN_THREE_TAX[]
  = "there are less than 3 taxa\n";

static const char F_MK_ADJ_IN_INIT_GTREE[]
  = "make adj failed in init gtree\n";

static const char F_ATTACH_LEAF[]
  = "attch leaf failed\n";

static const char F_SWAP_ADJ_IN_ATTACH_EDGE[]
  = "swap_adj failed in attach edge\n";

static const char F_MK_ADJ_IN_ATTACH_EDGE[]
  = "make adj failed in attach edge\n";

static const char F_CHECK_NODE_IN_SWAP_ADJ[]
  = "check node failed in swap adjacent\n";

static const char F_MK_ADJ_IN_MK_PAR[]
  = "make adjacent failed in make parent\n";

static const char F_MK_PAR_IN_RD_NW[]
  = "make parent failed in read_newick\n";

static const char F_SV_NAME_IN_RD_NW[]
  = "save name failed in read_newick\n"; 
#endif