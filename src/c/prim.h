// File in inc_ml, created by Thien Le in July 2018

#ifndef PRIM_H
#define PRIM_H

#include "c_inc.h"

extern int prim(INC_GRP * meta, MST_GRP * mst);
extern int prim_on_small_graph(int n, GRAPH * graph, MST_GRP * mst, DIST_MOD distance_model, char ** data);
extern int parse_initial_tree_as_mst(INC_GRP * meta, MST_GRP * mst);

static const char F_HEAD_POP_IN_PRIM[]
  = "head pop failed in prim\n";

static const char F_UPDATE_IN_PRIM[]
  = "update failed in prim\n";

static const char F_SWAP_NODE_IN_HEAPIFY[]
  = "swap node failed failed in heapify\n";

static const char F_HEAPIFY_IN_POP_HEAD[]
  = "heapify failed in pop head\n";

static const char F_SWAP_NODE_IN_UPDATE[]
  = "swap node failed in update\n";

static const char F_SWAP_NODE_IN_POP_HEAD[]
  = "swap node failed in pop head\n";

static const char F_REC_IN_HEAPIFY[]
  = "recursion failed in heapify\n";

#endif