// File in inc_ml, created by Thien Le in July 2018

#ifndef PRIM_H
#define PRIM_H

#include "c_inc.h"

extern int prim(INC_GRP * meta, MST_GRP * mst);
extern int prim_on_small_graph(int n, GRAPH * graph, MST_GRP * mst, DIST_MOD distance_model, char ** data);
extern int parse_initial_tree_as_mst(INC_GRP * meta, MST_GRP * mst);

#endif