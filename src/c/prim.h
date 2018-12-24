// File in inc_ml, created by Thien Le in July 2018

#ifndef PRIM_H
#define PRIM_H

#include "c_inc.h"

extern int prim(INC_GRP * meta, MST_GRP * mst);
extern int prim_on_small_graph(int n, GRAPH * graph, MST_GRP * mst, char * distance_model, char ** data);

#endif