#ifndef FAST_MST_H
#define FAST_MST_H

#include "c_inc.h"

extern int fast_mst(char ** data, int n, DIST_MOD distance_model, int seed, MST_GRP * mst, int *** disjoint_subset);

#endif 