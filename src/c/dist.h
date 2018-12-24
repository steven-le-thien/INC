// File in inc_ml, created by Thien Le in July 2018

#ifndef DIST_H
#define DIST_H

#include "c_inc.h"

extern int parse_distance_matrix(INC_GRP * meta, MAP_GRP * map, option_t * options);
extern double compute_logdet_distance (char ** data, int numSites);
extern double compute_jc_distance (char ** data, int numSites);

#endif