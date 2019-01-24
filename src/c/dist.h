// File in inc_ml, created by Thien Le in July 2018

#ifndef DIST_H
#define DIST_H

#include "c_inc.h"

extern int parse_distance_matrix(INC_GRP * meta, MAP_GRP * map, ml_options * options);
extern double compute_logdet_distance (char ** data, int numSites);
extern double compute_jc_distance (char ** data, int numSites);
extern double dist_from_msa(msa_t * msa, DIST_MOD distance_model, int i, int j, double correction);

#endif