// File in inc_ml, created by Thien Le in July 2018

#ifndef QUARTET_H
#define QUARTET_H
#include "c_inc.h"


extern int four_point_method(float ** d, int * u, int * res);
extern int new_quartets_raxml(char ** name_map, int * u, int * res, ml_options * master_ml_options);
extern int ml_quartet(char ** name_map, int * u, int * res, ml_options * master_ml_options, double * M);
#endif