// File in inc_ml, created by Thien Le in July 2018

#ifndef QUARTET_H
#define QUARTET_H
#include "c_inc.h"


extern int four_point_method(float ** d, int u1, int u2, int up, int x, int * res);
extern int new_quartets_raxml(char * up, char * u1, char * u2, char * x, int * res, ml_options * master_ml_options);
extern int ml_quartet(char * up, char * u1, char * u2, char * x, int * res, ml_options * master_ml_options, float * M,  int * revote_power);
#endif