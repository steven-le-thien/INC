// File in inc_ml, created by Thien Le in July 2018

#ifndef QUARTET_H
#define QUARTET_H
#include "c_inc.h"


static const int QUAD = 4;


extern int four_point_method(float ** d, int * u, int * res);

extern int new_quartets_raxml(
    char ** name_map, 
    int * u, 
    int * res, 
    ml_options * master_ml_options);

extern int ml_quartet(
    char ** name_map, 
    int * u, 
    int * res, 
    ml_options * master_ml_options, 
    double * M);

static const char F_MATCH_PRINT_IN_NEW_Q_RXML[] 
  = "matching print failed in new quartet raxml\n";

static const char F_MK_Q_IN_NEW_Q_RXML[]
  = "make quartet failed in new quartet raxml\n";

static const char F_PROCESS_Q_IN_NEW_Q_RXML[]
  = "process quartet failed in new quartet raxml\n";

static const char F_MATCH_PRINT_IN_ML_Q[]
  = "matching print failed in ml quartet\n";

static const char F_DO_Q_LL_IN_ML_Q[]
  = "do quartet ll failed in ml quartet\n";

static const char W_TMP_Q_EXIST[]
  = "write tmp q exists\n";

static const char F_RXML_JOB_IN_MAKE_Q[]
  = "raxml job failed in make quartet\n";

static const char F_RM_LBL_JOB_IN_MAKE_Q[]
  = "remove label job failed in make quartet\n";

static const char F_GET_LL_IN_DO_Q_LL[]
  = "get ll failed in do quartet ll\n";

#endif