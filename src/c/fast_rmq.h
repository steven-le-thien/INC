#ifndef FAST_RMQ_H
#define FAST_RMQ_H

#include "c_inc.h"

extern int fast_rmq_init(int n, int * a,  RMQ_T * RMQ);
extern int fast_rmq(int qs, int qe, int * min_idx, RMQ_T * RMQ);

static const char F_NUM_TO_CFG_IN_FAST_RMQ_INIT[]
  = "num_to_config failed in fast_rmq_init\n";

//static const char 

#endif