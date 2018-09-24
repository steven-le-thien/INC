#ifndef FAST_RMQ_H
#define FAST_RMQ_H

#include "c_inc.h"

extern int fast_rmq_init(int n, int * a,  RMQ_T * RMQ);
extern int fast_rmq(int qs, int qe, int * min_idx, RMQ_T * RMQ);

#endif