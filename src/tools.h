// File in HMMDecompositionDecision, created by Thien Le in July 2018

#ifndef TOOLS_H
#define TOOLS_H

#include "options.h"

const static int CMD_BUFFER_SIZE        = (int) 1e6;

extern int fastphylo_job(fp_options * options);
extern int constraint_inc(int argc);
extern int fasttree_job(option_t * options);
extern int upp_job(option_t * options);

#endif