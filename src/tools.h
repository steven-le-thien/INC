// File in inc_ml, created by Thien Le in July 2018

#ifndef TOOLS_H
#define TOOLS_H

#include "options.h"

const static int CMD_BUFFER_SIZE        = (int) 1e6;

extern int fastphylo_job(fp_options * options);
extern int constraint_inc(int argc, option_t * options);
extern int fasttree_job(option_t * options);
extern int upp_job(option_t * options);
extern int subset_job(option_t * options);
extern int nw_utils_job(option_t * options);
extern int rm_label_job(option_t * options);
extern int distance_matrix_job(option_t * options);
#endif