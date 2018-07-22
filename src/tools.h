// File in HMMDecompositionDecision, created by Thien Le in July 2018

#ifndef TOOLS_H
#define TOOLS_H

#include "options.h"

const static int CMD_BUFFER_SIZE        = (int) 1e6;

extern int hmmbuild_job(hmmbuild_option_t * hmm_options);
extern int fasttree_job(fasttree_options_t * fasttree_options);
extern int hmmsearch_job(hmmsearch_options_t * hmm_options);

#endif