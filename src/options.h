// File in HMMDecompositionDecision, created by Thien Le in July 2018

#ifndef OPTION_H
#define OPTION_H

#define NULL_OPTION             -1

#include "c_inc.h"

extern char DEFAULT_FP_OUTPUT_FORMAT[];
extern char DEFAULT_FP_OUTPUT_NAME[];
extern char DEFAULT_FP_INPUT_FORMAT[];
extern char DEFAULT_FP_STDOUT[];

extern fp_options default_fp_options;

// Functions
extern int read_cmd_arg(int argc,char ** argv, option_t * options);
extern int init_options(option_t * options);
extern void destroy_options(option_t * options);

#endif