// File in HMMDecompositionDecision, created by Thien Le in July 2018

#ifndef OPTION_H
#define OPTION_H

#define NULL_OPTION             -1

#include <stdlib.h>

extern int DEFAULT_NUM_OPTIONS;

typedef struct options{
    int num_options;
    int num_trees;

    int input_index;
    char * input_name;

    int output_index;
    char * output_name;

    int tree_index;
    char ** tree_names ;
} option_t;


typedef struct fp{
    char * output_format;
    char * output_name;
    char * input_format;
    char * input_name;
    char * stdout;
} fp_options;

extern char[]   DEFAULT_

extern fp_options default_fp_options;

// Functions
extern int read_cmd_arg(int argc,char ** argv, option_t * options);
extern int init_options(option_t * options);
extern void destroy_options(option_t * options);

#endif