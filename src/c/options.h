// File in inc_ml, created by Thien Le in July 2018

#ifndef OPTION_H
#define OPTION_H

#include "c_inc.h"

// Functions
extern int read_cmd_arg(int argc,char ** argv, ml_options * options);
extern int read_ml_cmd_arg(int argc,char ** argv, ml_options * options);
extern int init_ml_options(ml_options * options);


static const char WRONG_ARG_FORMAT[]
  = "wrong argument format, please check your commands\n";

static const char F_FIND_ARG_INC_OPT[]
  = "find_arg_idx failed in inc options\n";

static const char F_PARSE_ML_ARG[]
  = "parse_ml_arg failed in read_ml_cmd_arg\n";

static const char F_INIT_ML_OPT[] 
  = "init_ml_options failed\n";

static const char M_ERR_IN_FIND_ARG_IDX[]
  = "malloc error in find_arg_idx\n";

#endif