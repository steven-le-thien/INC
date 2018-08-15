// File in inc_ml, created by Thien Le in July 2018

#ifndef OPTION_H
#define OPTION_H

#include "c_inc.h"

// Functions
extern int read_cmd_arg(int argc,char ** argv, option_t * options);
extern int init_options(option_t * options);
extern void destroy_options(option_t * options);

#endif