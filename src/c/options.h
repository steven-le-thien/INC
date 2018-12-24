// File in inc_ml, created by Thien Le in July 2018

#ifndef OPTION_H
#define OPTION_H

#include "c_inc.h"

// Functions
extern int read_cmd_arg(int argc,char ** argv, option_t * options);
extern int read_ml_cmd_arg(int argc,char ** argv, ml_options * options);

#endif