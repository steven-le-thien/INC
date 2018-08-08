// File in inc_ml, created by Thien Le in July 2018

#ifndef TOOLS_H
#define TOOLS_H

#include "c_inc.h"

extern int constraint_inc(int argc, option_t * options);
extern int fasttree_job(option_t * options);
extern int upp_job(option_t * options);
extern int subset_job(option_t * options);
extern int nw_utils_job(option_t * options);
extern int rm_label_job(option_t * options);
extern int distance_matrix_job(option_t * options);
extern int raxml_job(option_t * options);

// Wrappers
extern int make_subtree_label(char * tree_name, char * out_name);
extern int make_subtree(char * label, char * outname, char * tree_name); 
extern int make_fasttree_constraint(char * msa_name, char * out_name);
extern int make_raxml_constraint(char * input_path, char * msa_name, char * out_name); 

#endif