// File in inc_ml, created by Thien Le in July 2018

#ifndef TOOLS_H
#define TOOLS_H

#include "c_inc.h"

extern int constraint_inc(int argc, ml_options * master_ml_options);
extern int fasttree_job(option_t * options, ml_options * master_ml_options);
extern int upp_job(option_t * options, ml_options * master_ml_options);
extern int subset_job(option_t * options, ml_options * master_ml_options);
extern int nw_utils_job(option_t * options);
extern int rm_label_job(option_t * options);
extern int distance_matrix_job(option_t * options, ml_options * master_ml_options);
extern int raxml_job(option_t * options, ml_options * master_ml_options);
extern int unweighted_job(option_t * options);
extern int raxml_with_initial_tree_job(option_t * options,  ml_options * master_ml_options);
extern int fasttree_initial_tree_job(option_t * tmp_options, ml_options * master_ml_options);
extern int get_ll_from_raxml(option_t * options,  ml_options * master_ml_options, char * constraint_quartet, double* lp, double* l1, double* l2);
extern int raxml_with_quartet_tree_job(option_t * options,  ml_options * master_ml_options, char * constraint_quartet);

// Wrappers
extern int make_subset_label(char * tree_name, char * out_name, ml_options * master_ml_options);
extern int make_subtree(char * label, char * outname, char * tree_name); 
extern int make_fasttree_constraint(char * msa_name, char * out_name, ml_options * master_ml_options);
extern int make_raxml_constraint(char * msa_name, char * out_name, ml_options * master_ml_options); 
extern int construct_unweighted_matrix_job(char * filename, char * output_prefix, float ** dm, char ** name_map, int * master_to_midx);
extern int make_constraint_trees_from_disjoint_subsets(int n, msa_t * msa,  int ** disjoint_subset, ml_options * master_ml_options);
#endif