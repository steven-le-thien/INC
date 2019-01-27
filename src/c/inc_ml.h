// File in inc_ml, created by Thien Le in Jan 2019

#ifndef INC_ML_H
#define INC_ML_H

// State strings
static const char STATE_INIT_TREE[] 
  = "cheking for initial tree...\n";

static const char STATE_PASTA_DEC[] 
  = "performing PASTA decomposition\n";

static const char STATE_DIST[]
  = "writing distance matrix using PAUP* to ";

static const char STATE_RD_OPT[]
  = "reading in options...\n";

static const char STATE_PARSE_MAT_TREE[]
  = "parsing initial matrix and tree...\n";

static const char STATE_NO_DIST[]
  = "without distance matrix...\n";

static const char STATE_MST[]
  = "computing the mst...\n";

static const char STATE_INIT_GTREE[]
  = "initializing the growing tree...\n"; 

static const char STATE_BUILD_TREE[]
  = "building the tree...\n";

static const char STATE_OUT_TREE[]
  = "outputing the tree...\n";

static const char STATE_CLEAN[]
  = "done, cleaning up\n";

static const char ITER_COUNT[] 
  = "current iteration is 2...";

// Failure strings
static const char F_MAKE_SS_LBL_IN_MAKE_CTREE[] 
  = "make_subset_label failed in make_constraint_trees\n";

static const char F_PARSE_IN_IN_MAKE_CTREE[]
  = "parse_input failed in make_constraint_trees\n";

static const char F_SS_MSA_IN_MAKE_CTREE[]
  = "subset_msa failed in make_constraint_trees\n";

static const char F_MAKE_ST_IN_MAKE_CTREE[]
  = "make_subtree failed in make_constraint_trees\n";

static const char F_MAKE_RAXML_IN_MAKE_CTREE[]
  = "make_raxml_constraint failed in make_constraint_trees\n";

static const char F_MAKE_FT_IN_MAKE_CTREE[]
  = "make_fasttree_constraint failed in make_constraint_trees\n";

static const char F_MAKE_NJ_IN_MAKE_CTREE[] 
  = "make_nj_constraint failed in make_constraint_trees\n";

static const char F_MAKE_FASTME_IN_MAKE_CTREE[]
  = "make_fastme_constraint failed in make_constraint_trees\n";

static const char F_CINC_IN_MAIN[] 
  = "constrained_inc failed in main\n";

static const char F_RD_CMD_ARG[]
  = "read_cmd_arg failed\n";

static const char F_SET_UP_TREES_IN_MAIN[]
  = "set_up_trees_for_in failed in main\n";

static const char F_SET_UP_DIST_IN_MAIN[] 
  = "set_up_dist failed in main\n";

static const char F_SET_UP_I_TREE_IN_SET_UP_TREE[]
  = "set_up_init_tree failed in set_up_trees\n";

static const char F_MAKE_CTREE_IN_SET_UP_TREE[]   
  = "make_ctree failed in set_up_trees\n";

static const char F_ITREE_IN_MAIN[] 
  = "itree failed in main\n";

static const char F_DIST_MAT_IN_MAIN[]
  = "distance matrix failed in main\n";

// static const char STATE_PASTA_DEC[]
//   = "starting PASTA decompositions...\n";

#endif