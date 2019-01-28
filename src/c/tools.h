// File in inc_ml, created by Thien Le in July 2018

#ifndef TOOLS_H
#define TOOLS_H

#include "c_inc.h"

extern int make_constraint_trees_from_disjoint_subsets(
    int n, 
    msa_t * msa,  
    int ** disjoint_subset, 
    ml_options * master_ml_options                    
);

extern int make_unweighted_matrix(
    char * filename, 
    char * output_prefix, 
    float ** dm, 
    char ** name_map, 
    int * master_to_midx              
);

extern int make_subset_label(
    char * tree_name, 
    char * out_name, 
    ml_options * master_ml_options
);

extern int make_subtree(
    char * label, 
    char * outname, 
    ml_options * master_ml_options,
    int clear_lab
);

extern int make_fasttree_constraint(
    char * msa_name, 
    char * out_name, 
    ml_options * master_ml_options,
    int clear_lab
);

extern int make_nj_constraint(
    char * msa_name, 
    char * out_name, 
    ml_options * master_ml_options,
    int clear_lab
);

extern int make_fastme_constraint(
    char * msa_name, 
    char * out_name,
    ml_options * master_ml_options,
    int clear_lab
);

extern int make_raxml_constraint(
    char * msa_name, 
    char * out_name, 
    ml_options * master_ml_options,
    int clear_lab
);

extern int constraint_inc(
    int argc, 
    ml_options * master_ml_option
);

extern int unweighted_job(
    char * in_tree, 
    char * out_path
);

extern int fasttree_job(
    DIST_MOD dist_model, 
    char * in_aln, 
    char * out_path
); 


extern int raxml_job(
    DIST_MOD dist_model, 
    char * out_pfx,    
    char * in_seq,     
    char * wd_sfx          
);

extern int fastme_job(
    char * in_aln, 
    char * out_path
);

extern int nj_job(
    DIST_MOD dist_model, 
    char * in_aln, 
    char * out_path
);

extern int fasttree_initial_tree_job(
    DIST_MOD dist_model, 
    char * in_aln, 
    char * out_path
); 

extern int upp_job(char * in_aln, char * in_tree, char * in_seq, int ss_size);

extern int subset_job(
    char * in_tree, 
    char * out_path, 
    int ss_size, 
    char * in_aln
); 

extern int nw_utils_job(char * in_tree, char * in_label, char * out_path);

extern int rm_label_job(char * in_tree, char * out_path);


extern int raxml_with_initial_tree_job(
    DIST_MOD dist_model,
    char * out_pfx,
    char * in_aln,
    char * out_dir
);

extern int get_ll_from_raxml(
    DIST_MOD dist_model,
    char * out_pfx,
    char * in_aln,
    char * out_dir, 
    char * constraint_quartet, 
    double* ll
);

extern int raxml_with_quartet_tree_job(
    DIST_MOD dist_model,
    char * out_pfx,
    char * in_aln,
    char * out_dir,
    char * constraint_quartet
);                        


extern int distance_matrix_job(
    DIST_MOD dist_model, 
    char * in_aln, 
    char * out_path
);                            

// Failure string
static const char F_RAXML_CONSTRAINT[]    
  = "make_raxml_constraint failed in make_constraint_trees_from_disjoint_"\
    "subsets\n"; 

static const char F_FASTTREE_CONSTRAINT[] 
  = "make_fasttree_constraint failed in make_constraint_trees_from_disjoint_"\
    "subsets\n";

static const char F_CONSTR_UNW_MAT[]    
  = "construct_unweighted_matrix_job failed\n";

static const char F_SUBSET_LABEL[] 
  = "subset job failed in main\n"; 

static const char F_RMV_LABEL_IN_MAKE_SUBTREE[] 
  = "remove label failed in make_subtree\n";

static const char F_NW_UTILS_IN_MAKE_SUBTREE[]
  = "nw_utils failed in make_subtree\n";

static const char F_FT_IN_MK_FT_CONSTRAINT[] 
  = "fasttree failed in make_fasttree_constraint\n";

static const char F_RM_LBL_IN_MK_FT_CONSTRAINT[] 
  = "remove label failed in make_fasttree_constraint\n";

static const char F_NJ_IN_MK_NJ_CONSTRAINT[] 
  = "nj failed in make_fasttree_constraint\n";

static const char F_RM_LBL_IN_MK_NJ_CONSTRAINT[]
  = "rm_label_job failed in make_nj_constraint\n";

static const char F_NJ_IN_MK_FASTME_CONSTRAINT[]
  = "nj failed in make_fastme_constraint\n"; 

static const char F_RM_LBL_IN_MK_FASTME_CONSTRAINT[]
  = "rm_label_job in make_fastme_constraint\n";

static const char F_RXML_IN_MK_RXML_CONSTRAINT[]
  = "raxml failed in make_raxml_constraint\n";

static const char F_RM_LBL_IN_MK_RXML_CONSTRAINT[] 
  = "rm_label_job failed in make_raxml_constraint\n"; 

static const char F_FFT_IN_RAXML_W_INIT[]
  = "fast fasttree failed in raxml_with_initial_tree_job\n";

static const char F_RXAML_QTREE_IN_GET_LL[]
  = "raxml_with_quartet_tree_job failed in get_ll_from_raxml\n";

// Error in system
static const char ERR_RXML[]    
  = "error in calling raxml job\n"; 
static const char ERR_RM_LBL[] 
  = "error in calling rm_label_job\n";
static const char ERR_NW_UTILS[] 
  = "error in calling nw_utils\n"; 
static const char ERR_SUBSET[]   
  = "error in subset_job \n";
static const char ERR_UPP[]       
  = "error in upp\n"; 
static const char ERR_FT[]      
  = "error in fasttree\n";
static const char ERR_FASTA_TO_PHYLIP[]
  = "error in fasta_to_phylip\n";

static const char ERR_FASTME[]  
  = "error in fastme\n";
static const char ERR_RM[]        
  = "error in rm\n";
static const char ERR_NJ_1[]      
  = "error in nj part 1\n";
static const char ERR_NJ_2[]        
  = "error in nj part 2\n"; 

static const char ERR_CINC[]      
  = "error in constraint_inc\n";
static const char ERR_UNWGHT[]    
  = "error in unweighted_job\n";

static const char ERR_FT_INIT[]     
  = "error in fasttree_initial_tree_job\n";
static const char ERR_RAXML_W_INIT[] 
  = "error in raxml_with_initial_tree_job\n";
static const char ERR_RAXML_QTREE[]
  = "error in raxml_with_quartet_tree_job\n";
static const char ERR_DIST_1[]    
  = "error in distance_matrix_job_part1\n";
static const char ERR_DIST_2[]    
  = "error in distance_matrix_job part2\n";

// Assertion
static const char N_PATH_IN_FIND_PFX_DIR[]
  = "null path in find prefix dir\n";

// Settings
static const char RAXML_GTRGAMMA[] 
  = "GTRGAMMA";
static const char RAXML_JC[]      
  = "GTRCAT -V --JC69";
static const char FT_GTRGAMMA[]     
  = "-gtr -gamma";
static const char FT_JC[]           
  = "";   
static const char NJ_LOGDET[]     
  = "logDet";
static const char NJ_JC[]           
  = "JC";
static const char IN_LOGDET[]     
  = "logDet";
static const char IN_JC[]           
  = "JC"; 


// Bins
#if 1
static const char RAxML_bin[]           
  = "raxmlHPC-AVX2";
static const char FastTree_bin[]        
  = "FastTree"; 
static const char PAUP_bin[]              
  = "paup4a164_osx";
static const char constraint_inc_bin[]    
  = "constraint_inc";
static const char tree_to_dist_bin[]    
  = "tree_to_dist.py";
static const char fasta_to_phylip_bin[]   
  = "fasta_to_phylip.py";
static const char fastme_bin[]            
  = "fastme";
static const char run_upp_bin[]         
  = "run_upp.py";
static const char build_subsets_bin[]    
  = "build_subsets_from_tree.py";
static const char nw_prune_bin[]        
  = "nw_prune";
static const char rm_lbl_bin[]          
  = "remove_internal_labels.py";
#else
static const char RAxML_bin[]    = "raxmlHPC-PTHREADS-SSE3";
static const char FastTree_bin[] = "FastTree";
static const char PAUP_bin[]     = "paup4a163_centos64";
static const char constraint_inc_bin[]    = "constraint_inc";
static const char tree_to_dist_bin[]    = "tree_to_dist.py";
static const char fasta_to_phylip_bin[]   = "fasta_to_phylip.py";
static const char fastme_bin[]            = "fastme";
static const char run_upp_bin[]         = "run_upp.py";
static const char build_subsets_bin[]     = "build_subsets_from_tree.py";
static const char nw_prune_bin[]        = "nw_prune";
static const char rm_lbl_bin[]          = "remove_internal_labels.py";
#endif

#endif
