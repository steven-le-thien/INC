// File in inc_ml, created by Thien Le in July 2018

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <unistd.h>

#include "tools.h"
#include "utilities.h"


// Private functions
int find_prefix_and_dir_from_path(char *, char *, char *, char *);

// Public functions 
int constraint_inc(int argc, ml_options * master_ml_options){
  int i, j, read_mode, len;
  char command[GENERAL_BUFFER_SIZE];
  char name[GENERAL_BUFFER_SIZE];

  int arg_count = 6 + argc + (master_ml_options->qtree_method == Q_SUBTREE); 
  char ** argv = malloc(arg_count * sizeof(char*));

  for(i = 0; i < arg_count; i++){
    argv[i] = malloc(10000 * sizeof(char));
    argv[i][0] = 0;
  }
  
  if(master_ml_options->qtree_method == Q_SUBTREE)
    sprintf(
        command, 
        "%s -i %s -o %s -t %s", 
        constraint_inc_bin,
        master_ml_options->use_distance_matrix ? 
            master_ml_options->init_d_name : 
            master_ml_options->input_alignment, 
        master_ml_options->output_prefix, 
        master_ml_options->guide_tree_name
    );
  else if (master_ml_options->ctree_method != C_NO && argc > 0)
    sprintf(
        command, 
        "%s -i %s -o %s -t",
        constraint_inc_bin, 
        master_ml_options->use_distance_matrix ? 
            master_ml_options->init_d_name : 
            master_ml_options->input_alignment, 
        master_ml_options->output_prefix
    );
  else
    sprintf(
        command, 
        "%s -o %s", 
        constraint_inc_bin,
        master_ml_options->output_prefix
    );

  for(i = 0; i < argc; i++){
    sprintf(name, " %s_ctree%d.tree", master_ml_options->output_prefix, i);
    strcat(command, name);
  }

  // Put commands into argv
  j = 0;
  len = 0;
  read_mode = 1;
  for(i = 0; i < strlen(command); i++){
    if(!read_mode){
      if(command[i] == ' ') continue; // still blank
      else{
        argv[j][len] = 0;
        len = 0;
        j++;
        argv[j][len++] = command[i];
        read_mode = 1;
      }
    } else { // still reading
      if(command[i] == ' ')
        read_mode = 0;
      else 
        argv[j][len++] = command[i]; 
    }   
    
  }
  argv[j][len] = 0;

  FCAL(GENERAL_ERROR,
        ERR_CINC,
        constraint_inc_main(arg_count, argv, master_ml_options)
  );

  return 0;
}

int make_constraint_trees_from_disjoint_subsets(
    int n, 
    msa_t * msa,  
    int ** disjoint_subset, 
    ml_options * master_ml_options)
{
  int i, j;
  FILE * f;

  char msa_name[GENERAL_BUFFER_SIZE];
  char out_name[GENERAL_BUFFER_SIZE];

  int sqrt_n = (int) sqrt(1.0 * n);

  for(i = 0; i < sqrt_n; i++){
    // Write the subset msa to a file 
    sprintf(out_name, "%s_ctree%d.tree", master_ml_options->output_prefix, i);
    sprintf(msa_name, "%s_ctree%d.msa", master_ml_options->output_prefix, i);

    // Read in MSA
    f = fopen(msa_name, "w");
    for(j = 0; j < n; j++)
      if(disjoint_subset[j][i])
        fprintf(f, ">%s\n%s\n", msa->name[j], msa->msa[j]);
    fclose(f); 

    // Call fasttree or raxml on the msa
    if(master_ml_options->ctree_method == C_RAXML)
      FCAL(
          GENERAL_ERROR,
          F_RAXML_CONSTRAINT,
          make_raxml_constraint(msa_name, out_name, master_ml_options, 1)
      );
    else if(master_ml_options->ctree_method == C_FASTTREE)
      FCAL(
          GENERAL_ERROR,
          F_FASTTREE_CONSTRAINT,
          make_fasttree_constraint(msa_name, out_name, master_ml_options, 1)
      );
  }

  return 0;
}

int make_unweighted_matrix(
    char * in_tree, 
    char * output_prefix, 
    float ** dm, 
    char ** name_map, 
    int * master_to_midx)                        
{
  FILE * f;
  char output_name[GENERAL_BUFFER_SIZE];
  char name[GENERAL_BUFFER_SIZE];
  int i, j, n;

  // Test output
  sprintf(output_name, "%ssecondary_matrix", output_prefix);
  f = fopen(output_name, "r");

  // Make if not present
  if(!f){

    FCAL(
        GENERAL_ERROR, 
        F_CONSTR_UNW_MAT, 
        unweighted_job(in_tree, output_name)
    );

    f = fopen(output_name, "r");
  }
  
  fscanf(f, "%d", &n);
  for(i = 0; i < n; i++){

    // Locate name and update the mapping to the master indexing scheme
    fscanf(f, "%s", name);
    for(j = 0; j < n; j++){
      if(strcmp(name, name_map[j]) == 0){
        master_to_midx[j] = i;
        break;
      }
    }

    // Copy the actual distance matrix
    for(j = 0; j < n; j++){
      fscanf(f, "%f", &dm[i][j]);
    }
  }
  fclose(f);

  return 0;
}

int make_subset_label(
    char * tree_name, 
    char * out_name, 
    ml_options * master_ml_options)                  
{
  FCAL(
      GENERAL_ERROR,
      F_SUBSET_LABEL,
      subset_job(
          tree_name,
          out_name,
          master_ml_options->ss_threshold,
          master_ml_options->input_alignment
      )
  );
  return 0;
}

int make_subtree(
    char * label, 
    char * out_path, 
    ml_options * master_ml_options, 
    int clear_lab)
{
  FCAL(
      GENERAL_ERROR, 
      F_RMV_LABEL_IN_MAKE_SUBTREE,
      rm_label_job(label, master_ml_options->init_tree_name)
  );

  if(clear_lab)
    FCAL(
        GENERAL_ERROR,
        F_RMV_LABEL_IN_MAKE_SUBTREE,
        nw_utils_job(master_ml_options->init_tree_name, label, out_path)
    );  
  return 0;
}

int make_fasttree_constraint(
    char * msa_name, 
    char * out_name, 
    ml_options * master_ml_options, 
    int clear_lab)
{
  FCAL(
      GENERAL_ERROR, 
      F_FT_IN_MK_FT_CONSTRAINT,
      fasttree_job(
          master_ml_options->distance_model,
          msa_name,
          out_name
      )
  );

  if(clear_lab)
    FCAL(
        GENERAL_ERROR, 
        F_RM_LBL_IN_MK_FT_CONSTRAINT,
        rm_label_job(out_name, out_name)
    );
  return 0;
}

int make_nj_constraint(
    char * msa_name, 
    char * out_name, 
    ml_options * master_ml_options, 
    int clear_lab)
{
  FCAL(
      GENERAL_ERROR, 
      F_NJ_IN_MK_NJ_CONSTRAINT,
      fasttree_job(
          master_ml_options->distance_model,
          msa_name,
          out_name
      )
  );
  if(clear_lab)
    FCAL(
        GENERAL_ERROR, 
        F_RM_LBL_IN_MK_NJ_CONSTRAINT,
        rm_label_job(out_name, out_name)
    );
  return 0;
}

int make_fastme_constraint(
    char * msa_name, 
    char * out_name, 
    ml_options * master_ml_options, 
    int clear_lab)
{
  FCAL(
      GENERAL_ERROR, 
      F_NJ_IN_MK_FASTME_CONSTRAINT,
      fasttree_job(
          master_ml_options->distance_model,
          msa_name,
          out_name
      )
  );
  if(clear_lab)
    FCAL(
        GENERAL_ERROR, 
        F_RM_LBL_IN_MK_FASTME_CONSTRAINT,
        rm_label_job(out_name, out_name)
    );
  return 0;
}

int make_raxml_constraint(
    char * msa_name, 
    char * out_name, 
    ml_options * master_ml_options, 
    int clear_lab)
{
  char        raxml_name[GENERAL_BUFFER_SIZE];
  char        raxml_out_name[GENERAL_BUFFER_SIZE];
  char        raxml_dir_name[GENERAL_BUFFER_SIZE];

  // Find the out_name and dir_name for RAxML (this is required only for RAxML)
  find_prefix_and_dir_from_path(out_name, 
                                raxml_out_name, 
                                raxml_dir_name, 
                                raxml_name);

  FCAL(
      GENERAL_ERROR, 
      F_RXML_IN_MK_RXML_CONSTRAINT,
      raxml_job(master_ml_options->distance_model,
                  raxml_out_name,
                  msa_name,
                  raxml_dir_name)
  );
  if(clear_lab)
    FCAL(
        GENERAL_ERROR, 
        F_RM_LBL_IN_MK_FASTME_CONSTRAINT,
        rm_label_job(raxml_name, out_name)
    );
  return 0;
}

// Jobs
int unweighted_job(char * in_tree, char * out_path){
  SYSCAL(
      GENERAL_ERROR,
      ERR_UNWGHT,
      "%s -t %s > %s",
      tree_to_dist_bin,
      in_tree,
      out_path
  );
  return 0;
}

int fasttree_job(DIST_MOD dist_model, char * in_aln, char * out_path)
{
  SYSCAL(
      GENERAL_ERROR,
      ERR_FT,
      "%s -nt %s -quiet < %s > %s",
      FastTree_bin,
      dist_model == D_JC ? FT_JC : FT_GTRGAMMA,
      in_aln,
      out_path
  );
  return 0;
}

int fastme_job(char * in_aln, char * out_path)
{
  SYSCAL(
      GENERAL_ERROR,
      ERR_FASTA_TO_PHYLIP,
      "%s %s -o > %s",
      fasta_to_phylip_bin,
      in_aln,
      TMP_FILE1
  );

  SYSCAL(
      GENERAL_ERROR,
      ERR_FASTME,
      "%s -i %s -o %s -m B -n -s -d",
      fastme_bin, 
      TMP_FILE1,
      out_path
  );

  SYSCAL(
      GENERAL_ERROR,
      ERR_RM,
      "rm %s", 
      TMP_FILE1
  );
  return 0;
}

int nj_job(DIST_MOD dist_model, char * in_aln, char * out_path)
{
  SYSCAL(
      GENERAL_ERROR,
      ERR_NJ_1,
      "echo \"ToNEXUS format=FASTA fromFile=%s toFile=%s; Quit;\"| %s -n",
      in_aln,
      TMP_FILE1,
      PAUP_bin
  );

  SYSCAL(
      GENERAL_ERROR,
      ERR_NJ_2,
      "echo \"exe %s; NJ distance=%s showtree=No;savetrees"\
        "file=%s format=newick; Quit;\" | %s -n",
      TMP_FILE1,
      dist_model == D_JC ? NJ_JC : NJ_LOGDET, 
      out_path, 
      PAUP_bin
  );
  return 0;
}

int fasttree_initial_tree_job(DIST_MOD dist_model, char * in_aln, char * out_path)
{
  FILE * f, * p;
  char buf[GENERAL_BUFFER_SIZE]; 

  SYSCAL(
      GENERAL_ERROR,
      ERR_FT_INIT,
      "%s -nt %s -quiet -noml -nni 0 -log %s < %s",
      FastTree_bin,
      dist_model == D_JC ? FT_JC : FT_GTRGAMMA,
      TMP_FILE1,
      in_aln
  );

  f = fopen(TMP_FILE1, "r");
  p = fopen(out_path, "w");
  while(fscanf(f, "%s", buf) >=0){
    if(strcmp("NJ", buf) == 0){
      fscanf(f, "%s", buf); 
      fprintf(p, "%s\n", buf);
      break;
    }
  }
  fclose(f);
  fclose(p);

  SYSCAL(
      GENERAL_ERROR, 
      ERR_RM, 
      "rm %s", 
      TMP_FILE1
  );
  return 0;
}

int upp_job(char * in_aln, char * in_tree, char * in_seq, int ss_size)
{
  // TODO: first_tree.tree
  SYSCAL(
      GENERAL_ERROR,
      ERR_UPP,
      "%s -a %s -t %s -s %s -A %d -S centroid", 
      run_upp_bin, 
      in_aln,
      in_tree,
      in_seq,
      ss_size
  );
  return 0;
}

int subset_job(char * in_tree, char * out_path, int ss_size, char * in_aln)
{ 
  SYSCAL(
      GENERAL_ERROR,
      ERR_SUBSET,
      "%s -t %s -o %s -n %d", 
      build_subsets_bin,
      in_tree,
      out_path,
      ss_size
  );
  return 0;
}

int nw_utils_job(char * in_tree, char * in_label, char * out_path)
{
  SYSCAL(
      GENERAL_ERROR,
      ERR_NW_UTILS,
      "%s -v %s $(cat %s) > %s",
      nw_prune_bin, 
      in_tree,
      in_label,
      out_path
  );
  return 0;
}

int rm_label_job(char * in_tree, char * out_path)
{
  SYSCAL(
      GENERAL_ERROR,
      ERR_RM_LBL,
      "%s -t %s -o %s", 
      rm_lbl_bin,
      in_tree, 
      out_path
  );
  return 0;
}

int raxml_job(
    DIST_MOD dist_model, 
    char * out_pfx,    
    char * in_seq,     
    char * wd_sfx)          
{
  // TODO: add length checking

  if(wd_sfx)
    SYSCAL(
        GENERAL_ERROR,
        ERR_RXML,
        "%s -T 12 --silent -m %s -j -n %s -s %s -w %s -p 1 > %s", 
        RAxML_bin, 
        dist_model == D_JC ? RAXML_JC : RAXML_GTRGAMMA,
        out_pfx,
        in_seq,
        wd_sfx,
        TMP_FILE1
    );
  else
    SYSCAL(
        GENERAL_ERROR,
        ERR_RXML,
        "%s -T 12 --silent -m %s -j -n %s -s %s -p 1 > %s",    
        RAxML_bin, 
        dist_model == D_JC ? RAXML_JC : RAXML_GTRGAMMA,
        out_pfx, 
        in_seq,
        TMP_FILE1
    );

  SYSCAL(GENERAL_ERROR, ERR_RM, "rm %s", TMP_FILE1); 
  return 0; 
}

// This also deletes the initial tree 
int raxml_with_initial_tree_job(
    DIST_MOD dist_model,
    char * out_pfx,
    char * in_aln,
    char * out_dir)                            
{
  FCAL(
      GENERAL_ERROR,
      F_FFT_IN_RAXML_W_INIT,
      fasttree_initial_tree_job(dist_model, in_aln, TMP_FILE1) 
  );

  if(out_dir)
    SYSCAL(
        GENERAL_ERROR,
        ERR_RAXML_W_INIT,
        "%s -T 12 --silent -m %s -j -n %s -s %s -w %s -t %s -p 1 > %s", 
        RAxML_bin, 
        dist_model == D_JC ? RAXML_JC : RAXML_GTRGAMMA,
        out_pfx,
        in_aln,
        out_dir,
        TMP_FILE1,
        TMP_FILE2
    );
  else 
    SYSCAL(
        GENERAL_ERROR,
        ERR_RAXML_W_INIT,
        "%s -T 12 --silent -m %s -j -n %s -s %s -t %s -p 1 > %s", 
        RAxML_bin, 
        dist_model == D_JC ? RAXML_JC : RAXML_GTRGAMMA,
        out_pfx,
        in_aln,
        TMP_FILE1,
        TMP_FILE2
    );

  SYSCAL(GENERAL_ERROR, ERR_RM, "rm %s %s", TMP_FILE1, TMP_FILE2); 
  return 0; 
}

int get_ll_from_raxml(
    DIST_MOD dist_model,
    char * out_pfx,
    char * in_aln,
    char * out_dir, 
    char * constraint_quartet, 
    double* ll)
{
  FILE * f; 
  char tmp[GENERAL_BUFFER_SIZE];
  double tmp_d;
  int tmp_i;
  int i;
  sprintf(tmp, "RAxML_info.%s", out_pfx);

  FCAL(
      GENERAL_ERROR,
      F_RXAML_QTREE_IN_GET_LL,
      raxml_with_quartet_tree_job(
          dist_model, 
          out_pfx,
          in_aln,
          out_dir,
          constraint_quartet
      )
  );

  f = fopen(tmp, "r");
  while(fgets(tmp, sizeof(tmp), f))
    if(sscanf(tmp, "%d %lf", &tmp_i, &tmp_d) == 2)
      for(i = 0; i < 3; i++)
        if(tmp_i == i) ll[i] = tmp_d;
    
  fclose(f);

  SYSCAL(GENERAL_ERROR, ERR_RM, "rm RAxML_* %s", "");
  return 0;
}

int raxml_with_quartet_tree_job(
    DIST_MOD dist_model,
    char * out_pfx,
    char * in_aln,
    char * out_dir,
    char * constraint_quartet)
{
  FILE * f;

  f = fopen(TMP_FILE1, "w"); // this will overwrite tmp.fast
  fprintf(f, "%s\n", constraint_quartet); // this may not be secure
  fclose(f);

  if(out_dir)
    SYSCAL(
        GENERAL_ERROR,
        ERR_RAXML_QTREE,
        "%s -T 12 --silent -m %s -j -n %s -s %s -w %s -f N -z %s -p 1 > %s",
        RAxML_bin,  
        dist_model == D_JC ? RAXML_JC : RAXML_GTRGAMMA,
        out_pfx,
        in_aln,
        out_dir,
        TMP_FILE1,
        TMP_FILE2
    );
  else 
    SYSCAL(
        GENERAL_ERROR,
        ERR_RAXML_QTREE,
        "%s -T 12 --silent -m %s -j -n %s -s %s -f N -z %s -p 1 > %s",
        RAxML_bin,
        dist_model == D_JC ? RAXML_JC : RAXML_GTRGAMMA,
        out_pfx,
        in_aln,
        TMP_FILE1,
        TMP_FILE2
    );

  SYSCAL(GENERAL_ERROR, ERR_RM, "rm %s %s", TMP_FILE1, TMP_FILE2); 
  return 0; 
}


int distance_matrix_job(
    char * tmp_folder, 
    DIST_MOD dist_model, 
    char * in_aln, 
    char * out_path)
{
  SYSCAL(
      GENERAL_ERROR,
      ERR_DIST_1,
      "echo \"ToNEXUS format=FASTA fromFile=%s toFile=%s/%s; Quit;\" | %s -n;",
      in_aln, 
      tmp_folder, 
      TMP_FILE1,
      PAUP_bin
  );

  SYSCAL(
      GENERAL_ERROR,
      ERR_DIST_2,
      "echo \"exe %s/%s; DSet distance=%s; SaveDist format=RelPHYLIP"\
                        " file=%s triangle=both diagonal=yes; Quit;\" | %s -n",
      tmp_folder,
      TMP_FILE1,
      dist_model == D_JC ? NJ_JC : NJ_LOGDET,
      out_path,
      PAUP_bin
  );

  SYSCAL(
      GENERAL_ERROR,
      ERR_RM,
      "rm %s/%s",
      tmp_folder,
      TMP_FILE1
  );
  return 0;
}

// Internal functions

int find_prefix_and_dir_from_path(
    char * path, 
    char * prefix, 
    char * dir, 
    char * name)
{                              
  int i, j;
  // Malloc sequence, assuming path is created
  ASSERT(GENERAL_ERROR, N_PATH_IN_FIND_PFX_DIR, path);

  // Init
  STR_CLR(prefix);
  STR_CLR(dir); 

  // Find the last backslash, this separates the dir from the prefix
  for(i = strlen(path) - 1; i >= 0; i--)
    if(path[i] == '/') break;
  
  // Do a string copy
  for(j = 0; j < strlen(path); j++){
    if(j <= i) dir[j] = path[j];
    else prefix[j - i - 1] = path[j];
  }
  dir[i + 1] = 0;
  prefix[strlen(path) - i - 1] = 0;

  sprintf(name, "%sRAxML_bestTree.%s", dir, prefix);

  return 0; 
}
