// File in inc_ml, created by Thien Le in July 2018

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "tools.h"
#include "utilities.h"

char distance_model[] = DISTANCE_MODEL;
int subset_threshold = SS_THRESHOLD;

// Private functions
void add_command(char * current_command, char * new_command);
int find_prefix_and_dir_from_path(char * path, char * prefix, char * dir);

int make_subset_label(char * tree_name, char * out_name){
    option_t tmp_options; 

    tmp_options.input_name = tree_name;
    tmp_options.output_name = out_name;
    if(subset_job(&tmp_options)                 != SUCCESS)         PRINT_AND_EXIT("subset job failed in main\n", GENERAL_ERROR);

    return 0;
}

int make_subtree(char * label, char * outname, char * tree_name){
    option_t tmp_options; 
    tmp_options.input_name = tree_name;
    tmp_options.output_name = tree_name; // replacing file
    if(rm_label_job(&tmp_options) != SUCCESS) PRINT_AND_RETURN("remove label failed in make_subtree\n", GENERAL_ERROR);

    tmp_options.input_name = label;
    tmp_options.output_name = outname;
    tmp_options.tree_names = &tree_name;
    if(nw_utils_job(&tmp_options) != SUCCESS) PRINT_AND_RETURN("nw_utils failed in make_subtree\n", GENERAL_ERROR);

    return 0;
}

int make_fasttree_constraint(char * msa_name, char * out_name){
    option_t tmp_options; 

    // Call fasttree
    tmp_options.input_name = msa_name;
    tmp_options.tree_names = &out_name;

    if(fasttree_job(&tmp_options) != SUCCESS) PRINT_AND_RETURN("fasttree failed in main", GENERAL_ERROR);

    tmp_options.input_name = out_name;
    tmp_options.output_name = out_name;
    if(rm_label_job(&tmp_options) != SUCCESS) PRINT_AND_RETURN("remove label failed in main", GENERAL_ERROR);

    return 0;
}

int make_raxml_constraint(char * input_path, char * msa_name, char * out_name){
    option_t tmp_options; 

    char        raxml_name[10000];
    char        raxml_out_name[10000];
    char        raxml_dir_name[10000];

    // Find the out_name and dir_name for RAxML (this is required only for RAxML)
    find_prefix_and_dir_from_path(out_name, raxml_out_name, raxml_dir_name);

    tmp_options.input_name = msa_name;
    tmp_options.output_name = raxml_dir_name;
    tmp_options.tree_names = malloc(sizeof(char *));
    tmp_options.tree_names[0] = raxml_out_name;
    if(raxml_job(&tmp_options) != SUCCESS) PRINT_AND_RETURN("raxml_job failed in main", GENERAL_ERROR);
    
    sprintf(raxml_name, "%sRAxML_bestTree.%s", raxml_dir_name, raxml_out_name);
    tmp_options.input_name = raxml_name;
    tmp_options.output_name = out_name;
    if(rm_label_job(&tmp_options) != SUCCESS) PRINT_AND_RETURN("remove label failed in main", GENERAL_ERROR);

    return 0;
}



int constraint_inc(int argc, option_t * options){
    int i;
    char command[GENERAL_BUFFER_SIZE];
    char name[GENERAL_BUFFER_SIZE];
#if use_induced_quartet && use_induced_fasttree
    sprintf(command, "constraint_inc -i %sc_inc_input -o %s -t %sfirst_tree.tree", options->output_name, options->output_name, options->output_name);
#else 
    sprintf(command, "constraint_inc -i %sc_inc_input -o %s -t", options->output_name, options->output_name);
#endif

    for(i = 0; i < argc; i++){
        sprintf(name, "%s_ctree%d.tree", options->output_name, i);
        add_command(command, name);
    }
    printf("%s\n", command);
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling constraint_inc\n", GENERAL_ERROR);

    return 0;
}

int fasttree_job(option_t * options){
    char command[GENERAL_BUFFER_SIZE];
    sprintf(command, "FastTree -nt -gtr -quiet < %s > %s", options->input_name, options->tree_names[0]);
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling fasttree_job\n", GENERAL_ERROR);
    return 0;
}

int upp_job(option_t * options){
    char command[GENERAL_BUFFER_SIZE];
    sprintf("run_upp.py -a %s -t %sfirst_tree.tree -s %s -A %d -S centroid", options->input_name, options->output_name, options->input_name, subset_threshold);
    
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling upp_job\n", GENERAL_ERROR);
    return 0;
}

int subset_job(option_t * options){ 
    char command[GENERAL_BUFFER_SIZE];
    sprintf(command, "build_subsets_from_tree.py -t %s -o %s -n %d", options->input_name, options->output_name, subset_threshold);

    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in subset_jon \n", GENERAL_ERROR);
    return 0;
}

int nw_utils_job(option_t * options){
    char command[GENERAL_BUFFER_SIZE];
    sprintf(command, "nw_prune -v %s $(cat %s) > %s", options->tree_names[0], options->input_name, options->output_name);
    printf("%s\n", command);
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling nw_utils\n", GENERAL_ERROR);
    return 0;
}

int rm_label_job(option_t * options){
    char command[GENERAL_BUFFER_SIZE];
    sprintf(command, "remove_internal_labels.py -t %s -o %s", options->input_name, options->output_name);
    printf("%s\n", command);
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling rm_label_job\n", GENERAL_ERROR);
    return 0;
}

int raxml_job(option_t * options){
    char command[GENERAL_BUFFER_SIZE];
    sprintf(command, "raxmlHPC-AVX2 -T 12 -m GTRGAMMA -j -n %s -s %s -w %s -p 1",  options->tree_names[0], options->input_name,options->output_name);
    printf("%s\n", command);
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling raxml job\n", GENERAL_ERROR);
    return 0; 
}

int distance_matrix_job(option_t * options){
    char command[GENERAL_BUFFER_SIZE];

    // This is too long
    sprintf(command, "echo \"ToNEXUS format=FASTA fromFile=%s toFile=~/nexus; Quit;\" | paup4a163_osx -n;", options->input_name);
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling distance matrix job\n", GENERAL_ERROR);

    sprintf(command, "echo \"exe ~/nexus; DSet distance=%s; SaveDist format=RelPHYLIP file=%sc_inc_input triangle=both diagonal=yes; Quit;\" | paup4a163_osx -n", distance_model, options->output_name);
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling distance matrix job\n", GENERAL_ERROR);
    return 0;
}

// Internal functions

int find_prefix_and_dir_from_path(char * path, char * prefix, char * dir){
    int i, j;
    // Malloc sequence, assuming path is created
    if(!path)                               PRINT_AND_RETURN("path is NULL in find_prefix_and_dir_from_path\n", GENERAL_ERROR);

    // Init
    prefix[0] = 0;
    dir[0] = 0;

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
    return 0; 
}


/* Helper function to add a command into a string of commands
 * Input:   current string of command and the new command to be added 
 * Output:  none
 * Effect   none
 */
void add_command(char * current_command, char * new_command){
    str_add_space(current_command);
    strcat(current_command, new_command);
}

