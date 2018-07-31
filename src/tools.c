// File in inc_ml, created by Thien Le in July 2018

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "tools.h"
#include "utilities.h"

int use_true_tree_for_constraint_trees = 1;
int subset_threshold = 200;

// Private functions
void add_command(char * current_command, char * new_command);

// Function to create command for hmmsearch. 
int fastphylo_job(fp_options * fp_options){
    char command[CMD_BUFFER_SIZE];
    sprintf(command, "fastdist %s %s %s -i %s > %s", fp_options->output_format, fp_options->output_name, fp_options->input_format, fp_options->input_name, fp_options->stdout);

    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling fastphylo job\n", GENERAL_ERROR);
    return 0;
}

int constraint_inc(int argc, option_t * options){
    int i;
    char command[CMD_BUFFER_SIZE];
    char name[CMD_BUFFER_SIZE];
    sprintf(command, "constraint_inc -i %sc_inc_input -o %s -t", options->output_name, options->output_name);

    for(i = 0; i < argc; i++){
        sprintf(name, "%s_ctree%d.tree", options->output_name, i);
        add_command(command, name);
    }
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling constraint_inc\n", GENERAL_ERROR);

    return 0;
}

int fasttree_job(option_t * options){
    char command[CMD_BUFFER_SIZE];
    sprintf(command, "FastTree -nt -gtr -quiet < %s > %s%s", options->input_name, options->output_name, options->tree_names[0]);

    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling fasttree_job\n", GENERAL_ERROR);
    return 0;
}

int upp_job(option_t * options){
    char command[CMD_BUFFER_SIZE];
    sprintf("run_upp.py -a %s -t %sfirst_tree.tree -s %s -A %d -S centroid", options->input_name, options->output_name, options->input_name, subset_threshold);
    
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling upp_job\n", GENERAL_ERROR);
    return 0;
}

int subset_job(option_t * options){
    char command[CMD_BUFFER_SIZE];
    if(use_true_tree_for_constraint_trees)
        sprintf(command, "build_subsets_from_tree.py -t %s -o %s -n %d", options->tree_names[1], options->output_name, subset_threshold);
    else 
        sprintf(command, "build_subsets_from_tree.py -t %sfirst_tree.tree -o %s -n %d", options->output_name, options->output_name, subset_threshold);

    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in subset_jon \n", GENERAL_ERROR);
    return 0;
}

int nw_utils_job(option_t * options){
    char command[CMD_BUFFER_SIZE];
    sprintf(command, "nw_prune -v %s $(cat %s) > %s", options->tree_names[0], options->input_name, options->output_name);
        printf("%s\n", command);

    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling nw_utils\n", GENERAL_ERROR);
    return 0;
}

int rm_label_job(option_t * options){
    char command[CMD_BUFFER_SIZE];
    sprintf(command, "remove_internal_labels.py -t %s -o %s", options->input_name, options->output_name);

    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling rm_label_job\n", GENERAL_ERROR);
    return 0;
}

int raxml_job(option_t * options){
    char command[CMD_BUFFER_SIZE];
    sprintf(command, "raxmlHPC-AVX -T 12 -m GTRGAMMA -j -n %s -s %s -w %s -p 1", options->output_name, options->input_name, options->tree_names[0]);

    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling raxml job\n", GENERAL_ERROR);
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

