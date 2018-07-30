// File in inc_ml, created by Thien Le in July 2018

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "tools.h"
#include "utilities.h"

// Private functions
void add_command(char * current_command, char * new_command);

// Function to create command for hmmsearch. 
int fastphylo_job(fp_options * fp_options){
    char command[CMD_BUFFER_SIZE];

    strclr(command);
    strcat(command, "fastdist");
    add_command(command, fp_options->output_format);
    add_command(command, fp_options->output_name);
    add_command(command, fp_options->input_format);
    add_command(command, "-i");
    add_command(command, fp_options->input_name);
    add_command(command, ">");
    add_command(command, fp_options->stdout);

    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling fastphylo job\n", GENERAL_ERROR);

    return 0;
}

int constraint_inc(int argc, option_t * options){
    int i;
    char command[CMD_BUFFER_SIZE];
    char name[CMD_BUFFER_SIZE];
    char num[1000];

    strclr(command);
    strcat(command, "constraint_inc");
    add_command(command, "-i");
    add_command(command, "c_inc_input");
    add_command(command, "-o");
    add_command(command, options->output_name);
    add_command(command, "-t");
    for(i = 0; i < argc; i++){
        strclr(name);
        strclr(num);
        sprintf(num, "%d", i);
        strcat(name, options->
            output_name);
        strcat(name, "_ctree");
        strcat(name, num);
        strcat(name, ".tree");
        add_command(command, name);
    }
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling constraint_inc\n", GENERAL_ERROR);

    return 0;
}

int fasttree_job(option_t * options){
    char command[CMD_BUFFER_SIZE];

    strclr(command);
    strcat(command, "FastTree -nt -gtr");
    add_command(command, "<");
    add_command(command, options->input_name);
    add_command(command, ">");
    add_command(command, options->tree_names[0]);

    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling fasttree_job\n", GENERAL_ERROR);

    return 0;
}

int upp_job(option_t * options){
    char command[CMD_BUFFER_SIZE];

    strclr(command);
    // printf("%s\n", options->input_name);
    strcat(command, "run_upp.py");
    add_command(command, "-a");
    add_command(command, options->input_name);
    add_command(command, "-t");
    add_command(command, "first_tree.tree");
    add_command(command, "-s");
    add_command(command, options->input_name);
    add_command(command, "-A");
    add_command(command, "500");
    add_command(command, "-S");
    add_command(command, "centroid");
    
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling upp_job\n", GENERAL_ERROR);

    return 0;
}

int subset_job(option_t * options){
    char command[CMD_BUFFER_SIZE];

    strclr(command);
    strcat(command, "build_subsets_from_tree.py -t first_tree.tree");
    add_command(command, "-o");
    add_command(command, options->output_name);
    add_command(command, "-n");
    add_command(command, "200");


    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in subset_jon \n", GENERAL_ERROR);
    return 0;
}

int nw_utils_job(option_t * options){
    char command[CMD_BUFFER_SIZE];

    strclr(command);
    strcat(command, "nw_prune");
    add_command(command, options->tree_names[0]);
    add_command(command, "$(cat ");
    strcat(command, options->input_name);
    strcat(command, ")");
    add_command(command, ">");
    add_command(command, options->output_name);
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling nw_utils\n", GENERAL_ERROR);

    return 0;
}

int rm_label_job(option_t * options){
    char command[CMD_BUFFER_SIZE];

    strclr(command);
    strcat(command, "remove_internal_labels.py");
    add_command(command, "-t");
    add_command(command, options->input_name);
    add_command(command, "-o");
    add_command(command, options->output_name);
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling rm_label_job\n", GENERAL_ERROR);

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

