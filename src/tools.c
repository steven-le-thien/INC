#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "tools.h"
#include "utilities.h"

// Private functions
void add_command(char * current_command, char * new_command);

// Function to create command for hmmsearch. 
int fastphylo_job(fastphylo_job * fp_options){
    char command[CMD_BUFFER_SIZE];

    strclr(command);
    strcat(command, "fastdist");
    add_command(command, fp_options->output_format);
    add_command(command, fp_options->output_name);
    add_command(command, fp_options->input_format);
    add_command(command, '-i');
    add_command(command, fp_options->input_name);
    add_command(command, ">");
    add_command(command, fp_options->stdout);

    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling hmmsearch job\n", GENERAL_ERROR);

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

