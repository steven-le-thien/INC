// File in HMMDecompositionDecision, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "options.h"
#include "utilities.h"


// Constants
int DEFAULT_NUM_OPTIONS = 3;

char * DEFAULT_OUTPUT_FORMAT    = "-O phylip";
char * DEFAULT_OUTPUT_NAME      = "defaultjob.distance_matrix.phy";
char * DEFAULT_INPUT_FORMAT     = "-I fasta";
char * DEFAULT_STDOUT           = "defaultjob.fastphylo.stdout"

fp_options default_fp_options = {DEFAULT_OUTPUT_FORMAT, DEFAULT_OUTPUT_NAME, DEFAULT_INPUT_FORMAT, NULL, DEFAULT_STDOUT};

/* Function to check the flag tag and assign its content to the appropriate field in option structure
 * This function also stores the index in argv to the content
 * Input    the flag, the content and pointer to an option struct as well as index to argv
 * Output   0 on success, ERROR otherwise
 * Effect   set fields in option struct
 */
int find_arg_index(char * flag, char * content, option_t * options, int i, int argc, char ** argv){
    int j;
    if(strcmp(flag, "-i") == 0){//reading input name
        options->input_index = i;
        options->input_name = malloc(strlen(content));  

        if(!options->input_name)        
            PRINT_AND_RETURN("malloc failure for input name in find_arg_index",     MALLOC_ERROR);
        else
            strcpy(options->input_name, content);

    } else if(strcmp(flag, "-o") == 0){
        options->output_index = i;
        options->output_name = malloc(strlen(content));

        if(!options->output_name)       
            PRINT_AND_RETURN("malloc failure for output name in find_arg_index",    MALLOC_ERROR);
        else
            strcpy(options->output_name, content);

    } else if(strcmp(flag, "-t") == 0){
        options->tree_index = i;
        options->num_tree = 0;

        while(1){
            if(i + options->num_tree + 1 >= argc || argv[i + options->num_tree + 1][0] == '-') break;
            else options->num_tree++;
        }

        options->tree = malloc(options->num_tree * sizeof(char*));
        for(j = 0; j < options->num_tree; j++){
            options->tree[j] = malloc(strlen(argv[i + j + 1]));
            strcpy(options->tree[j], &argv[i + j + 1]);
        }

    } else PRINT_AND_RETURN("unrecognized argument", GENERAL_ERROR); 

    return 0;
}

/* Wrapper function to check for command length and loop through command to parse options
 * Input    argc, argv from command line and an option struct pointer
 * Output   0 on success, ERROR otherwise
 * Effect   set fields in option struct
 */
int read_cmd_arg(int argc, char ** argv, option_t * options){
    if(!options) 
        PRINT_AND_RETURN("options is null in read_cmd_arg", GENERAL_ERROR);

    // Loop through command
    for(int i = 0; i < argc; i++)
        if(i != 0 && i % 2) //currently reading a flag
            if(find_arg_index(argv[i], argv[i + 1], options, i, argc, argv) != SUCCESS)
                PRINT_AND_RETURN("failed reading of the arguments", GENERAL_ERROR);

    return 0;
}

/* Constructor for option struct
 * Input    pointer to option struct
 * Output   0 on success
 * Effect   set fields in option struct to default value
 */
int init_options(option_t * options){
    options->num_options = DEFAULT_NUM_OPTIONS;
    options->num_tree = 0;

    options->input_index = -1;
    options->output_index = -1;
    options->tree_index = -1;

    options->input_name = NULL;
    options->output_name = NULL;
    options->tree_name = NULL;

    return 0;
}

/* Destructor for option struct
 * Input    pointer to option struct
 * Output   0 on success, ERROR otherwise
 * Effect   set fields in option struct
 */
void destroy_options(option_t * options){
    int i;

    // If the struct passed in is not allocated then nothing happens
    if(!options) return;
    if(options->input_name)     free(options->input_name);
    if(options->output_name)    free(options->output_name);
    for(i = 0; i < num_tree; i++){
        free(options->tree_name[i])
    }
    if(options->tree_name)        free(options->symfrac);

    options->input_index = options->output_index = options->symfrac_index = 0;
}