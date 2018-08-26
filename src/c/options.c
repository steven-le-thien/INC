// File in inc_ml, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "options.h"
#include "utilities.h"

// Internal implementations
int find_arg_index(char * flag, char * content, option_t * options, int * i, int argc, char ** argv);
int init_options(option_t * options);
int init_ml_options(ml_options * options);
int parse_ml_arg(char ** argv, int flag_index, int content_s, int content_e, ml_options * ml_options);
void destroy_options(option_t * options);

/* This function takes a point to some ml_options, popularizes its fields and set the corresponding values depending on the user input
 * It will also catch errors in incompatible calls
 * Input:       user inputs and an empty ml_options struct
 * Output:      0 on success, ERROR otherwise
 * Effect:      popularizes fields in ml_options, may call malloc on stuff
 */ 
int read_ml_cmd_arg(int argc, char ** argv, ml_options * ml_options){
    int i;
    int is_flag;
    int flag_index;
    int content_s, content_e; // these are supposed to be left close, right open 

    // Initializing the options to default fields (which assumes that nothing has been computed before and this is the very first run)
    init_ml_options(ml_options);

    // State machine to read in user input
    if(argv[1][0] != '-')       PRINT_AND_RETURN("the first argument to the binary file should be a flag, please check your input", GENERAL_ERROR);
    is_flag = 1; // always start with a flag
    flag_index = 1;
    for(i = 2; i < argc; i++){
        if(is_flag){ // was reading a flag from the previous round, actual flag is stored in flag_index
            if(argv[i][0] == '-') // read another flag.. something is wrong
                PRINT_AND_RETURN("reading 2 consecutive flags in a row, make sure that your contents does not start with a -", GENERAL_ERROR);
            else{ // reading the first content, note that there may be a lot more contents to read
                content_s = i;
                is_flag = 0;
            }
        } else { // have been reading content
            if(argv[i][0] == '-'){ // read a flag here, so no more content to read
                content_e = i;
                is_flag = 1;

                // Read the previous flag into option
                parse_ml_arg(argv, flag_index, content_s, content_e, ml_options);

                // Reset the machine
                flag_index = i;
                content_s = -1;
                content_e = -1;
            } 
            // Else, reading another content, so we don't need to do anything
        }
    }
    parse_ml_arg(argv, flag_index, content_s, argc, ml_options);

    // May want to double check things here, but for now assume that they are parsed correctly
    return 0; 
}

/* Wrapper function to check for command length and loop through command to parse options
 * Input    argc, argv from command line and an option struct pointer
 * Output   0 on success, ERROR otherwise
 * Effect   set fields in option struct
 */
int read_cmd_arg(int argc, char ** argv, option_t * options){
    if(init_options(options)               != SUCCESS)         PRINT_AND_RETURN("init_options failed in main\n", GENERAL_ERROR);

    if(!options) 
        PRINT_AND_RETURN("options is null in read_cmd_arg", GENERAL_ERROR);

    printf("argc is %d\n", argc);
    for(int i  = 0; i < argc; i++){
        printf("argv is %s\n", argv[i]);
    }
    // Loop through command
    for(int i = 0; i < argc; i++)
        if(i != 0 && i % 2) //currently reading a flag
            if(find_arg_index(argv[i], argv[i + 1], options, &i, argc, argv) != SUCCESS)
                PRINT_AND_RETURN("failed reading of the arguments", GENERAL_ERROR);

    return 0;
}


int parse_ml_arg(char ** argv, int flag_index, int content_s, int content_e, ml_options * ml_options){
    char * flag = argv[flag_index];
    // Switch case for the flag index, which should be a single character, following a dot
    // First check that there is exactly one character that follows the dash 
    if(strlen(flag) != 2 || flag[0] != '-' || flag[1] < 'A' || flag[1] > 'z') PRINT_AND_RETURN("flag name is not correct, please double check your argument into the binary file", GENERAL_ERROR);

    // printf("flag is %s, contetn s is %s content e is %s\n", argv[flag_index], argv[content_s], argv[content_e]);

    // Actually look up what exact field show we be setting and set said field
    // -i: input alignment              (expect the full path - does not check)
    // -o: output prefix                (expect the full path - does not check)
    // -t: initial tree                 (expect the full path - does not check)
    // -g: guide tree                   (expect the full path - does not check)
    // -m: initial distance matrix      (expect the full path - does not check)
    // -r: recompute constraint trees   (expect 0/1) (this make the next option do nothing)
    // -c: constraint trees type        (expect "no"/"subtree"/"fasttree"/"raxml" for no constraint/subtree contraint/fasttree constraint/raxml constraint
    // -q: quartet tree type            (expect "fpm"/"subtree"/"raxml"/"ml" for four point method on original distance matrix/induced quartet from initial tree/raxml to compute each quartet on the fly/other ml methods)
    // -d: distance model               (expect "JC"/"logDet"/"K2P"/"P")
    // -s: max subset size              (expect a positive number - does not check)
    switch(flag[1]){   
        case 'i':
            if(content_e - content_s == 1){
                ml_options->input_alignment = argv[content_s];
            } else PRINT_AND_RETURN("wrong number of argument for input alignment, please check to make sure that you have exactly 1 string following -i", GENERAL_ERROR);
            break;
        case 'o':
            if(content_e - content_s == 1){
                ml_options->output_prefix = argv[content_s];
            } else PRINT_AND_RETURN("wrong number of argument for output_prefix, please check to make sure that you have exactly 1 string following -o", GENERAL_ERROR);
            break;
        case 't':
            if(content_e - content_s == 1){
                ml_options->init_tree_name = argv[content_s];
            } else PRINT_AND_RETURN("wrong number of argument for initial tree, please check to make sure that you have exactly 1 string following -t", GENERAL_ERROR);
            break;
        case 'g':
            if(content_e - content_s == 1){
                ml_options->guide_tree_name = argv[content_s];
            } else PRINT_AND_RETURN("wrong number of argument for guide tree, please check to make sure that you have exactly 1 string following -t", GENERAL_ERROR);
            break;    
        case 'm':
            if(content_e - content_s == 1){
                ml_options->init_d_name = argv[content_s];
            } else PRINT_AND_RETURN("wrong number of argument for init_tree_name, please check to make sure that you have exactly 1 string following -d", GENERAL_ERROR);
            break;
        case 'r':
            if(content_e - content_s == 1 && (argv[content_s][0] == '0' || argv[content_s][0] == '1')){
                ml_options->recompute_constraint_trees = argv[content_s][0] - '0';
            } else PRINT_AND_RETURN("wrong number of argument or wrong type of argument for recompute_constraint_trees, please check to make sure that you have exactly 1 string of 0 or 1 following -r", GENERAL_ERROR);
            break;
        case 'c':
            if(content_e - content_s == 1 && (!strcmp(argv[content_s], "no") || !strcmp(argv[content_s], "subtree") || !strcmp(argv[content_s], "fasttree") || !strcmp(argv[content_s], "raxml"))) {
                switch(argv[content_s][0]){
                    case 'n':
                        ml_options->use_constraint = 0;
                        ml_options->use_subtree_for_constraint_trees    = 0;
                        ml_options->use_raxml_for_constraint_trees      = 0;
                        ml_options->use_fasttree_for_constraint_trees   = 0;
                        break;
                    case 's':
                        ml_options->use_constraint = 1;
                        ml_options->use_subtree_for_constraint_trees    = 1;
                        ml_options->use_raxml_for_constraint_trees      = 0;
                        ml_options->use_fasttree_for_constraint_trees   = 0;
                        break;
                    case 'f':
                        ml_options->use_constraint = 1;
                        ml_options->use_subtree_for_constraint_trees    = 0;
                        ml_options->use_raxml_for_constraint_trees      = 0;
                        ml_options->use_fasttree_for_constraint_trees   = 1;
                        break;
                    case 'r':
                        ml_options->use_constraint = 1;
                        ml_options->use_subtree_for_constraint_trees    = 0;
                        ml_options->use_raxml_for_constraint_trees      = 1;
                        ml_options->use_fasttree_for_constraint_trees   = 0;
                        break;
                    default: 
                        break;
                }
            } else PRINT_AND_RETURN("wrong number of argument or wrong type of argument for constraint trees type, please check to make sure that you have exactly 1 string following -c", GENERAL_ERROR);
            break;
        case 'q':
            if(content_e - content_s == 1 && (!strcmp(argv[content_s], "fpm") || !strcmp(argv[content_s], "subtree") || !strcmp(argv[content_s], "raxml") || !strcmp(argv[content_s], "ml"))) {
                switch(argv[content_s][0]){
                    case 'f':
                        ml_options->use_four_point_method_with_distance     = 1;
                        ml_options->use_four_point_method_with_tree         = 0;
                        ml_options->use_new_quartet_raxml                   = 0;
                        ml_options->use_ml_method                           = 0;
                        break;
                    case 's':
                        ml_options->use_four_point_method_with_distance     = 0;
                        ml_options->use_four_point_method_with_tree         = 1;
                        ml_options->use_new_quartet_raxml                   = 0;
                        ml_options->use_ml_method                           = 0;
                        break;
                    case 'r':
                        ml_options->use_four_point_method_with_distance     = 0;
                        ml_options->use_four_point_method_with_tree         = 0;
                        ml_options->use_new_quartet_raxml                   = 1;
                        ml_options->use_ml_method                           = 0;
                        break;
                    case 'm':
                        ml_options->use_four_point_method_with_distance     = 0;
                        ml_options->use_four_point_method_with_tree         = 0;
                        ml_options->use_new_quartet_raxml                   = 0;
                        ml_options->use_ml_method                           = 1;
                        break;
                    default: 
                        break;
                }
            } else PRINT_AND_RETURN("wrong number of argument or wrong type of argument for quartet tree type, please check to make sure that you have exactly 1 string following -q", GENERAL_ERROR);
            break;
        case 'd':
            if(content_e - content_s == 1 && (!strcmp(argv[content_s], "JC") || !strcmp(argv[content_s], "P") || !strcmp(argv[content_s], "K2P") || !strcmp(argv[content_s], "logDet"))) {
                ml_options->distance_model = argv[content_s];
            } else PRINT_AND_RETURN("wrong number of argument for distance model, please check to make sure that you have exactly 1 string following -d", GENERAL_ERROR);
            break;
        case 's':
            if(content_e - content_s == 1){
                ml_options->ss_threshold = atoi(argv[content_s]);
            } else PRINT_AND_RETURN("wrong number of argument for max subset size, please check to make sure that you have exactly 1 string following -s", GENERAL_ERROR);
            break;
        default:
            PRINT_AND_RETURN("does not recognize flag, please check to make sure that there are only valid flag", GENERAL_ERROR);
    } 
    return 0;
}

/* Function to check the flag tag and assign its content to the appropriate field in option structure
 * This function also stores the index in argv to the content
 * Input    the flag, the content and pointer to an option struct as well as index to argv
 * Output   0 on success, ERROR otherwise
 * Effect   set fields in option struct
 */
int find_arg_index(char * flag, char * content, option_t * options, int * i, int argc, char ** argv){
    int j;
    if(strcmp(flag, "-i") == 0){//reading input name
        options->input_index = (*i);
        options->input_name = malloc(strlen(content) + 1);  

        if(!options->input_name)        
            PRINT_AND_RETURN("malloc failure for input name in find_arg_index",     MALLOC_ERROR);
        else{
            strcpy(options->input_name, content);
            options->input_name[strlen(content)] = 0;
        }

    } else if(strcmp(flag, "-o") == 0){
        options->output_index = (*i);
        options->output_name = malloc(strlen(content) + 1);

        if(!options->output_name)       
            PRINT_AND_RETURN("malloc failure for output name in find_arg_index",    MALLOC_ERROR);
        else{
            strcpy(options->output_name, content);
            options->output_name[strlen(content)] = 0;
        }

    } else if(strcmp(flag, "-t") == 0){
        options->tree_index = (*i);
        options->num_trees = 0;

        while(1){
            if((*i) + options->num_trees + 1 >= argc || argv[(*i) + options->num_trees + 1][0] == '-') break;
            else options->num_trees++;
        }

        options->tree_names = malloc(options->num_trees * sizeof(char*));
        for(j = 0; j < options->num_trees; j++){
            options->tree_names[j] = malloc(strlen(argv[(*i) + j + 1]) + 1);
            strcpy(options->tree_names[j], argv[(*i) + j + 1]);
            options->tree_names[j][strlen(argv[(*i) + j + 1])] = 0;
        }
        (*i) += options->num_trees ;

    } else {
        printf("%s\n", flag);
        PRINT_AND_RETURN("unrecognized argument", GENERAL_ERROR); 
    }

    return 0;
}

/* Constructor for option struct
 * Input    pointer to option struct
 * Output   0 on success
 * Effect   set fields in option struct to default value
 */
int init_options(option_t * options){
    options->num_options = 3; // there are 3 main blocks in the option_t struct
    options->num_trees = 0;

    options->input_index = -1;
    options->output_index = -1;
    options->tree_index = -1;

    options->input_name = NULL;
    options->output_name = NULL;
    options->tree_names = NULL;

    return 0;
}

int init_ml_options(ml_options * ml_options){
    // The first two should be initialized properly
    ml_options->input_alignment         = NULL;
    ml_options->output_prefix           = NULL;
    ml_options->init_tree_name          = NULL;
    ml_options->init_d_name             = NULL;

    ml_options->use_constraint                      = 1;
    ml_options->recompute_constraint_trees          = 1;

    ml_options->use_subtree_for_constraint_trees    = 0;
    ml_options->use_raxml_for_constraint_trees      = 0;
    ml_options->use_fasttree_for_constraint_trees   = 1;

    ml_options->use_four_point_method_with_distance     = 1;
    ml_options->use_four_point_method_with_tree         = 0;
    ml_options->use_new_quartet_raxml                   = 0;
    ml_options->use_ml_method                           = 0;

    // These 2 needs to be initialized properly
    ml_options->distance_model                          = NULL;
    ml_options->ss_threshold                            = -1;

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
    for(i = 0; i < options->num_trees; i++){
        free(options->tree_names[i]);
    }
    // if(options->tree_names)        free(options->symfrac);

    // options->input_index = options->output_index = options->symfrac_index = 0;
}