// File in inc_ml, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "options.h"
#include "utilities.h"

// Internal implementations
int find_arg_index(
    char * flag, 
    char * content, 
    ml_options * options, 
    int * i, 
    int argc, 
    char ** argv);

int parse_ml_arg(
    char ** argv, 
    int flag_index, 
    int content_s, 
    int content_e, 
    ml_options * ml_options);

int pass_vote_options(ml_options * master_ml_options, VOTE_GRP * vote){
  vote->is_revoting = master_ml_options->is_revoting;
  vote->revoting_weight_power = master_ml_options->revoting_weight_power;
  vote->is_revoting_all_quartet = master_ml_options->is_revoting_all_quartet;
  vote->weight_power = master_ml_options->weight_power;
  vote->is_all_quartet = master_ml_options->is_all_quartet;

  return 0;
}

/* This function takes a point to some ml_options, popularizes its fields and 
 *    set the corresponding values depending on the user input
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

  // Initializing the options to default fields (which assumes that nothing has 
  // been computed before and this is the very first run)
  FCAL(
      GENERAL_ERROR,
      F_INIT_ML_OPT,
      init_ml_options(ml_options)
  ); 
  ASSERT(
      GENERAL_ERROR,
      WRONG_ARG_FORMAT,
      argv[1][0] == '-'
  );

  // State machine to read in user input
  is_flag = 1; // always start with a flag
  flag_index = 1;
  for(i = 2; i < argc; i++){
    if(is_flag){  // was reading a flag from the previous round, actual flag is 
                  // stored in flag_index
      ASSERT(
          GENERAL_ERROR,
          WRONG_ARG_FORMAT,
          argv[1][0] == '-'
      );
      

      // reading the first content, note that there may be a lot more contents 
      // to read
      content_s = i;
      is_flag = 0;

    } else { // have been reading content
      if(argv[i][0] == '-'){ // read a flag here, so no more content to read
        content_e = i;
        is_flag = 1;

        // Read the previous flag into option
        FCAL(
            GENERAL_ERROR,
            F_PARSE_ML_ARG,
            parse_ml_arg(argv, flag_index, content_s, content_e, ml_options)
        );
      
        // Reset the machine
        flag_index = i;
        content_s = -1;
        content_e = -1;
      } 
      // Else, reading another content, so we don't need to do anything
    }
  }

  // Read last flag
  FCAL(
      GENERAL_ERROR,
      F_PARSE_ML_ARG,
      parse_ml_arg(argv, flag_index, content_s, argc, ml_options)
  );

  return 0; 
}

/* Wrapper function to check for command length and loop through command to 
 *  parse options
 * Input    argc, argv from command line and an option struct pointer
 * Output   0 on success, ERROR otherwise
 * Effect   set fields in option struct
 */
int read_cmd_arg(int argc, char ** argv, ml_options * ml_options){
  int i;

  // Loop through command
  for(i = 1; i < argc; i++)
    if(argv[i][0] == '-') //currently reading a flag
      FCAL(
          GENERAL_ERROR, 
          F_FIND_ARG_INC_OPT, 
          find_arg_index(argv[i], argv[i + 1], ml_options, &i, argc, argv)
      );
  return 0;
}


int parse_ml_arg(char ** argv, 
                  int flag_index, 
                  int content_s, 
                  int content_e, 
                  ml_options * ml_options)
{
  int i;
  char * flag = argv[flag_index];
  // Switch case for the flag index, which should be a single character, 
  // following a dot. 
  // First check that there is exactly one character that follows the dash 
  ASSERT(
      GENERAL_ERROR,
      WRONG_ARG_FORMAT,
      strlen(flag) == 2 && 
          flag[0] == '-' && 
          flag[1] >= 'A' && 
          flag[1] <= 'z' &&
          content_e - content_s == 1
  );


  // Actually look up what exact field show we be setting and set said field
  // -a: input alignment              (expect the full path - does not check)
  // -i: initial tree method          (expect "fasttree/raxml/nj/me")
  // -o: output prefix                (expect the full path - does not check)
  // -t: initial tree                 (expect the full path - does not check)
  // -g: guide tree                   (expect the full path - does not check)
  // -m: initial distance matrix      (expect the full path - does not check)
  // -r: recompute constraint trees   (expect 0/1) 
  // -c: constraint trees type        (expect "no"/"subtree"/"fasttree"/
  //                                      "raxml"/"nj"/"me")
  // -q: quartet tree type            (expect "fpm"/"subtree"/"raxml"/"ml")
  // -d: distance model               (expect "JC"/"logDet"/"K2P"/"P")
  // -s: max subset size              (expect a positive number)
  // -p: tmp folder location          (expect a string)
  switch(flag[1]){   
    case 'p':
        ml_options->tmp_folder = argv[content_s];
        break;
    case 'a':
        ml_options->input_alignment = argv[content_s];
      break;
    case 'o':
        ml_options->output_prefix = argv[content_s];
      break;
    case 't':
        ml_options->init_tree_name = malloc(GENERAL_BUFFER_SIZE * sizeof(char));
        strcpy(ml_options->init_tree_name, argv[content_s]);
      break;
    case 'g':
        ml_options->guide_tree_name = malloc(GENERAL_BUFFER_SIZE * sizeof(char));
        strcpy(ml_options->guide_tree_name, argv[content_s]);
      break;    
    case 'm':
        ml_options->init_d_name = argv[content_s];
      break;
    case 'r': 
        ml_options->recompute_constraint_trees = argv[content_s][0] - '0';
      break;
    case 'c':
      if(content_e - content_s == 1)
        for(i = 0; i < N_CTREE_METHOD; i++)
          if(STR_EQ(argv[content_s], CTREE_FLAG[i]))
            ml_options->ctree_method = i;
      break;
    case 'i':
      if(content_e - content_s == 1)
        for(i = 0; i < N_ITREE_METHOD; i++)
          if(STR_EQ(argv[content_s], ITREE_FLAG[i]))
            ml_options->itree_method = i;
      break;
    case 'q':
      if(content_e - content_s == 1)
        for(i = 0; i < N_QTREE_METHOD; i++)
          if(STR_EQ(argv[content_s], QTREE_FLAG[i]))
            ml_options->qtree_method = i; 
      break;
    case 'd':
      if(content_e - content_s == 1)
        for(i = 0; i < N_DIST_MOD; i++)
          if(STR_EQ(argv[content_s], DIST_MOD_FLAG[i]))
            ml_options->distance_model = i;
      break; 
    case 's':
      ml_options->ss_threshold = atoi(argv[content_s]);
      break;
    case 'v':
        switch(argv[content_s][0]){
          case '1':
            ml_options->is_revoting = 0; 
            ml_options->weight_power = 0; 
            ml_options->is_all_quartet = 0; 
            break;

          case '2':
            ml_options->is_revoting = 1;
            ml_options->weight_power = 0; 
            ml_options->is_all_quartet = 0; 

            ml_options->revoting_weight_power = argv[content_s][2] == '1' ? 1 : 2;
            ml_options->is_revoting_all_quartet = 0;
            break;

          case '3':
            ml_options->is_revoting = 0; 
            ml_options->weight_power = 2; 
            ml_options->is_all_quartet = 0;
            break;

          case '4':
            ml_options->is_revoting = 1;
            ml_options->weight_power = 2; 
            ml_options->is_all_quartet = 0; 

            ml_options->revoting_weight_power = 2;
            ml_options->is_revoting_all_quartet = 1;
            break;

          case '5':
            ml_options->is_revoting = 0; 
            ml_options->weight_power = 2; 
            ml_options->is_all_quartet = 1;
            break; 
        }
      break;
  } 
  return 0;
}

/* Function to check the flag tag and assign its content to the appropriate 
 *  field in option structure
 * This function also stores the index in argv to the content
 * Input    the flag, the content and pointer to an option struct as well as
 *   index to argv
 * Output   0 on success, ERROR otherwise
 * Effect   set fields in option struct
 */
int find_arg_index(
    char * flag, 
    char * content, 
    ml_options * options, 
    int * i, 
    int argc, 
    char ** argv)                
{
  int j;

  ASSERT(
      GENERAL_ERROR,
      WRONG_ARG_FORMAT,
      flag[0] == '-'
  );

  switch(flag[1]){
    case 'i':
      options->init_d_name = content;  
      break;

    case 'o':
      options->output_prefix = content;
      break;

    case 't':
      options->num_trees = 0;

      while(1){
        if((*i) + options->num_trees + 1 >= argc || 
            argv[(*i) + options->num_trees + 1][0] == '-') 
          break;
        else options->num_trees++;
      }

      options->tree_names = malloc(options->num_trees * sizeof(char*));
      ASSERT(
          MALLOC_ERROR,
          M_ERR_IN_FIND_ARG_IDX,
          options->tree_names
      );

      for(j = 0; j < options->num_trees; j++){
        options->tree_names[j] = malloc(GENERAL_BUFFER_SIZE * sizeof(char));
        ASSERT(
            MALLOC_ERROR,
            M_ERR_IN_FIND_ARG_IDX,
            options->tree_names[j]
        );
        strcpy(options->tree_names[j], argv[(*i) + j + 1]);
      }
      
      (*i) += options->num_trees;
      break;

    default:
      ASSERT(GENERAL_ERROR, WRONG_ARG_FORMAT, 0);
      break;
  }

  return 0;
}


int init_ml_options(ml_options * ml_options){
  // The first two should be initialized properly
  ml_options->input_alignment         = NULL;
  ml_options->output_prefix           = NULL;
  ml_options->init_tree_name          = NULL;
  ml_options->init_d_name             = NULL;
  ml_options->guide_tree_name         = NULL; 

  ml_options->num_trees               = 0;
  ml_options->tree_names              = NULL;

  ml_options->is_revoting             = 0;
  ml_options->revoting_weight_power   = 0;
  ml_options->is_revoting_all_quartet = 0;
  ml_options->weight_power            = 0; 
  ml_options->is_all_quartet          = 0;


  ml_options->use_distance_matrix                 = 1;

  ml_options->recompute_constraint_trees          = 1;
  ml_options->itree_method                        = I_FASTTREE;
  ml_options->ctree_method                        = C_FASTTREE;
  ml_options->qtree_method                        = Q_FPM;
  ml_options->distance_model                      = D_LOGDET;

  ml_options->ss_threshold                            = -1;

  ml_options->use_initial_tree_as_spanning_tree       = 0;

  return 0;
}
