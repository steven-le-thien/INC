#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "utilities.h"
#include "msa.h"
#include "inc_ml.h"

// Internal functions
int setup_name(char ** name_holder, char* name);
int setup_sequence(char ** sequence_holder, char* sequence);
char * find_sequence_by_name(msa_t * msa, char * name);
int init_msa(msa_t* msa, int N, int num_seq, char** msa_core, char** msa_name);
int write_sub_msa(char * infile, char * taxon_name, msa_t * msa);
int rm_empty_col(msa_t * msa, msa_t * sub_msa, char * outfile);

int subset_msa(char * infile, char * outfile, msa_t * msa){
  char taxon_name[MAX_NAME_SIZE];

  msa_t sub_msa;

  sprintf(taxon_name, "%s_tmp", outfile);
  FCAL(
      GENERAL_ERROR,
      F_WRITE_SUB_MSA_IN_SUBSET_MSA,
      write_sub_msa(infile, taxon_name, msa)
  );

  FCAL(
      GENERAL_ERROR,
      F_PARSE_INPUT_IN_SUBSET_MSA,
      parse_input(&sub_msa, taxon_name)
  );

  FCAL(
      GENERAL_ERROR,
      F_RM_EMPTY_COL_IN_SUBSET_MSA,
      rm_empty_col(msa, &sub_msa, outfile)
  );

  SYSCAL(
      GENERAL_ERROR,
      RM_ERR,
      "rm %s",
      taxon_name
  );
  return 0;
}

/* Parse FASTA file into msa_t struct. This function assumes a maximum length 
 *    of a sequence and number of sequences currently
 * Input:   pointer to the msa and name of the FASTA file
 * Output:  0 on sucess, GENERAL_ERROR otherwise
 * Effect:  may print onto outstream if error occurs, set fields in the msa 
 *              struct, open instream, assuming filename is in the same 
 *              directory as the binary file calls malloc
 */
int parse_input(msa_t * msa_ptr, char * filename){
  // Placeholder for contents of the msa
  char** msa_core;
  char** msa_name;

  // Input FASTA file
  FILE * f;

  // File reader variables
  int sequence_counter;
  char line[MAX_SEQUENCE_LENGTH];
  char seq[MAX_SEQUENCE_LENGTH];

  // Safe checking
  ASSERT(GENERAL_ERROR, N_MSA_PTR_IN_PARSE_IN, msa_ptr);
  ASSERT(GENERAL_ERROR, N_FILENAME_IN_PARSE_IN, filename);

  // Allocate space
  msa_core = SAFE_MALLOC(MAX_NUM_SEQUENCE * sizeof(char*));
  msa_name = SAFE_MALLOC(MAX_NUM_SEQUENCE * sizeof(char*)); 

  // Clear strings    
  sequence_counter = 0;
  STR_CLR(seq);
  STR_CLR(line);

  // Redirecting stdin
  f = SAFE_FREOPEN_RD(filename, stdin);

  // Ignoring comments
  while(1){
    ASSERT(GENERAL_ERROR, NO_SEQUENCE, scanf("%s", line) >= 0);
    if(str_start_with(line, '>')) break; 
  }

  // Set up name for the first sequence
  FCAL(
      GENERAL_ERROR,
      F_SETUP_NAME_IN_PARSE_INPUT,
      setup_name(&msa_name[sequence_counter], &line[1])
  );

  // Read the stdin line by line until EOF signal
  while(scanf("%s", line) >= 0){
    switch(*line){ //read the first character of the word
      case ';': break;
      case '>': // Finished previous seuqence
        FCAL(
            GENERAL_ERROR,
            F_SETUP_SEQ_IN_PARSE_INPUT,
            setup_sequence(&msa_core[sequence_counter], seq)
        );
        sequence_counter++;
        FCAL(
            GENERAL_ERROR,
            F_SETUP_NAME_IN_PARSE_INPUT,
            setup_name(&msa_name[sequence_counter], &line[1])
        );
        STR_CLR(seq);
        break;
      default:
        strcat(seq, line);
    }
  }

  // Final iteration for the last sequence
  FCAL(
      GENERAL_ERROR,
      F_SETUP_SEQ_IN_PARSE_INPUT,
      setup_sequence(&msa_core[sequence_counter], seq)
  );
  sequence_counter++;

  // Close file
  fclose(f);

  // Fill in the fields
  return(init_msa(msa_ptr, strlen(seq), sequence_counter, msa_core, msa_name));
}

// Very naive algorithm (time can at least be halvede easily)
int compute_k2p_distance(msa_t * msa, float *** d){
  int i, j, k;
  int p, q;
  char * s1;
  char * s2;

  (*d) = malloc(msa->num_seq * sizeof(float *));

  for(i = 0; i < msa->num_seq; i++)
    (*d)[i] = malloc(msa->num_seq * sizeof(float));

  for(i = 0; i < msa->num_seq; i++){
    for(j = 0; j < msa->num_seq; j++){
      p = 0;
      q = 0;
      s1 = msa->msa[i];
      s2 = msa->msa[j];
      for(k = 0; k < msa->N; k++){
        p += (s1[k] == 'A' && s2[k] == 'G') || 
            (s1[k] == 'G' && s2[k] == 'A') ||
            (s1[k] == 'T' && s2[k] == 'C') ||
            (s1[k] == 'C' && s2[k] == 'T');

        q += (s1[k] == 'A' && s2[k] == 'T') || 
            (s1[k] == 'T' && s2[k] == 'A') ||
            (s1[k] == 'G' && s2[k] == 'C') ||
            (s1[k] == 'C' && s2[k] == 'G') ||

            (s1[k] == 'A' && s2[k] == 'C') || 
            (s1[k] == 'C' && s2[k] == 'A') ||
            (s1[k] == 'G' && s2[k] == 'T') ||
            (s1[k] == 'T' && s2[k] == 'G');
      }

      (*d)[i][j] = -0.5 * log2((1.0 - 2.0 * p / msa->N - 1.0 * q / msa->N) *
          sqrt(1 - 2 * q / msa->N));
    }
  }

  return 0;
}

int write_distance_matrix(float ** d, ml_options * options, msa_t * msa){
  FILE * f;
  int i, j;
  char filename[MAX_BUFFER_SIZE];

  STR_CLR(filename);
  strcat(filename, options->output_prefix);
  strcat(filename, DEFAULT_DIST_SUF);

  f = fopen(filename, "r");
  if(f){  // such file already exists, then we just skip and assume 
          // that we have the correct distance matrix
    fclose(f);
    return 0;
  }

  f = fopen(filename, "w");
  fprintf(f, "%d\n", msa->num_seq);
  for(i = 0; i < msa->num_seq; i++){
    fprintf(f, "%s\n", msa->name[i]);
    for(j = 0; j < msa->num_seq; j++)
      fprintf(f, "%f ", d[i][j]);
    
    fprintf(f, "\n");
  }
  fclose(f);
  return 0;
}

// INTERNAL IMPLEMENTATION

/* Wrapper to malloc name holder and copy name into that place
 * Input:   the name and the location to be copied to
 * Output:  0 on success, ERROR otherwise   
 * Effect:  calls malloc
 */
int setup_name(char ** name_holder, char* name){
  ASSERT(GENERAL_ERROR, N_NAME_IN_SETUP_NAME, name);
  *name_holder = (char *) SAFE_MALLOC(strlen(name) + 1);
  strcpy(*name_holder, name);
  name_holder[0][strlen(name)] = 0;
  return 0;
}

/* Wrapper to malloc sequence holder and copy sequence into that place
 * Input:   the sequence and the location to be copied to
 * Output:  0 on success, ERROR otherwise   
 * Effect:  calls malloc
 */
int setup_sequence(char ** sequence_holder, char* sequence){
  // The syntax is exactly the same as setup_name so we will just call that 
  // function instead
  return setup_name(sequence_holder, sequence);
}


/* Brute force algorithm to find sequence based on name
 * Input:   msa structure and the name of the sequence
 * Output:  pointer to the sequence on success, NULL otherwise
 * Effect:  none
 */
char * find_sequence_by_name(msa_t * msa, char * name){
  int i; //loop variable

  ASSERT(NULL, N_MSA_IN_FIND_SEQ_BY_NAME, msa);
  ASSERT(NULL, N_NAME_IN_FIND_SEQ_BY_NAME, name);

  for(i = 0; i < msa->num_seq; i++)
    if(STR_EQ(name, msa->name[i]))
      return msa->msa[i];

  return NULL;
}

/* Initialize fields in the MSA to predetermined values
 * Input:   msa         pointer to the MSA srtuct
 *          N           length of each sequence
 *          num_seq     number of sequences in the MSA
 *          msa_core    array of sequences 
 *          msa_name    array of name for each sequence
 * Output:  0 on success, ERROR otherwise
 * Effect   set fields in the input msa
 */
int init_msa(msa_t* msa, int N, int num_seq, char** msa_core, char** msa_name){
  // Safety check
  ASSERT(GENERAL_ERROR, N_MSA_IN_INIT_MSA, msa);
  ASSERT(GENERAL_ERROR, N_MSA_CORE_IN_INIT_MSA, msa_core);
  ASSERT(GENERAL_ERROR, N_MSA_NAME_IN_INIT_MSA, msa_name);

  // Set the fields
  msa->N = N;
  msa->num_seq = num_seq;
  msa->msa = msa_core;
  msa->name = msa_name;
  return 0;
}

int write_sub_msa(char * infile, char * taxon_name, msa_t * msa){
  FILE *f, *p;
  int i;

  f = SAFE_FOPEN_RD(infile);
  p = fopen(taxon_name, "w"); 

  while(fscanf(f, "%s", taxon_name) >= 0)
    for(i = 0; i < msa->num_seq; i++)
      if(strcmp(taxon_name, msa->name[i]) == 0)
        fprintf(p, ">%s\n%s\n", taxon_name, msa->msa[i]);
  fclose(f);
  fclose(p);


  return 0;
}

int rm_empty_col(msa_t * msa, msa_t * sub_msa, char * outfile){
  int i, j;
  int * masked = SAFE_MALLOC(msa->N * sizeof(int)); 
  FILE * f = fopen(outfile, "w");
  for(i = 0; i < sub_msa->N; i++){
    masked[i] = 1;
    for(j = 0; j < sub_msa->num_seq; j++)
      masked[i] = masked[i] && sub_msa->msa[j][i] == '-';
  }

  for(i = 0; i < sub_msa->num_seq; i++){
    fprintf(f, ">%s\n", sub_msa->name[i]);
    for(j = 0; j < sub_msa->N; j++)
      if(!masked[j]) 
        fprintf(f, "%c", sub_msa->msa[i][j]);
    
    if(i != sub_msa->num_seq - 1) 
      fprintf(f, "\n");
  }
  fclose(f);
  return 0;
}