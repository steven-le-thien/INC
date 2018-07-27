#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utilities.h"
#include "options.h"
#include "tools.h"

// Constants
const static int MAX_SEQUENCE_LENGTH    = (int) 1e6;
const static int MAX_NUM_SEQUENCE       = (int) 1e6;
const static int MAX_NAME_LENGTH        = (int) 1e3;

// Structure for the multiple sequence alignment
typedef struct msa {
    int     N;          // size of one sequence
    int     num_seq;    // number of sequences
    char**  msa;
    char**  name;
} msa_t;

/* Wrapper to malloc name holder and copy name into that place
 * Input:   the name and the location to be copied to
 * Output:  0 on success, ERROR otherwise   
 * Effect:  calls malloc
 */
int setup_name(char ** name_holder, char* name){
    if(!name)   PRINT_AND_RETURN("name is NULL when trying to parse",   GENERAL_ERROR);

    *name_holder = (char *) malloc(strlen(name));
    if(!*name_holder) PRINT_AND_RETURN("malloc failed for name holder", MALLOC_ERROR);
    strcpy(*name_holder, name);
    return 0;
}

/* Wrapper to malloc sequence holder and copy sequence into that place
 * Input:   the sequence and the location to be copied to
 * Output:  0 on success, ERROR otherwise   
 * Effect:  calls malloc
 */
int setup_sequence(char ** sequence_holder, char* sequence){
    // The syntax is exactly the same as setup_name so we will just call that function instead
    return setup_name(sequence_holder, sequence);
}


/* Brute force algorithm to find sequence based on name
 * Input:   msa structure and the name of the sequence
 * Output:  pointer to the sequence on success, NULL otherwise
 * Effect:  none
 */
char * find_sequence_by_name(msa_t * msa, char * name){
    int i; //loop variable

    if(!msa)        PRINT_AND_RETURN("msa is NULL in find_sequence_by_name",    NULL);
    if(!name)       PRINT_AND_RETURN("name is NULL in find_sequence_by_name",   NULL);
    for(i = 0; i < msa->num_seq; i++)
        if(strcmp(name, msa->name[i]) == 0)
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
    if(!msa)        PRINT_AND_RETURN("msa is NULL from init_msa",       GENERAL_ERROR);
    if(!msa_core)   PRINT_AND_RETURN("msa_core is NULL from init_msa",  GENERAL_ERROR);
    if(!msa_name)   PRINT_AND_RETURN("msa_name is NULL from init_msa",  GENERAL_ERROR);

    // Set the fields
    msa->N = N;
    msa->num_seq = num_seq;
    msa->msa = msa_core;
    msa->name = msa_name;
    return 0;
}

/* Parse FASTA file into msa_t struct. This function assumes a maximum length of a sequence and number of sequences currently
 * Input:   pointer to the msa and name of the FASTA file
 * Output:  0 on sucess, GENERAL_ERROR otherwise
 * Effect:  may print onto outstream if error occurs, set fields in the msa struct, 
 *          open instream, assuming filename is in the same directory as the binary file
 *          calls malloc
 */
int parse_input(msa_t * msa_ptr, char * filename){
    // Placeholder for contents of the msa
    char** msa_core;
    char** msa_name;

    int i; //loop variable

    // Input FASTA file
    FILE * f;

    // File reader variables
    int sequence_counter;
    char line[MAX_SEQUENCE_LENGTH];
    char seq[MAX_SEQUENCE_LENGTH];

    // Safe checking
    if(!msa_ptr)    PRINT_AND_RETURN("msa_ptr is NULL in parse input",          GENERAL_ERROR);
    if(!filename)   PRINT_AND_RETURN("filename is NULL in parse input",         GENERAL_ERROR);

    // Allocate space
    msa_core = malloc(MAX_NUM_SEQUENCE * sizeof(char*));
    msa_name = malloc(MAX_NUM_SEQUENCE * sizeof(char*)); 
    if(!msa_core)   PRINT_AND_RETURN("malloc for msa_core failed in parse_input",   MALLOC_ERROR);
    if(!msa_name){
        free(msa_core);
        PRINT_AND_RETURN("malloc for msa_name failed in parse_input",   MALLOC_ERROR);
    }

    // Clear strings    
    sequence_counter = 0;
    strclr(seq);
    strclr(line);

    // Redirecting stdin
    f = freopen(filename, "r", stdin);
    if(!f) PRINT_AND_RETURN("fail to open input file",  OPEN_ERROR);    

    // Ignoring comments
    while(1){
        if(scanf("%s", line) < 0) PRINT_AND_RETURN("input file contains no sequence",   GENERAL_ERROR);
        if(str_start_with(line, '>')) break; 
    }

    // Set up name for the first sequence
    if(setup_name(&msa_name[sequence_counter], &line[1]) != SUCCESS){ //roll back and deallocate
        free(msa_core);
        free(msa_name);
    }

    // Read the stdin line by line until EOF signal
    while(scanf("%s", line) >= 0){
        switch(*line){ //read the first character of the word
            case ';': break;
            case '>': // Finished previous seuqence
                if(setup_sequence(&msa_core[sequence_counter], seq) != SUCCESS){ // roll back and deallocate
                    for(i = 0; i < sequence_counter; i++){
                        if(i != sequence_counter - 1) free(msa_core[i]);
                        free(msa_name[i]);
                    }
                    free(msa_core);
                    free(msa_name);
                    PRINT_AND_RETURN("sequence allocation faiiled in no sequence", MALLOC_ERROR);
                }

                sequence_counter++;
                if(setup_name(&msa_name[sequence_counter], &line[1]) != SUCCESS){ // roll back and deallocate
                    for(i = 0; i < sequence_counter; i++){
                        free(msa_core[i]);
                        free(msa_name[i]);
                    }
                    free(msa_core);
                    free(msa_name);
                    PRINT_AND_RETURN("name allocation faiiled in no sequence", MALLOC_ERROR);
                }
                strclr(seq);
                break;
            default:
                strcat(seq, line);
        }
    }

    // Final iteration for the last sequence
    if(setup_sequence(&msa_core[sequence_counter], seq) != SUCCESS){ // roll back and deallocate
        for(i = 0; i < sequence_counter; i++){
            if(i != sequence_counter - 1) free(msa_core[i]);
            free(msa_name[i]);
        }
        free(msa_core);
        free(msa_name);
        PRINT_AND_RETURN("sequence allocation faiiled in no sequence", MALLOC_ERROR);
    }
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
		// printf("i is %d\n", i);
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
						(s1[k] == 'C' && s2[k] == 'G');
			}

			(*d)[i][j] = -0.5 * log2((1.0 - 2.0 * p / msa->N - 1.0 * q / msa->N) * sqrt(1 - 2 * q / msa->N));
		}
	}

	return 0;
}

int write_distance_matrix(float ** d, char * filename, msa_t * msa){
	FILE * f;
	int i, j;

	f = fopen(filename, "w");
	fprintf(f, "%d\n", msa->num_seq);
	for(i = 0; i < msa->num_seq; i++){
		// printf("i is %d\n", i);
		fprintf(f, "%s\n", msa->name[i]);
		for(j = 0; j < msa->num_seq; j++){
					// printf("j is %d\n", j);

			fprintf(f, "%f ", d[i][j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
	return 0;
}


int main(int argc, char ** argv){
	// Pipe in a bunch of programs together
	option_t    options;
	msa_t 		msa;
	FILE * 		f;
	float ** 	d;
	int 		num_ctree;
	char ** 	ctree_name;
	char * 		input_name;
	char *		output_dir;
	char 		name[CMD_BUFFER_SIZE];
	int i;
	char 		num[100];

	if(init_options(&options)               != SUCCESS)         PRINT_AND_EXIT("init_options failed in main\n", GENERAL_ERROR);
    if(read_cmd_arg(argc, argv, &options)   != SUCCESS)         PRINT_AND_EXIT("read_cmd_arg failed in main\n", GENERAL_ERROR);

 //    // Piping into fasttree 
    if(fasttree_job(&options) 				!= SUCCESS) 		PRINT_AND_EXIT("fasttree_job failed in main\n", GENERAL_ERROR);

 //    // Piping into modified UPP code
    if(upp_job(&options)					!= SUCCESS) 		PRINT_AND_EXIT("upp job failed in main\n", GENERAL_ERROR);

 //    // Count the number of c trees
    num_ctree = 0;
    while(1){
    	strclr(num);
    	sprintf(num, "%d", num_ctree);

    	strclr(name);
    	strcat(name, "small_ctree");
    	strcat(name, num);
    	strcat(name, ".tree");
    	f = fopen(name, "r");
    	if(!f) break;
    	else num_ctree++;
    }

    printf("parsing input.. \n");
    parse_input(&msa, options.input_name);

    printf("finish computing distance ...\n");
    compute_k2p_distance(&msa, &d);

    printf("finish writing matrix ...\n");
    write_distance_matrix(d, "c_inc_input", &msa);

    // Piping into constrained_inc code
    constraint_inc(num_ctree);

    return 0; 
}	

