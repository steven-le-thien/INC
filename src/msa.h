#ifndef MSA_H
#define MSA_H

#include "c_inc.h"

// Constants
const static int MAX_SEQUENCE_LENGTH    = (int) 1e6;
const static int MAX_NUM_SEQUENCE       = (int) 1e6;

// Structure for the multiple sequence alignment
typedef struct msa {
    int     N;          // size of one sequence
    int     num_seq;    // number of sequences
    char**  msa;
    char**  name;
} msa_t;

extern int parse_input(msa_t * msa_ptr, char * filename);
extern int compute_k2p_distance(msa_t * msa, float *** d);
extern int write_distance_matrix(float ** d, option_t * options, msa_t * msa);
extern int subset_msa(char * infile, char * outfile, msa_t * msa);

#endif