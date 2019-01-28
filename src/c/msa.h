#ifndef MSA_H
#define MSA_H

#include "c_inc.h"

// Constants
const static int MAX_SEQUENCE_LENGTH    = (int) 1e6;
const static int MAX_NUM_SEQUENCE       = (int) 1e6;

extern int parse_input(msa_t * msa_ptr, char * filename);
extern int compute_k2p_distance(msa_t * msa, float *** d);
extern int write_distance_matrix(float ** d, ml_options * options, msa_t * msa);
extern int subset_msa(char * infile, char * outfile, msa_t * msa);

const static char F_WRITE_SUB_MSA_IN_SUBSET_MSA[]
  = "write_sub_msa failed in subset_msa\n";

const static char F_PARSE_INPUT_IN_SUBSET_MSA[]
  = "parse_input failed in subset_msa\n";

const static char F_RM_EMPTY_COL_IN_SUBSET_MSA[]
  = "remove empty col failed in subset msa\n";

const static char RM_ERR[] = "remove error\n";

const static char N_MSA_PTR_IN_PARSE_IN[]
  = "null msa ptr in parse_input\n";

const static char N_FILENAME_IN_PARSE_IN[]
  = "null filanme in parse_input\n";

const static char NO_SEQUENCE[]
  = "no sequence in file\n";

const static char F_SETUP_NAME_IN_PARSE_INPUT[]
  = "set up name failed in parse_input\n";

const static char F_SETUP_SEQ_IN_PARSE_INPUT[]
  = "set up seq failed in parse_input\n";

const static char N_NAME_IN_SETUP_NAME[]
  = "null name in setup_name\n";

const static char N_MSA_IN_FIND_SEQ_BY_NAME[]
  = "null msa in find_seq_by_name\n";

const static char N_NAME_IN_FIND_SEQ_BY_NAME[]
  = "null name in find seq by name\n";

const static char N_MSA_IN_INIT_MSA[]
  = "null msa in init_msa\n";

const static char N_MSA_CORE_IN_INIT_MSA[]
  = "null msa_core in init_msa\n";

const static char N_MSA_NAME_IN_INIT_MSA[]
  = "null_msa_name in init_msa\n";

#endif