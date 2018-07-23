#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dist.h"
#include "utilities.h"
#include "options.h"

#define is_fasta(string) (string[0] == '>')
#define is_phylip(string) ('0' <= string[0] && string[0] <= '9')

int parse_distance_matrix(char ** master_name_map, int ** adj_mat, option_t * options, int * n){
	FILE * f;
	char buf[1000050];
	char command[10000];
	strclr(buf);
	strclr(command);

	f = fopen(options->input_name, "r");
	if(!f) PRINT_AND_RETURN("cannot open input file in parse_distance_matrix", OPEN_ERROR);

	// Read the first word of the file to see whether it is a FASTA file or a PHYLIP distance matrix
	if(fscanf(f, "%s", buf) < 0) PRINT_AND_RETURN("empty input file", GENERAL_ERROR); 
	fclose(f);


	if(is_fasta(buf)){
		// Route through a PHYLIP creater (for now)
		default_fp_options.input_name = options->input_name;
		fastphylo_job(&default_fp_options);
		return read_phylip(master_name_map, adj_mat, options->input_name, n);
	} else if(is_phylip(buf)) {
		return read_phylip(master_name_map, adj_mat, options->input_name, n);
	} else {
		PRINT_AND_RETURN("invalid input format", GENERAL_ERROR);
	}
}

int read_phylip(char ** master_name_map, float ** adj_mat, char * filename, int * n){
	FILE * f;
	int num_sequence;
	int i, j;

	f = fopen(filename, "r");
	if(!f) PRINT_AND_RETURN("cannot open input file in parse_distance_matrix", OPEN_ERROR);
	fscanf(f, "%d", num_sequence);

	// Mallocation sequence
	master_name_map 	= malloc(num_sequence * sizeof(char*));
	adj_mat 			= malloc(num_sequence * sizeof(float *));

	for(i = 0; i < num_sequence; i++){
		master_name_map[i] = malloc(maxNameSize * sizeof(char));
		strclr(master_name_map[i]);

		adj_mat[i] = malloc(num_sequence * sizeof(char));
	}

	// Copy sequence (this could have been combined with the previous section for cache optimization but for now we separate them for clarity)
	for(i = 0; i < num_sequence; i++){
		fscanf(f, "%s", master_name_map[i]);

		for(j = 0; j < num_sequence; j++)
			fscanf(f, "%f", adj_mat[i][j]);
		
	}

	*n = num_sequence;
	return 0;
}