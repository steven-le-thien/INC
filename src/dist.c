#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dist.h"
#include "utilities.h"
#include "options.h"

#define IS_FASTA(string) 	(string[0] == '>')
#define IS_PHYLIP(string) 	('0' <= string[0] && string[0] <= '9')

// Internal functions
int read_phylip(INC_GRP * meta, MAP_GRP * map, option_t * options);

/* Initialization matrix-related fields (number of taxon, etc)
 * Input: 		meta 		meta-variables including the distance matrix
 *				map 		maps, including a naming map (since the taxa are digitized later)
 * 				options 	options set by users, including input file
 * Output: 		0 on success, ERROR otherwise
 * Effect: 		open a file, try to mallocate some structure
 */
int parse_distance_matrix(INC_GRP * meta, MAP_GRP * map, option_t * options){
	FILE * f;
	char buf[1000050];
	char command[10000];
	strclr(buf);
	strclr(command);

	f = fopen(options->input_name, "r");
	if(!f) 							PRINT_AND_RETURN("cannot open input file in parse_distance_matrix", OPEN_ERROR);

	// Read the first word of the file to see whether it is a FASTA file or a PHYLIP distance matrix
	if(fscanf(f, "%s", buf) < 0) 	PRINT_AND_RETURN("empty input file", GENERAL_ERROR); 
	fclose(f);

	if(IS_FASTA(buf)){ // Route through a PHYLIP creater (for now)
		default_fp_options.input_name = options->input_name;
		fastphylo_job(&default_fp_options);

		return read_phylip(meta, map, default_fp_options->output_name);

	} else if(IS_PHYLIP(buf)) {

		return read_phylip(meta, map, options->input_name);

	} else {

		PRINT_AND_RETURN("invalid input format", GENERAL_ERROR);
	}
}

// INTERNAL IMPLEMENTATION
/* Reading a phylip distance matrix files and initialize meta variables
 * Input: 		meta 		meta-variables including the distance matrix
 *				map 		maps, including a naming map (since the taxa are digitized later)
 * 				filename 	name of the phylip distance matrix file
 * Output: 		0 on success, ERROR otherwise
 * Effect: 		open a file, try to mallocate some structure, modify meta and map
 */
int read_phylip(INC_GRP * meta, MAP_GRP * map, char * filename){
	FILE * f; 	// file object
	int n; 		// number of sequence
	int i, j;	// loop counters

	// Open file
	f = fopen(filename, "r");
	if(!f) 							PRINT_AND_RETURN("cannot open input file in parse_distance_matrix", OPEN_ERROR);

	// Record number of sequence
	if(fscanf(f, "%d", n) < 0) 		PRINT_AND_RETURN("input phylip file is empty", GENERAL_ERROR);
	meta->n_taxa 	= n;
	map->n_taxa 	= n; 

	// Mallocation sequence
	map->master_to_name			= malloc(n * sizeof(char*));
	meta->d 					= malloc(n * sizeof(float *));

	if(!map->master_to_name) 		PRINT_AND_RETURN("malloc failed to allocate name map in read_phylip", MALLOC_ERROR);
	if(!meta->d)					PRINT_AND_RETURN("malloc failed to allocate distance matrix in read_phylip", MALLOC_ERROR);


	for(i = 0; i < num_sequence; i++){
		// Allocation
		map->master_to_name[i] 	= malloc(MAX_NAME_SIZE * sizeof(char));
		meta->d[i] 				= malloc(n * sizeof(char));

		if(!map->master_to_name[i]) PRINT_AND_RETURN("malloc failed to allocate name map in read_phylip", MALLOC_ERROR);
		if(!meta->d[i]) 			PRINT_AND_RETURN("malloc failed to allocate distance matrix in read_phylip", MALLOC_ERROR);

		// Initialization
		fscanf(f, "%s", map->master_to_name[i]);
		for(j = 0; j < num_sequence; j++)
			fscanf(f, "%f", &(meta->d[i][j]));

	}
	return 0;
}