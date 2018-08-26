// File in inc_ml, created by Thien Le in July 2018


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dist.h"
#include "utilities.h"
#include "options.h"
#include "tools.h"

#define IS_FASTA(string)    (string[0] == '>')
#define IS_PHYLIP(string)   ('0' <= string[0] && string[0] <= '9')

// Internal functions
int read_phylip(INC_GRP * meta, MAP_GRP * map, char * filename);

/* Initialization matrix-related fields (number of taxon, etc)
 * Input:       meta        meta-variables including the distance matrix
 *              map         maps, including a naming map (since the taxa are digitized later)
 *              options     options set by users, including input file
 * Output:      0 on success, ERROR otherwise
 * Effect:      open a file, try to mallocate some structure
 */
int parse_distance_matrix(INC_GRP * meta, MAP_GRP * map, option_t * options){
    FILE * f;
    char buf[1000050];
    strclr(buf);

    f = fopen(options->input_name, "r");
    if(!f)                          PRINT_AND_RETURN("cannot open input file in parse_distance_matrix", OPEN_ERROR);

    // Read the first word of the file to see whether it is a FASTA file or a PHYLIP distance matrix
    if(fscanf(f, "%s", buf) < 0)    PRINT_AND_RETURN("empty input file", GENERAL_ERROR); 
    fclose(f);

    if(IS_FASTA(buf)){ // Route through a PHYLIP creater (for now) THIS FEATURE IS NOT TESTED
        // default_fp_options.input_name = options->input_name;
        // fastphylo_job(&default_fp_options);

        //                                                                                     #if DEBUG
        //                                                                                         printf("debug: input file is not a  distance matrix\n"); 
        //                                                                                     #endif

        // return read_phylip(meta, map, default_fp_options.output_name);
        return -1;

    } else if(IS_PHYLIP(buf)) {
                                                                                            #if DEBUG
                                                                                                printf("debug: input file is a distance matrix\n"); 
                                                                                            #endif
        return read_phylip(meta, map, options->input_name);

    } else {

        PRINT_AND_RETURN("invalid input format", GENERAL_ERROR);
    }
}

// INTERNAL IMPLEMENTATION
/* Reading a phylip distance matrix files and initialize meta variables
 * Input:       meta        meta-variables including the distance matrix
 *              map         maps, including a naming map (since the taxa are digitized later)
 *              filename    name of the phylip distance matrix file
 * Output:      0 on success, ERROR otherwise
 * Effect:      open a file, try to mallocate some structure, modify meta and map
 */
int read_phylip(INC_GRP * meta, MAP_GRP * map, char * filename){
    FILE * f;   // file object
    int n = 0;      // number of sequence
    int i, j;   // loop counters
    float min;

    printf("%s\n", filename);
    // Open file
    f = fopen(filename, "r");
    if(!f)                          PRINT_AND_RETURN("cannot open input file in parse_distance_matrix", OPEN_ERROR);

    // Record number of sequence
    if(fscanf(f, "%d", &n) < 0)         PRINT_AND_RETURN("input phylip file is empty", GENERAL_ERROR);
    meta->n_taxa    = n;
    map->n_taxa     = n; 
                                                                                            #if DEBUG 
                                                                                                printf("debug: number of taxa is %d %d (these 2 should be the same)\n", meta->n_taxa, map->n_taxa); 
                                                                                            #endif
    // Mallocation sequence
    map->master_to_name         = malloc(n * sizeof(char*));
    meta->d                     = malloc(n * sizeof(float *));

    if(!map->master_to_name)        PRINT_AND_RETURN("malloc failed to allocate name map in read_phylip", MALLOC_ERROR);
    if(!meta->d)                    PRINT_AND_RETURN("malloc failed to allocate distance matrix in read_phylip", MALLOC_ERROR);


    for(i = 0; i < n; i++){
        // Allocation
        map->master_to_name[i]  = malloc(MAX_NAME_SIZE * sizeof(char));
        meta->d[i]              = malloc(n * sizeof(float));

        if(!map->master_to_name[i]) PRINT_AND_RETURN("malloc failed to allocate name map in read_phylip", MALLOC_ERROR);
        if(!meta->d[i])             PRINT_AND_RETURN("malloc failed to allocate distance matrix in read_phylip", MALLOC_ERROR);

        // Initialization
        fscanf(f, "%s", map->master_to_name[i]);
        for(j = 0; j < n; j++){
            fscanf(f, "%f", &(meta->d[i][j]));
        }
    }
    fclose(f);

    // Corrected 2kp (for 0 entries)
    min = 1.0 * INT_MAX; 
    for(i = 0; i < meta->n_taxa; i++)
        for(j = 0; j < meta->n_taxa; j++)
            if(min - meta->d[i][j] > EPS && meta->d[i][j] > EPS)
                min = meta->d[i][j];

    min /= 2;
    for(i = 0; i < meta->n_taxa; i++)
        for(j = 0; j < meta->n_taxa; j++)
            if(meta->d[i][j] < EPS && i != j)
                meta->d[i][j] = min;


                                                                                                #if LARGE_DEBUG 
                                                                                                    printf("debug: this prints out the list of name\n");
                                                                                                    for(i = 0; i < meta->n_taxa; i++)
                                                                                                        printf("%d:%s ", i, map->master_to_name[i]); 
                                                                                                    printf("\n");
                                                                                                #endif 
                                                                                                #if DEBUG
                                                                                                    printf("debug: this tests whether all the distances are positive\n");
                                                                                                    int k;
                                                                                                    k = 0;
                                                                                                    for(i = 0; i < meta->n_taxa; i++)
                                                                                                        for(j = 0; j < meta->n_taxa; j++)
                                                                                                            if(meta->d[i][j] > 1e-10)
                                                                                                                k++;
                                                                                                            else if(i != j) printf("debug: problem is %f, at %d %d, name is %s %s\n", meta->d[i][j], i, j, map->master_to_name[i], map->master_to_name[j]);
                                                                                                    printf("debug: double check, num sequence is %d\n", meta->n_taxa);
                                                                                                    printf("debug: there are %d positive entries\n", k);
                                                                                                    // for(i = 0; i < meta->n_taxa; i++){
                                                                                                    //     for(j = 0; j < meta->n_taxa; j++){
                                                                                                    //         printf("%f ", meta->d[i][j]);
                                                                                                    //     }
                                                                                                    //     printf("\n");
                                                                                                    // }
                                                                                                    // while(1);
                                                                                                    // printf("1 45 %f 1 57 %f 1 27 %f 45 57 %f 45 27 %f 57 27 %f\n", meta->d[0][44], meta->d[0][56], meta->d[0][26], meta->d[44][56], meta->d[44][26],meta->d[56][26]);
                                                                                                #endif
    return 0;
}
