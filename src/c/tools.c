// File in inc_ml, created by Thien Le in July 2018

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>

#include "tools.h"
#include "utilities.h"

#if 1
char RAxML_bin[]    = "raxmlHPC-AVX2";
char FastTree_bin[] = "FastTree"; 
char PAUP_bin[]     = "paup4a163_osx";
#else
char RAxML_bin[]    = "raxmlHPC-PTHREADS-AVX2";
char FastTree_bin[] = "FastTree";
char PAUP_bin[]     = "paup4a163_centos64";
#endif

// Private functions
void add_command(char * current_command, char * new_command);
int find_prefix_and_dir_from_path(char * path, char * prefix, char * dir);

int construct_unweighted_matrix_job(char * filename, char * output_prefix, float ** dm, char ** name_map, int * master_to_midx){
    option_t tmp_options;
    FILE * f;
    char output_name[GENERAL_BUFFER_SIZE];
    char name[GENERAL_BUFFER_SIZE];
    int i, j, n;

    sprintf(output_name, "%ssecondary_matrix", output_prefix);
    f = fopen(output_name, "r");

    if(!f){
        tmp_options.input_name = filename; 
        tmp_options.output_name = output_name;
        // Call the actual maker
        if(unweighted_job(&tmp_options) != SUCCESS)             PRINT_AND_RETURN("construct_unweighted_matrix_job failed\n", GENERAL_ERROR);

        f = fopen(output_name, "r");
    }
    
    fscanf(f, "%d", &n);
    for(i = 0; i < n; i++){
        // Locate name and update the mapping to the master indexing scheme
        fscanf(f, "%s", name);
        for(j = 0; j < n; j++){
            if(strcmp(name, name_map[j]) == 0){
                master_to_midx[j] = i;
                break;
            }
        }

        // Copy the actual distance matrix
        for(j = 0; j < n; j++){
            fscanf(f, "%f", &dm[i][j]);
        }
    }
    fclose(f);

    return 0;
}

int make_subset_label(char * tree_name, char * out_name, ml_options * master_ml_options){
    option_t tmp_options; 

    tmp_options.input_name = tree_name;
    tmp_options.output_name = out_name;
    if(subset_job(&tmp_options, master_ml_options)                 != SUCCESS)         PRINT_AND_RETURN("subset job failed in main\n", GENERAL_ERROR);

    return 0;
}

int make_subtree(char * label, char * outname, char * tree_name){
    option_t tmp_options; 
    tmp_options.input_name = tree_name;
    tmp_options.output_name = tree_name; // replacing file
    if(rm_label_job(&tmp_options) != SUCCESS) PRINT_AND_RETURN("remove label failed in make_subtree\n", GENERAL_ERROR);

    tmp_options.input_name = label;
    tmp_options.output_name = outname;
    tmp_options.tree_names = &tree_name;
    if(nw_utils_job(&tmp_options) != SUCCESS) PRINT_AND_RETURN("nw_utils failed in make_subtree\n", GENERAL_ERROR);

    return 0;
}

int make_fasttree_constraint(char * msa_name, char * out_name, ml_options * master_ml_options){
    option_t tmp_options; 

    // Call fasttree
    tmp_options.input_name = msa_name;
    tmp_options.tree_names = &out_name;

    if(fasttree_job(&tmp_options, master_ml_options) != SUCCESS) PRINT_AND_RETURN("fasttree failed in main", GENERAL_ERROR);

    tmp_options.input_name = out_name;
    tmp_options.output_name = out_name;
    if(rm_label_job(&tmp_options) != SUCCESS) PRINT_AND_RETURN("remove label failed in main", GENERAL_ERROR);

    return 0;
}

int make_raxml_constraint(char * msa_name, char * out_name, ml_options * master_ml_options){
    option_t tmp_options; 

    char        raxml_name[10000];
    char        raxml_out_name[10000];
    char        raxml_dir_name[10000];

    // Find the out_name and dir_name for RAxML (this is required only for RAxML)
    find_prefix_and_dir_from_path(out_name, raxml_out_name, raxml_dir_name);

#if 1
    tmp_options.input_name = msa_name;
    tmp_options.output_name = raxml_dir_name;
    tmp_options.tree_names = malloc(sizeof(char *));
    tmp_options.tree_names[0] = raxml_out_name;
    if(raxml_job(&tmp_options, master_ml_options) != SUCCESS) PRINT_AND_RETURN("raxml_job failed in main", GENERAL_ERROR);
#else
    tmp_options.input_name = msa_name;
    tmp_options.output_name = raxml_dir_name;
    tmp_options.tree_names = malloc(sizeof(char *));
    tmp_options.tree_names[0] = raxml_out_name;
    if(raxml_with_initial_tree_job(&tmp_options, master_ml_options) != SUCCESS) PRINT_AND_RETURN("raxml_job failed in main", GENERAL_ERROR);
#endif

    sprintf(raxml_name, "%sRAxML_bestTree.%s", raxml_dir_name, raxml_out_name);
    tmp_options.input_name = raxml_name;
    tmp_options.output_name = out_name;
    if(rm_label_job(&tmp_options) != SUCCESS) PRINT_AND_RETURN("remove label failed in main", GENERAL_ERROR);
    system("rm RAxML*");

    return 0;
}



int constraint_inc(int argc, ml_options * master_ml_options){
    int i, j, read_mode, len;
    char command[GENERAL_BUFFER_SIZE];
    char name[GENERAL_BUFFER_SIZE];

    int arg_count = 6 + argc + master_ml_options->use_four_point_method_with_tree; 
    char ** argv = malloc(arg_count * sizeof(char*));

    for(i = 0; i < arg_count; i++){
        argv[i] = malloc(10000 * sizeof(char));
        argv[i][0] = 0;
    }
    
    if(master_ml_options->use_four_point_method_with_tree){
        sprintf(command, "constraint_inc -i %s -o %s -t %s", master_ml_options->init_d_name, master_ml_options->output_prefix, master_ml_options->guide_tree_name);
    } else {
        sprintf(command, "constraint_inc -i %s -o %s -t", master_ml_options->init_d_name, master_ml_options->output_prefix);
    }

    for(i = 0; i < argc; i++){
        sprintf(name, "%s_ctree%d.tree", master_ml_options->output_prefix, i);
        add_command(command, name);
    }

    printf("constraint_inc called with\n%s\n", command);

    j = 0;
    len = 0;
    read_mode = 1;
    for(i = 0; i < strlen(command); i++){
        if(!read_mode){
            if(command[i] == ' ') continue; // still blank
            else{
                argv[j][len] = 0;
                len = 0;
                j++;
                argv[j][len++] = command[i];
                read_mode = 1;
            }
        } else { // still reading
            if(command[i] == ' ')
                read_mode = 0;
            else 
                argv[j][len++] = command[i]; 
        }   
        
    }
    argv[j][len] = 0;
// printf("argc is %d\n", arg_count);
//     for(int i  = 0; i < argc; i++){
//         printf("argv is %s\n", argv[i]);
//     }
    constraint_inc_main(arg_count, argv, master_ml_options);

    return 0;
}

int unweighted_job(option_t * options){
    char command[GENERAL_BUFFER_SIZE];
    sprintf(command, "tree_to_dist.py -t %s > %s", options->input_name, options->output_name);
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling unweighted job\n", GENERAL_ERROR);
    return 0;
}

int fasttree_job(option_t * options,  ml_options * master_ml_options){
    char command[GENERAL_BUFFER_SIZE];
    char * gtrgamma_str = "-gtr -gamma";
    char * jc_str = " ";
    sprintf(command, "%s -nt %s -quiet < %s > %s", FastTree_bin, strcmp(master_ml_options->distance_model, "JC") ? gtrgamma_str : jc_str, options->input_name, options->tree_names[0]);
    printf("fasttree was called as %s\n", command);
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling fasttree_job\n", GENERAL_ERROR);
    return 0;
}

int fasttree_initial_tree_job(option_t * options,  ml_options * master_ml_options){
    FILE * f, * p;
    // int i;
    char command[100 * GENERAL_BUFFER_SIZE];
    char * gtrgamma_str = "-gtr -gamma";
    char * jc_str = " ";
    sprintf(command, "%s -nt %s -quiet -noml -nni 0 -log tmp.log < %s", FastTree_bin, strcmp(master_ml_options->distance_model, "JC") ? gtrgamma_str : jc_str, options->input_name);
    printf("fasttree was called as %s\n", command);
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling fasttree_job\n", GENERAL_ERROR);

    f = fopen("tmp.log", "r");
    p = fopen(options->tree_names[0], "w");
    while(fscanf(f, "%s", command) >=0){
        if(strcmp("NJ", command) == 0){
            fscanf(f, "%s", command); 
            // printf("%s %d %d\n", command, command[0], command[1]);
            // for(i = 0; i < )
            fprintf(p, "%s\n", command);
            break;
        }
    }
    fclose(f);
    fclose(p);
    // while(1);
    system("rm tmp.log");
    return 0;
}

int upp_job(option_t * options,  ml_options * master_ml_options){
    char command[GENERAL_BUFFER_SIZE];
    sprintf("run_upp.py -a %s -t %sfirst_tree.tree -s %s -A %d -S centroid", options->input_name, options->output_name, options->input_name, master_ml_options->ss_threshold);
    
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling upp_job\n", GENERAL_ERROR);
    return 0;
}

int subset_job(option_t * options,  ml_options * master_ml_options){ 
    char command[GENERAL_BUFFER_SIZE];
    sprintf(command, "build_subsets_from_tree.py -t %s -o %s -n %d", options->input_name, options->output_name, master_ml_options->ss_threshold);
    printf("command is %s\n", command);
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in subset_jon \n", GENERAL_ERROR);
    return 0;
}

int nw_utils_job(option_t * options){
    char command[GENERAL_BUFFER_SIZE];
    sprintf(command, "nw_prune -v %s $(cat %s) > %s", options->tree_names[0], options->input_name, options->output_name);
    printf("%s\n", command);
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling nw_utils\n", GENERAL_ERROR);
    return 0;
}

int rm_label_job(option_t * options){
    char command[GENERAL_BUFFER_SIZE];
    sprintf(command, "remove_internal_labels.py -t %s -o %s", options->input_name, options->output_name);
    // printf("%s\n", command);
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling rm_label_job\n", GENERAL_ERROR);
    return 0;
}

int raxml_job(option_t * options,  ml_options * master_ml_options){
    char command[GENERAL_BUFFER_SIZE];
    char * gtrgamma_str = "GTRGAMMA";
    char * jc_str = "GTRCAT -V --JC69";

    if(options->output_name){
        sprintf(command, "%s -T 12 --silent -m %s -j -n %s -s %s -w %s -p 1 > rubbish",     RAxML_bin, strcmp(master_ml_options->distance_model, "JC") ? gtrgamma_str : jc_str, options->tree_names[0], options->input_name, options->output_name);
    }
    else{
        sprintf(command, "%s -T 12 --silent -m %s -j -n %s -s %s -p 1 > rubbish",    RAxML_bin, strcmp(master_ml_options->distance_model, "JC") ? gtrgamma_str : jc_str, options->tree_names[0], options->input_name);
    }
    
    //printf("raxml was called as %s\n", command);
    // while(1);
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling raxml job\n", GENERAL_ERROR);
    return 0; 
}

// This also deletes the initial tree 
int raxml_with_initial_tree_job(option_t * options,  ml_options * master_ml_options){
    char command[GENERAL_BUFFER_SIZE];
    char * gtrgamma_str = "GTRGAMMA";
    char * jc_str = "GTRCAT -V --JC69";

    option_t tmp_options;
    tmp_options.input_name = options->input_name;
    tmp_options.tree_names = malloc(sizeof(char*));
    tmp_options.tree_names[0] = malloc(1000);
    tmp_options.tree_names[0] = "tmp.fast";

    if(fasttree_initial_tree_job(&tmp_options, master_ml_options) != SUCCESS) PRINT_AND_RETURN("fast fasttre job gfailed in raxml job with initial tree\n", GENERAL_ERROR);

    if(options->output_name){
        sprintf(command, "%s -T 12 --silent -m %s -j -n %s -s %s -w %s -t tmp.fast -p 1 > rubbish", RAxML_bin, strcmp(master_ml_options->distance_model, "JC") ? gtrgamma_str : jc_str, options->tree_names[0], options->input_name, options->output_name);
    }
    else{
        sprintf(command, "%s -T 12 --silent -m %s -j -n %s -s %s -t tmp.fast -p 1 > rubbish",       RAxML_bin, strcmp(master_ml_options->distance_model, "JC") ? gtrgamma_str : jc_str, options->tree_names[0], options->input_name);
    }
    
    //printf("raxml was called as %s\n", command);
    //while(1);
    system(command);// != SUCCESS)          PRINT_AND_RETURN("error in calling raxml job\n", GENERAL_ERROR);
    system("rm tmp.fast");
    return 0; 
}

int get_ll_from_raxml(option_t * options,  ml_options * master_ml_options, char * constraint_quartet, double* lp, double* l1, double* l2){
    FILE * f; 
    char tmp[GENERAL_BUFFER_SIZE];
    double tmp_d;
    int tmp_i;
    int counter = 0;
    sprintf(tmp, "RAxML_info.%s", options->tree_names[0]);

    raxml_with_quartet_tree_job(options, master_ml_options, constraint_quartet);

    // while(1);

    f = fopen(tmp, "r");
    while(fgets(tmp, sizeof(tmp), f)){
        // printf("%s\n", tmp);
        if(sscanf(tmp, "%d %lf", &tmp_i, &tmp_d) == 2){
            if(tmp_i == 0) *lp = tmp_d;
            if(tmp_i == 1) *l1 = tmp_d;
            if(tmp_i == 2) *l2 = tmp_d;
            counter ++;
        }
    }

    // f = fopen(tmp ,"r");
    // while(fscanf(f, "%lf", &tmp_d) >=0){
    //     counter++;
    //     if(counter % 3 == 2) *ll = tmp_d;
    // }
    fclose(f);


    system("rm RAxML_*");
    return 0;
}

int raxml_with_quartet_tree_job(option_t * options,  ml_options * master_ml_options, char * constraint_quartet){
    FILE * f;
    char command[GENERAL_BUFFER_SIZE];
    char * gtrgamma_str = "GTRGAMMA";
    char * jc_str = "GTRCAT -V --JC69";

// printf("%s\n", constraint_quartet);
    f = fopen("tmp.fast", "w"); // this will overwrite tmp.fast
    fprintf(f, "%s\n", constraint_quartet); // this may not be secure
    fclose(f);

    if(options->output_name){
        sprintf(command, "%s -T 12 --silent -m %s -j -n %s -s %s -w %s -f N -z tmp.fast -e 0.0001 -p 1 > rubbish", RAxML_bin, strcmp(master_ml_options->distance_model, "JC") ? gtrgamma_str : jc_str, options->tree_names[0], options->input_name, options->output_name);
    }
    else{
        sprintf(command, "%s -T 12 --silent -m %s -j -n %s -s %s -f N -z tmp.fast -e 0.0001 -p 1 > rubbish",       RAxML_bin, strcmp(master_ml_options->distance_model, "JC") ? gtrgamma_str : jc_str, options->tree_names[0], options->input_name);
    }
    // printf("here\n");
    // printf("raxml was called as %s\n", command);
    // while(1);
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling raxml in raxml_with quartet tree\n", GENERAL_ERROR);
    // while(1);
    system("rm tmp.fast");
    return 0; 
}


int distance_matrix_job(option_t * options, ml_options * master_ml_options){
    char command[GENERAL_BUFFER_SIZE];
    char cwd[PATH_MAX];

    if (getcwd(cwd, sizeof(cwd)) != NULL) {
       printf(stdout, "Current working dir: %s\n", cwd);
    }

    // This is too long
    sprintf(command, "echo \"ToNEXUS format=FASTA fromFile=%s toFile=%s/nexus; Quit;\" | %s -n;", options->input_name, PAUP_bin, cwd);
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling distance matrix job\n", GENERAL_ERROR);

    sprintf(command, "echo \"exe nexus; DSet distance=%s; SaveDist format=RelPHYLIP file=%s triangle=both diagonal=yes; Quit;\" | %s -n", master_ml_options->distance_model, master_ml_options->init_d_name, PAUP_bin);
    if(system(command) != SUCCESS)          PRINT_AND_RETURN("error in calling distance matrix job\n", GENERAL_ERROR);
    return 0;
}

// Internal functions

int find_prefix_and_dir_from_path(char * path, char * prefix, char * dir){
    int i, j;
    // Malloc sequence, assuming path is created
    if(!path)                               PRINT_AND_RETURN("path is NULL in find_prefix_and_dir_from_path\n", GENERAL_ERROR);

    // Init
    prefix[0] = 0;
    dir[0] = 0;

    // Find the last backslash, this separates the dir from the prefix
    for(i = strlen(path) - 1; i >= 0; i--)
        if(path[i] == '/') break;
    
    // Do a string copy
    for(j = 0; j < strlen(path); j++){
        if(j <= i) dir[j] = path[j];
        else prefix[j - i - 1] = path[j];
    }
    dir[i + 1] = 0;
    prefix[strlen(path) - i - 1] = 0;
    return 0; 
}

// int extract_starting_tree_from_fasttree(char * filename){
//     // Call FastTree2
//     option_t options;   
//     FILE * f, * p;
//     char name[100 * GENERAL_BUFFER_SIZE];

//     options->input_name = NULL;     // task holder here
//     options->output_name = NULL;    // task holder here
//     options->tree_names = malloc(sizeof(char*));
//     options->tree_names[0] = malloc(100);

//     fasttree_job(&options);

//     f = fopen(options->output_name, "r");
//     p = fopen(filename, "w");
//     if(!f) PRINT_AND_RETURN("file open failure for extract_starting_tree_from_fasttree", GENERAL_ERROR);

//     while(fscanf(f, "%s", name) >= 0){
//         if(strcmp(name, "NJ")){
//             fscanf(f, "%s", name);
//             fprintf(p, "%s\n", name);
//             fclose(f); 
//             fclose(p);
//         }
//         break;
//     }
//     return 0;
// }

// int extract_raxml_


/* Helper function to add a command into a string of commands
 * Input:   current string of command and the new command to be added 
 * Output:  none
 * Effect   none
 */
void add_command(char * current_command, char * new_command){
    str_add_space(current_command);
    strcat(current_command, new_command);
}

