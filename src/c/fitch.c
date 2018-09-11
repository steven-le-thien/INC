// // Code written by Thien Le in September 2018 implementing Fitch algorithm. The main algorithm is implemented in the first part,
// // followed by a unit test suite. For more information on how to run the unit test, refer to that section.

// // Some limits: 
// //      The maximum number of nodes (not just leaves) are hardcoded to 1M. This is to facillitate building arrays on the stack 
// // at compile time. If more than 1M nodes is required, may need to switch to heap implementation (which is slower), or have the compiler
// // extend the stack space.
// //      The maximum number of letters in a leaf name is hardcoded to 10. Again, this is to facillitate arrays on stack. 
// //      Labels must be integer no larger than 31. This is to facillitate bitwise operations in the implementation. Extension is trivial 
// // but still require some work.

// // Compilation: you may need to use -lm flag for some version of GCC since some of the math.h functions would not be imported otherwise. 

// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <math.h>

// #define UNIT_TEST   0

// // Global variable to quickly have a tree on the stack 
// // If the tree is huge (> 10^6) then may want to switch to heap implementation; but will take more time since malloc is slow
// const int N_LIM = (int) 10e5;       // maximum number of nodes supported is 100K, which corresponds to 25K-taxon tree
// const int NAME_LIM = 10;            // maximum name length supported is 10
// const int LABEL_LIM = (int) 10e3;   // maximum number of distinct labels supported is 31K, represented by 1000 ints

// int root;

// int left    [N_LIM];    // left child node
// int right   [N_LIM];    // right child node
// int parent  [N_LIM];    // parent node
// int label   [N_LIM];    // digitized taxon id
// int set     [N_LIM][LABEL_LIM];    // bitset parsimony where bit i is set if label i is in the parsimony set   

// char name_map   [N_LIM][NAME_LIM];   
// double weight   [N_LIM][2]; // adjacency list style weight, the left edge is the first entry, the right edge is the second (currently not used for anything)

// // Set a bit in the concatenated array at the position specified in `index'
// int set_bit(int * array, int index){
//     if(index / 32 >= LABEL_LIM) // index exceeds array capacity
//         return -1;
//     else 
//         array[index / 32] |= (1 << index);
//     return 0;
// }

// // Fitch algorithm
// int post_order_dfs(int c){
//     int i;
//     int cap[LABEL_LIM], cup[LABEL_LIM];

//     if(left[c] == -1){ // at leaves
//         set_bit(set[c], label[c]);
//         return 0;
//     }

//     post_order_dfs(left[c]);
//     post_order_dfs(right[c]);
    
//     for(i = 0; i < LABEL_LIM; i++){ // this is easily parallelizable
//         cap[i] = set[left[c]][i] & set[right[c]][i];
//         cup[i] = set[left[c]][i] | set[right[c]][i];
//         set[c][i] = cap[i] ? cap[i] : cup[i];
//     }

//     return 0;
// }

// int pre_order_dfs(int c){
//     int cap = 0;

//     if(left[c] == -1) // at leaves
//         return 0;

//     if(parent[c] == -1 ||                           // at root and short circuit 
//             !(cap = set[c] & set[parent[c]]))       // or empty intersection
//         set[c] &= -set[c];                              // clear all bits but the LSB
//     else
//         set[c] = cap;                               // else, set it as the intersection

//     label[c] = (int) (log(set[c]) / log(2.0));

//     pre_order_dfs(left[c]);
//     pre_order_dfs(right[c]);

//     return 0;
// }

// int fitch_algorithm(){
//     post_order_dfs(root);
//     pre_order_dfs(root);
//     return 0;
// }

// // Unit test suite.
// // Input:   pipe into stdin a Newick string with leaf node of the form <name>?<label>:<edge_weight> 
// //              and internal node of the form <subtree>:<edge_weight>
// // Output:  the same tree with internal node labelled by the Fitch algorithm is printed to stdout.
// // Example input:   ((A?10:0.1,B?5:1.2):10,(C?30:0.9,(D?2:0.0,ADF?30:1.2):0.6):0.6);
// // Example output:  ((A?10:0.100000,B?5:1.200000)?5:10.000000,(C?30:0.900000,(D?2:0.000000,ADF?30:1.200000)?30:0.600000)?30:0.600000)?5;
// //                         notice that internal nodes, as well as the root node now have labels. 
// // 

// //--------------------- UNIT TEST CODE
// #if UNIT_TEST

// #define DEBUG       0

// // Unit test initialization
// const int READ_WEIGHT_STATE = 2;
// const int READ_LABEL_STATE  = 1;
// const int READ_NAME_STATE   = 0;
// const int OTHER_STATE       = 3;

// void write_to_mem(int parent, char * buf, int state){
//     int node;

//     if(cur_state != OTHER_STATE) return;

//     node = right[parent] == -1 ? left[parent] : right[parent];

//     if(buf[0] == 0) return;

//     if(state == READ_NAME_STATE)
//         strcpy(name_map[node], buf);
//     else if (state == READ_LABEL_STATE)
//         sscanf(buf, "%d", &label[node]);
//     else if (state == READ_WEIGHT_STATE)
//         sscanf(buf, "%lf", &weight[parent][right[parent] != -1]); // if the right parent exists then write to the right
    
//     buf[0] = 0;
// }

// void inc_node(int * node_counter, int * cur_state){
//     (*node_counter)++;
//     *cur_state = READ_NAME_STATE;
// }

// void set_parent(int par, int chi){
//     parent[chi] = par;
//     if(left[par] == -1)
//         left[par] = chi;
//     else 
//         right[par] = chi;
// }

// void init(){
//     char cur_char;
//     char cur_mem[3][10000]; // 1. name; 2. label; 3. weight
//     int node_counter    = 0;
//     int cur_parent      = -1;
//     int cur_state       = OTHER_STATE;

//     // Clear all variables
//     root = 0;
//     cur_mem[0][0] = 0;
//     cur_mem[1][0] = 0;
//     cur_mem[2][0] = 0;
//     memset(left,    -1, N_LIM * sizeof(int));
//     memset(right,   -1, N_LIM * sizeof(int));
//     memset(parent,  -1, N_LIM * sizeof(int));
//     memset(label,   -1, N_LIM * sizeof(int));
//     memset(set,     0,  N_LIM * sizeof(int));


//     // Read from stdin a Newick string a Newick string
//     while(scanf("%c", &cur_char) == 1){
//                                                                                             #if DEBUG
//                                                                                                  printf("debug: in loop, cur_label is %s, char is %c\n", cur_mem[2], cur_char); 
//                                                                                             #endif
//         switch(cur_char){
//             case '(': 
//                 cur_parent = node_counter;              // increment level 
//                 inc_node(&node_counter, &cur_state);
//                 set_parent(cur_parent, node_counter);
//                 break;
//             case ',': 
//                 write_to_mem(cur_parent, cur_mem[cur_state], cur_state);
//                 inc_node(&node_counter, &cur_state);
//                 set_parent(cur_parent, node_counter);
//                 break;
//             case ')':
//                 write_to_mem(cur_parent, cur_mem[cur_state], cur_state);
//                 cur_parent = parent[cur_parent];        // decrement level
//                 cur_state = OTHER_STATE;
//                 break;
//             case ':': 
//                 write_to_mem(cur_parent, cur_mem[cur_state], cur_state);
//                 cur_state = READ_WEIGHT_STATE;
//                 break;
//             case '?':
//                 write_to_mem(cur_parent, cur_mem[cur_state], cur_state);
//                 cur_state = READ_LABEL_STATE;
//                 break;
//             case ';':
//                 break;
//             default:
//                 if(cur_state != OTHER_STATE){
//                     cur_mem[cur_state][strlen(cur_mem[cur_state]) + 1] = 0;
//                     cur_mem[cur_state][strlen(cur_mem[cur_state])] = cur_char;
//                 }
//         }
//     }
// }

// void write_dfs(int c, char * builder){
//     char buffer[10000];

//     if(left[c] == -1){
//         strcat(builder, name_map[c]);
//         return;
//     }

//     // Left child
//     strcat(builder, "(");

//     write_dfs(left[c], builder);
//     sprintf(buffer, "?%d:%lf", label[left[c]], weight[c][0]);
//     strcat(builder, buffer);

//     // Right child
//     strcat(builder, ",");
//     write_dfs(right[c], builder);
//     sprintf(buffer, "?%d:%lf", label[right[c]], weight[c][1]);
//     strcat(builder, buffer);

//     strcat(builder, ")");
// }

// int write_newick(){
//     char buffer[100000];
//     buffer[0] = 0;

//     write_dfs(0, buffer);
//     printf("%s?%d;\n", buffer, label[0]);
//     return 0;
// }


// int main(){
//     init();
//                                                                                             #if DEBUG
//                                                                                                 int i;
//                                                                                                 int n = 50;

//                                                                                                 printf("debug: this print out the parent, left, right, label, name_map and set\n");
//                                                                                                 for(i = 0; i < n; i++){
//                                                                                                     printf("%d ", parent[i]);
//                                                                                                 }
//                                                                                                 printf("\n"); 
//                                                                                                 for(i = 0; i < n; i++){
//                                                                                                     printf("%d ", left[i]);
//                                                                                                 }
//                                                                                                 printf("\n"); 
//                                                                                                 for(i = 0; i < n; i++){
//                                                                                                     printf("%d ", right[i]);
//                                                                                                 }
//                                                                                                 printf("\n"); 
//                                                                                                 for(i = 0; i < n; i++){
//                                                                                                     printf("%d ", label[i]);
//                                                                                                 }
//                                                                                                 printf("\n"); 
//                                                                                                 for(i = 0; i < n; i++){
//                                                                                                     printf("%s ", name_map[i]);
//                                                                                                 }
//                                                                                                 printf("\n"); 
//                                                                                                 for(i = 0; i < n; i++){
//                                                                                                     printf("%d ", set[i]);
//                                                                                                 }
//                                                                                                 printf("\n"); 
//                                                                                             #endif
//     fitch_algorithm();
//                                                                                             #if DEBUG
//                                                                                                 // int i;
//                                                                                                 // int n = 50;

//                                                                                                 printf("debug: this print out the parent, left, right, label, name_map and set\n");
//                                                                                                 for(i = 0; i < n; i++){
//                                                                                                     printf("%d ", parent[i]);
//                                                                                                 }
//                                                                                                 printf("\n"); 
//                                                                                                 for(i = 0; i < n; i++){
//                                                                                                     printf("%d ", left[i]);
//                                                                                                 }
//                                                                                                 printf("\n"); 
//                                                                                                 for(i = 0; i < n; i++){
//                                                                                                     printf("%d ", right[i]);
//                                                                                                 }
//                                                                                                 printf("\n"); 
//                                                                                                 for(i = 0; i < n; i++){
//                                                                                                     printf("%d ", label[i]);
//                                                                                                 }
//                                                                                                 printf("\n"); 
//                                                                                                 for(i = 0; i < n; i++){
//                                                                                                     printf("%s ", name_map[i]);
//                                                                                                 }
//                                                                                                 printf("\n"); 
//                                                                                                 for(i = 0; i < n; i++){
//                                                                                                     printf("%d ", set[i]);
//                                                                                                 }
//                                                                                                 printf("\n"); 
//                                                                                             #endif

//     write_newick();
// }
// #endif
