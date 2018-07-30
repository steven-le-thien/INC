### Constraint INC

The INC tree building technique is described in `R. Zhang, S. Rao, and T. Warnow (2018). New absolute fast converging phylogeny estimation methods with improved scalability and accuracy. To appear, Proceedings of WABI (Workshop on Algorithms for Bioinformatics) 2018'. It is an absolute fast converging (AFC) method for the problem of estimating the model tree given a distance matrix between the sequences. 

Constraint INC is described in the same paper. In addition to a distance matrix input, it takes in a set of leaf-disjoint constraint trees and output a tree over the union of all leafsets that agrees with the constraint trees. When the constraint trees are derived under AFC methods (such as maximum likelihood), constraint INC is also AFC. 

The implementation design can be summarized as. 1/ Build an MST from the distance matrix and obtain the Prim ordering and the maximum edge weight, 2/ Build the tree incrementally following the Prim ordering: 2.1/ Start with a 3-leaf tree that is the first 3 taxa in the Prim ordering. 2.2/ For each new leaf added, identify the constraint tree it comes from, then use information in the constraint tree to limit the subtree on the growing tree that it can be added to, in order to agree with the constraint tree. 2.3/ Use the maximum edge weight of the MST to perform quartet `voting' and select the edge with the highest vote to add the new taxon into.  

## Requirement
Constraint INC was built from standard C packages 
1. Make
2. GCC or other C compile (code was developed and tested with GCC) 

If you want to use the maximum likelihood extension of constraint INC then the following is also currently required:
1. PASTA (working to remove this dependency) (https://github.com/smirarab/pasta) 
2. FastTree2 (http://www.microbesonline.org/fasttree/FastTree.c)
3. Newick Utils (https://github.com/tjunier/newick_utils)
4. Extra scripts in the tools folder (currently not published for liscensing reasons)

## Installation
Run `make' to generate the binaries `constraint_inc' and `ml', used for the constraint INC and its ML extension respectively. The command for constraint INC is 
```
constraint_inc -i <input_distance_matrix> -o <output_prefix> -t <constraint_tree1> <constraint_tree2> ... 
'''
The command for INC-ML is 
```
ml -i <input_alignment> -o <output_prefix> -t <initial_tree>
'''
The `-t' flag to INC-ML is optional
