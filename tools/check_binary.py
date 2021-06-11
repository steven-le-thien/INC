#!/usr/bin/python

import argparse
import dendropy
import sys

def main(args):
    taxa = dendropy.TaxonNamespace()
    tree = None
    proceed = True
    for i in range(len(args.input)):
        if proceed == False:
            break 
        tree = dendropy.Tree.get(path=args.input[i],
                                schema='newick',
                                taxon_namespace=taxa)
        tree.collapse_basal_bifurcation()
        proceed = True
        for n in tree.preorder_node_iter():
            degree = len(n.adjacent_nodes()) 
            if degree == 2:
                print("Error: There is a node with degree 2 that is not the root")
                proceed = False
                break

            elif degree > 3:
                print("Error: There is a node with degree more than 3")
                proceed = False
                break
                
    with open(args.output, 'w') as f:
        if(proceed):
            f.write("1")
        else:
            f.write("0")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-t", "--input", type=str,
                        nargs='+', required=True,
                        help="Input tree files")

    parser.add_argument("-o", "--output", type=str,
                        required=True,
                        help="Output prefix")

    main(parser.parse_args())
