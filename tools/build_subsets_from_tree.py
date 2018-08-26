#!/Users/lethien_96/anaconda/bin/python

"""
Build subsets using centroid decomposition from PASTA

Written by EKM (molloy.erin.k@gmail.com) in Spring 2018.
"""
import dendropy
import argparse
from pasta.pastaalignerjob import bisect_tree
from pasta.tree import PhylogeneticTree
import sys


def decompose_trees(tree, max_subset_size):
    """
    """
    next_trees = [tree]
    done_trees = []
    while len(next_trees) > 0:
        trees = next_trees
        next_trees = []
        for tree in trees:
            t1, t2 = bisect_tree(tree,
                                 breaking_edge_style="centroid")
            n1 = t1.n_leaves
            n2 = t2.n_leaves
            if n1 > max_subset_size:
                next_trees.append(t1)
            else:
                if n1 >= 5:
                    done_trees.append(t1)
                else:
                    sys.exit("T1 has fewer than 5 leaves!")
            if n2 > max_subset_size:
                next_trees.append(t2)
            else:
                if n2 >= 5:
                    done_trees.append(t2)
                else:
                    sys.exit("T2 has fewer than 5 leaves!")
    return done_trees


def main(args):
    # Step 1: Decompose tree
    tree = dendropy.Tree.get(path=args.input_tree_file,
                             schema="newick")
    tree.resolve_polytomies(limit=2,
                            update_bipartitions=True)
    tree = PhylogeneticTree(tree)
    trees = decompose_trees(tree,
                            args.max_subset_size)
    # Step 2: Write out leaf subsets
    n = len(trees)
    pad = len(str(n))
    # i = 1
    i = 0
    for tree in trees:
        keep = tree.leaf_node_names()
        j = str(i).zfill(pad)
        # output_file = args.output + "-" + j + "-outof-" + str(n) + ".txt"
        with open(args.output + "_ctree" + str(i) + ".lab", 'w') as f:
            f.write("\n".join(keep))
        i = i + 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-n", "--max_subset_size", type=int,
                        required=True,
                        help="Maximum subset size")

    parser.add_argument("-t", "--input_tree_file", type=str,
                        required=True,
                        help="Input tree file")

    parser.add_argument("-o", "--output", type=str,
                        required=True,
                        help="Output prefix")

    main(parser.parse_args())
