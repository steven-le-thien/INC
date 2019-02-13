#!/usr/bin/python
import argparse
import dendropy
import sys

def main(args):
	tree = dendropy.Tree.get_from_path(src=args.input_tree_file, schema="newick")
	for e in tree.postorder_edge_iter():
        	if e.length is not None:
                	e.length = 1

	for root in tree.preorder_node_iter():
		if len(root.child_edges()) == 2:
			for e in root.child_edge_iter():
				if e.length is not None:
					e.length = 0.5
		break

	pdc = tree.phylogenetic_distance_matrix()

#	print(pdm.as_string("phylip"))
	sys.stdout.write("%d\n" % len(tree.leaf_nodes())) 
	for i, t1 in enumerate(tree.taxon_namespace):
		sys.stdout.write("%s " % t1.label)
		for t2 in tree.taxon_namespace:
        		sys.stdout.write("%d " % pdc(t1, t2))# print("Distance between '%s' and '%s': %s" % (t1.label, t2.label, pdc(t1, t2)))
		sys.stdout.write("\n")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-t", "--input_tree_file", type=str,
                        required=True,
                        help="Input tree file")

    main(parser.parse_args())
