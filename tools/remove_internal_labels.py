#!/Users/lethien_96/anaconda/bin/python

import argparse
import dendropy
import sys

def get_newick_string(tree):
    """
    Creates a newick string with  branch lengths in decimal notation.
    This function was taken from Phybase.py, see
    github.com/ngcrawford/CloudForest/blob/master/cloudforest/phybase.py

    Parameters
    ----------
    tree : dendropy tree

    Returns
    -------
    new_tree : string
               newick representation of tree
    """
    from decimal import Decimal

    str_newick_tree = tree.as_string(schema='newick')

    colon_s = 0
    comma_back_paren_s = 0
#    num = ''
    new_tree = ''
    for count, char in enumerate(str_newick_tree):
        if char == ':':
            colon_s = count
            continue
        if char in (')', ','):
            comma_back_paren_s = 1
#            num = '%1.12f' % Decimal(num)
#            new_tree += ":" + num
            colon_s = 0
#            num = ''
#        if colon_s != 0:
#            num = num + char
        if colon_s == 0:
            new_tree += char
    new_tree = new_tree.strip('\'').strip('\"').strip('\'') + '\n'
    return new_tree

def main(args):
    taxa = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get(path=args.input,
                         schema='newick',
                         taxon_namespace=taxa)
    tree.resolve_polytomies(limit=2)

    # Remove branch lengths
    for e in tree.preorder_edge_iter():
        e.length = None
    for inode in tree.internal_nodes():
        inode.label = None



    if args.output is not None:
        with open(args.output, 'w') as f:
            f.write(get_newick_string(tree))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-t", "--input", type=str,
                        required=True,
                        help="Input tree file")

    parser.add_argument("-o", "--output", type=str,
                        required=True,
                        help="Output prefix")

    main(parser.parse_args())