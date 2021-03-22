import dendropy

def RF_error(pred_tree, true_tree):
    """
        Given PRED_TREE and TRUE_TREE as newick strings, calculate the RF error
        
        but this does include the trivial cuts, so we get like 50% for free

        let FP be the number of false positives
        let FN be the number of false negatives
        let N be the total number of taxa
        return (FP + FN)/(2N-6)
    """
    pred_tree = dendropy.Tree.get(data=pred_tree,schema='newick')
    true_tree = dendropy.Tree.get(data=true_tree, schema='newick')

    pred_N, true_N = len(pred_tree.leaf_nodes()), len(true_tree.leaf_nodes())
    assert pred_N == true_N
    
    # migrate taxon namespace so nodes in each tree are labeled in the same way
    # so that the biparitions later code the same thing
    pred_tree.migrate_taxon_namespace(true_tree.taxon_namespace)

    pred_cuts = set(pred_tree.encode_bipartitions())
    true_cuts = set(true_tree.encode_bipartitions())

    FP, FN = 0, len(true_tree.edges())
    for cut in pred_cuts:
        FP += (cut not in true_cuts)
        FN -= (cut in true_cuts)
    return (FP + FN) / (2 * true_N - 6)

def calc_tpr(pred_tree, true_tree):
    """
        Given PRED_TREE and TRUE_TREE as newick strings, calculate the true positive rate
        of the predicted nontrivial edges out of the total true edges
    """

    pred_tree = dendropy.Tree.get(data=pred_tree,schema='newick')
    true_tree = dendropy.Tree.get(data=true_tree, schema='newick')
    # cannot work with empty or trivial topology trees
    assert len(pred_tree.leaf_nodes()) > 3
    assert len(true_tree.leaf_nodes()) > 3

    # migrate taxon namespace so nodes in each tree are labeled in the same way
    # so that the biparitions later code the same thing
    pred_tree.migrate_taxon_namespace(true_tree.taxon_namespace)

    pred_cuts = set(pred_tree.encode_bipartitions())
    true_cuts = set(true_tree.encode_bipartitions())
    assert true_cuts

    num_correct = 0 # count the number of nontrivial edges correct
    total = 0 # count to total nontrivial edges
    for cut in true_cuts:
        # ignore trivial cuts, which are always shared
        if not cut.is_trivial():
            total += 1
            if cut in pred_cuts:
                num_correct += 1
    return num_correct/total

def calc_accuracy(names, noisy_matrix, method, true_tree_newick, return_est=False):
    """ Calculate the accuracy of METHOD on the noisy_matrices, given that
            the newick string true_tree_newick
            
        The metric is the RF error rate (see above)        
        METHOD must take a numpy matrix and return a newick string
        """
    est_string = method(names, noisy_matrix)
    if return_est:
    	return RF_error(est_string, true_tree_newick), est_string
    return RF_error(est_string, true_tree_newick)



def extract_tree(tree, taxa_names):
    """
        Given a newick string tree and list of taxa strings,
        return a tree with leaves only in the list, that agrees with the topology of the input tree
    """
    taxa_names = set(taxa_names)
    dendro_tree = dendropy.Tree.get(data=tree,schema='newick')
    for node in dendro_tree.leaf_nodes():
        if node.taxon.label in taxa_names:
            break
        node = None

    if not node: # no taxa in tree with desired taxa names
        return "();"
    dendro_tree.reroot_at_node(node.adjacent_nodes()[0],
                               update_bipartitions=True,
                              suppress_unifurcations=False)
    taxa_to_retain = set([taxon for taxon in dendro_tree.taxon_namespace
        if taxon.label in taxa_names])
    sub_tree = str(dendro_tree.extract_tree_with_taxa(taxa=taxa_to_retain,suppress_unifurcations=False)) + ";"
    return sub_tree

def satisfies_constraint(tree_str, constraint_str):
    """
        Return whether TREE satisfies CONSTRAINT, both input as newick strings
    """
    tree = dendropy.Tree.get(data=tree_str, schema='newick')
    tree_taxa = set(taxon.label for taxon in tree.taxon_namespace)
    constraint = dendropy.Tree.get(data=constraint_str, schema='newick')
    constraint_taxa = set(taxon.label for taxon in constraint.taxon_namespace)
    if len(tree_taxa.intersection(constraint_taxa)) < 4:
        # constraint satisfied
        return True
    sub_tree = extract_tree(tree_str, constraint_taxa)
    sub_const = extract_tree(constraint_str, tree_taxa)
    return calc_tpr(sub_tree, sub_const) == 1

# unused because satisfies_constraint needs to consider the size of the constraints
def matching_sub_trees(tree1_str,tree2_str):
    tree1 = dendropy.Tree.get(data=tree1_str, schema='newick')
    tree1_taxa = [leaf.taxon.label for leaf in tree1.leaf_node_iter()]
    tree2 = dendropy.Tree.get(data=tree2_str, schema='newick')
    tree2_taxa = [leaf.taxon.label for leaf in tree2.leaf_node_iter()]
    return extract_tree(tree1_str, tree2_taxa),extract_tree(tree2_str, tree1_taxa)