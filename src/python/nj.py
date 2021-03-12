from skbio import DistanceMatrix
from skbio.tree import nj

def NJ(names, matrix):
    """ 
        input: a numpy matrix 
        return: newick string corresponding to neighbor joining
    """
    dm = DistanceMatrix(matrix, names)
    newick_str = nj(dm, result_constructor=str)
    return newick_str