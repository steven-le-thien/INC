import numpy as np

def parse_distance_matrix(dist_file):
    """
        Read the PHYLIP distance matrix DIST_FILE
        return:
            names   - an ordered list of the names in the file
            matrix  - a numpy distance matrix containing the distances
    """
    with open(dist_file) as f:
        n = int(f.readline().strip())
        names = []
        matrix = np.zeros((n,n))
        for i in range(n):
            line = f.readline().split()
            names.append(line[0])
            matrix[i,:] = line[1:]
    return names, matrix

def mat_to_func(names, matrix):
    """
        Given 0-indexed list of NAMES and corresponding 0-indexed numpy matrix D
        return a function mapping pairs of names to their corresponding distance
    """
    name2Ind = dict()
    for i in range(len(names)):
        name2Ind[names[i]] = i
    def D(n1=None,n2=None, get_mat=False):
        if not (n1 or n2) and get_mat:
            return matrix
        i1,i2 = name2Ind[n1], name2Ind[n2]
        return matrix[i1,i2]
    return D
