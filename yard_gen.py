from ve.clustering.lmatrix import LMatrix

def load_groups(path):
    """
    produce the relationship matrix
    
    >>> m = load_groups("data/manual_classification_result/166.tab")
    >>> m['3LHP_S', '3LH2_S']
    1.0
    >>> m['3LH2_S', '3LHP_S']
    1.0
    >>> m['3LH2_S', '2YBR_C']
    0.0
    """
    with open(path) as f:
        mapping = dict( map(tuple, map(lambda s: s.split(), f.readlines())) )
    
    pdbs = mapping.keys()

    mat = LMatrix(pdbs)
    
    from collections import defaultdict
    d = defaultdict(list)
    
    for pdb in pdbs:
        d[mapping[pdb]].append(pdb)

    from itertools import permutations
    
    for items in d.values():
        for pairs in permutations(items, 2):
            l,r = pairs
            mat[l,r]=1
            mat[r,l]=1
            
    return mat

    
    
def print_yard_file(gmat, mat):
    print("output\tmethod1")

    for pdb1 in gmat.rlabels:
        for pdb2 in gmat.clabels:
            print("%f\t%f" %(gmat[pdb1, pdb2] if pdb1 != pdb2 else 1.0 , mat[pdb1, pdb2]))

def main():
    import sys
    from ve.clustering.load_data import load_mat

    mat_path = sys.argv[1]

    gmat = load_groups("data/manual_classification_result/166.tab")
    mat = load_mat(mat_path)
    
    print_yard_file(gmat, mat)
    
if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    main()