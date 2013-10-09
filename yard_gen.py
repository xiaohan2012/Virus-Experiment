from ve.clustering.lmatrix import LMatrix
from ve.clustering.load_data import load_mat

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

    
    
def print_yard_file(gmat, mat, method):
    print("output\t%s" %method)

    for pdb1 in gmat.rlabels:
        for pdb2 in gmat.clabels:
            print("%f\t%f" %(gmat[pdb1, pdb2] if pdb1 != pdb2 else 1.0 , mat[pdb1, pdb2]))

def print_overall_yards(gmat, mats, methods):
    print("output\t%s" %"\t".join(methods))
        
    for pdb1 in gmat.rlabels:
        for pdb2 in gmat.clabels:
            print("%f\t%s" %(gmat[pdb1, pdb2] if pdb1 != pdb2 else 1.0,
                                 "\t".join(map(lambda mat: "%f" %mat[pdb1,pdb2], mats))) )
    
def single_yard():
    import sys

    mat_path, method = sys.argv[1:]

    gmat = load_groups("data/manual_classification_result/166.tab")
    mat = load_mat(mat_path)

    if "RMSD" in mat_path:
        mat = 1 - mat / mat.max(1)
    
    print_yard_file(gmat, mat, method)

def main():
    gmat = load_groups("data/manual_classification_result/166.tab")

    mat_ids = ["CEpiMatch", "MATT", "Multiprot", "SPa", "SPb", "SPe", "TMa", "TMb", "TMc", "RMSD"]
    mats = map(lambda mat_id: load_mat("clustering/data/%s.csv" %mat_id), mat_ids)

    mats[-1] = 1 - mats[-1] / mats[-1].max(1)
    
    print_overall_yards(gmat, mats, mat_ids)
    
if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    #main()
    single_yard()