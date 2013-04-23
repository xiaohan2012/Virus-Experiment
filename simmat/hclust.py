"""
Perform hierarchical clustering on similarity matrix
"""
from source import load_simmat

def sim_to_dist(mat):
    """
    (matrix) -> matrix
    invert the similarity matrix to distance matrix
    """
    return 1 - mat

from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib import pyplot as plt

def main():
    mat = load_simmat("fp_370_atg.txt")
    dist_mat = sim_to_dist(mat)
    link_mat = linkage(mat)
    dendrogram(link_mat)

if __name__  == "__main__":
    main()
