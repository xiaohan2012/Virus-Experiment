"""
Perform hierarchical clustering on similarity matrix
"""
from common import load_simmat

def sim_to_dist(mat):
    """
    (matrix) -> matrix
    invert the similarity matrix to distance matrix
    """
    from ve.simmat.lmatrix import lmatrix
    return lmatrix(labels=mat.labels, data= 1 - mat)

from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib import pyplot as plt

import pylab 

def main():
    fig = pylab.figure(figsize = (80,40))
    
    mat = load_simmat("fp_370_atg.txt")
    print mat.labels
    dist_mat = sim_to_dist(mat)
    link_mat = linkage(mat)

    llf = lambda i: dist_mat.get_label(i)

    dendrogram(link_mat, leaf_label_func=llf)
    
    plt.tick_params(axis='both', which='major', labelsize=18)

    fig.show()
    fig.savefig("atg.png")

    """
    #show and save figure
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages('atg.pdf')
    pp.savefig(fig)
    pp.close()
    """
if __name__  == "__main__":
    main()
