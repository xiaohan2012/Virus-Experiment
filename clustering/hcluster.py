"""
using hierarchical clustering to classify protein complexes into clusters
"""
from scipy.cluster.hierarchy import linkage,fcluster,dendrogram,fclusterdata
import matplotlib.pyplot as plt

from load_data import load_mat
from categories import categorize, data, data1

import numpy as np

def compare():
    """
    plot a series of hcluster dendragram for various methods on various sets of virus
    """
    
    paths = ['CEpiMatch.csv', 'Multiprot.csv', 'SPa.csv', 'SPb.csv', 'SPe.csv', 'TMa.csv', 'TMb.csv', 'TMc.csv', 'MATT.csv', "RMSD.csv"]

    for path in paths:
        mats = categorize(load_mat("data/%s" %path), data)
        
        for name, mat in mats.items():
            if path not in ("RMSD.csv", ):
                mat = 1 / mat
            
            Z = linkage(mat)

            plt.figure()
            
            dendrogram(Z, labels=mat.rlabels, orientation="right")

            type_name = path.split(".")[0]

            plt.title("%s  %s" %(name, type_name))
            plt.ylabel("PDB ID")
            plt.savefig("img/%s/%s.png" %( name, type_name))


def main():
    mat = load_mat("data/CEpiMatch.csv")
    mat = 1 / mat
    
    Z = linkage(mat)
    
    #plt.figure()

    dendrogram(Z, labels=mat.rlabels, orientation="right", leaf_font_size=5)

    plt.title("Hierarchical clustering result")
    plt.ylabel("PDB ID")
    plt.plot()
    #plt.show()
    plt.savefig("img/166.png")

def for_lysozyme():
    """
    try various version of distance matrix for only lysozyme
    """
    
    import sys
    method = sys.argv[1]
    
    paths = ['CEpiMatch.csv', 'Multiprot.csv', 'SPa.csv', 'SPb.csv', 'SPe.csv', 'TMa.csv', 'TMb.csv', 'TMc.csv', 'MATT.csv', 'RMSD.csv']

    for path in paths:
        mats = categorize(load_mat("data/%s" %path), data)
        
        mat = mats["lysozyme"]

        if path not in ("RMSD.csv", ):
            mat = 1 / mat
        
        if method == "upper":
            mat = np.triu(mat)
        elif method == "lower":
            mat = np.transpose(np.tril(mat))
        elif  method == "average":
            mat = (np.tril(mat) + np.triu(mat)) / 2
            
                
        plt.figure()
        
        Z = linkage(mat)
        
        dendrogram(Z, labels=mat.rlabels, orientation="right")
        
        plt.title("Hierarchical clustering for lysozyme")
        plt.ylabel("PDB ID")
        plt.savefig("img/for_lysozyme_%s/%s.png" %(method, path.split(".")[0]))

if __name__ == '__main__':
    #main()
    compare()
    #for_lysozyme()