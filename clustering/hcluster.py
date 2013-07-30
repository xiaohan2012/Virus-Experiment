"""
using hierarchical clustering to classify protein complexes into clusters
"""
from scipy.cluster.hierarchy import linkage,fcluster,dendrogram,fclusterdata
import matplotlib.pyplot as plt

from load_data import load_mat
from categories import categorize, data

def main():

    paths = ['CEpiMatch.csv', 'Multiprot.csv', 'SPa.csv', 'SPb.csv', 'SPe.csv', 'TMa.csv', 'TMb.csv', 'TMc.csv', 'MATT.csv']

    for path in paths:
        mats = categorize(load_mat("data/%s" %path), data)
        
        for name, mat in mats.items():
            mat = 1 / mat
            
            Z = linkage(mat)

            plt.figure()
            
            dendrogram(Z, labels=mat.rlabels, orientation="right")

            type_name = path.split(".")[0]

            plt.title("%s  %s" %(name, type_name))
            plt.ylabel("PDB ID")
            plt.savefig("img/%s/%s.png" %( name, type_name))


if __name__ == '__main__':
    main()