"""
using hierarchical clustering to classify protein complexes into clusters
"""
from scipy.cluster.hierarchy import linkage,fcluster,dendrogram,fclusterdata
from pickle import load
from numpy import array,min,ceil
import matplotlib.pyplot as pp

sim_mat = array(load(open("dist_mat_epi_166.dat")))

sim_mat += ceil(min(sim_mat))#prevent negative

feature_vec= sim_mat / sim_mat.diagonal()

feature_vec = 1 / feature_vec#comment this out to use the feature vector directly

Z = linkage(feature_vec)

dendrogram(Z)

pp.show()

