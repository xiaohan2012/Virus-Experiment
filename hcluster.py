from scipy.cluster.hierarchy import linkage,fcluster,dendrogram,fclusterdata
from pickle import load
from numpy import array

sim_mat = array(load(open("dist_mat_402.dat")))
sim_mat = sim_mat / sim_mat.diagonal()
with open("normalized_402.dat",'w') as f:
    for row in sim_mat:
        f.write(','.join(map(str,row)) + "\n")
        

