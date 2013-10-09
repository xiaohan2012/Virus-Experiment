"""
Perform hierarchical clustering on similarity matrix
"""
from source import load_simmat

def sim_to_dist(mat):
    """
    (matrix) -> matrix
    invert the similarity matrix to distance matrix
    """
    from ve.simmat.lmatrix import lmatrix
    return lmatrix(labels=mat.labels, data= 1 - mat)

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from matplotlib import pyplot as plt

import pylab 


def get_groups(lk, llf, cutoff=1.1):
    """(np.array) => list of list of str
    given the linkage matrix, return a labeled groups
    """
    assignments = fcluster(lk, cutoff)
    cluster_ids = set(assignments)
    
    groups = []
    for c_id  in cluster_ids:
        indices = np.where(assignments == c_id)[0]
        groups.append(map(llf, indices))
    return groups

def gen_dendrogram(path, out):
    """(str, str) => None
    
    from complex pair-wise distance text file to dendagram
    """
    fig = pylab.figure(figsize = (80,40))
    
    mat = load_simmat(path)
    print mat.labels
    dist_mat = sim_to_dist(mat)
    link_mat = linkage(mat)

    llf = lambda i: dist_mat.get_label(i)

    dendrogram(link_mat, leaf_label_func=llf)
    
    plt.tick_params(axis='both', which='major', labelsize=18)

    fig.savefig(out)

def a_bunch_of_dendrograms():
    gen_dendrogram("data/first_110.atg.txt", "plot/first_110.atg.png")
    gen_dendrogram("data/second_110.atg.txt", "plot/second_110.atg.png")
    gen_dendrogram("data/last_150.atg.txt", "plot/last_150.atg.png")
    gen_dendrogram("data/fp_370_atg.txt", "plot/fp_370_atg.png")

    gen_dendrogram("data/first_110.atb.txt", "plot/first_110.atb.png")
    gen_dendrogram("data/second_110.atb.txt", "plot/second_110.atb.png")
    gen_dendrogram("data/last_150.atb.txt", "plot/last_150.atb.png")
    gen_dendrogram("data/fp_370_atb.txt", "plot/fp_370_atb.png")


def disp_linakge_mat(lk, llf, c_count):
    print "Cluster ID, Subset1 ID, Subset2 ID, Distance, Cluster Size"
    for i, l in enumerate(lk):
        id1, id2, dist, count = l
        id1, id2 = int(id1), int(id2)
        id1 = str(llf(id1) if id1 < c_count-1 else id1)
        id2 = str(llf(id2) if id2 < c_count-1 else id2)
        
        print "%d,%s,%s,%.5f,%d" %(c_count+i, id1, id2, dist, int(count))

def main(simmat_path="data/first_110.atg.txt"):
    #load similairy matrix
    mat = load_simmat(simmat_path)
    
    #convert to distance matrix
    dist_mat = sim_to_dist(mat)
    
    #the idx to label function
    llf = lambda i: dist_mat.get_label(i)
    
    link_mat = linkage(dist_mat)
    
    disp_linakge_mat(link_mat, llf, len(dist_mat))

    grps = get_groups(link_mat, llf)

    for g in grps:
        print ",".join(g)

def test():
    import doctest
    doctest.testmod()
    
        
if __name__  == "__main__":
    #test()
    main("data/fp_370_atb.txt")
    
