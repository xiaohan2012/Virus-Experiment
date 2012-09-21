from nltk import FreqDist
from sim_mat import load_sim_mat
import math

def get_freq_dist(matrix , step=0.05):
    return FreqDist(math.floor(c/step) * step for r in matrix for c in r)
            


if __name__ == "__main__":
    mat = load_sim_mat("WILM950103_dist_mat" )
    fd = get_freq_dist(mat)
    fd.tabulate([0.0,0.25,0.5,0.75,1.0])
