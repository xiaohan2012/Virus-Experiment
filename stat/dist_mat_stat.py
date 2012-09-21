from nltk import FreqDist
from sim_mat import load_sim_mat

def get_freq_dist(matrix , step=0.05):
    for r in matrix:
        print r

if __name__ == "__main__":
    mat = load_sim_mat("WILM950103_dist_mat" )
    get_freq_dist(mat)
