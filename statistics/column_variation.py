from pickle import load
import os

import numpy as np
import pylab 

from config import *

hydro_id = "ARGP820101_dist_mat.mat"
mat = load(open(os.path.join(sim_dist_pickle_dir , hydro_id)))
avg = np.average(mat,1)

sorted_avg = sorted(avg,reverse = True)
top_threshold = sorted_avg[2]
bottom_threshold = sorted_avg[-3]
top_col = np.where(avg >= top_threshold)[0]
bottom_col = np.where(avg <= bottom_threshold)[0]


print mat[:,top_col]
print mat[:,bottom_col]
