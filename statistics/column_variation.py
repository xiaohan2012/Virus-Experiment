import os
from pickle import load
from itertools import chain

import numpy as np
import pylab as pl

from dist_mat_stat import *
from util.aa2code import *
from sim_mat import load_sim_mat
from rmsd_166_correlation import load_epi_mat

from config import *
def gen_fds(step):
    mat = load_epi_mat()
    avg = np.average(mat,0)

    sorted_avg = sorted(avg,reverse = True)
    top_threshold = sorted_avg[1]
    bottom_threshold = sorted_avg[-2]
    #print top_threshold,bottom_threshold
    #print sorted_avg

    top_col_ind = np.where(avg >= top_threshold)[0]
    bottom_col_ind = np.where(avg <= bottom_threshold)[0]

    codemap = gen_166_codemap()
    inv_codemap = gen_166_inv_codemap()
    ind = inv_codemap["1FBI_X"]

    for ind in chain(top_col_ind, bottom_col_ind, np.array([ind])):
        col = np.sort(mat[:,ind])
        yield codemap[ind],get_freq_dist(col,step = step)
   
def plot_and_save(fd, step, path = "/home/xiaohan/Desktop/col_var"):
    cfd.save_fig(title = codemap[ind] , save_path= path, step= step)

def save_data2js(fds, path):
    with open(path,"w") as f:
        names = []
        for name, fd in fds:
            items = [[str(k),v] for k,v in fd.sequential_items()]
            f.write("var $%s= %s;\n" %(name,repr(items) ))
            names.append(name)
        f.write("\n")
        #join data into one chunk
        f.write("var data= [%s];\n" %','.join(map(lambda s:"$"+ s, names)))
        
        #write column name
        f.write("var names= [%s];\n" %','.join(map(lambda s:"'%s'" %s, names)))
    
    #write other stuff

if __name__ == "__main__":
    fds= gen_fds(0.1)
    path = os.path.join(base,"js_plot/rmsd/data.js")
    save_data2js(fds, path)
