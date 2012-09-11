"""
using k-means clustering to classify protein complexes into clusters
"""
import glob
import os
import pymongo 
from numpy import array,zeros
from scipy.cluster.vq import vq, kmeans, whiten
from pickle import dump,load
from collections import defaultdict

from util import *
from sim_mat import get_sim_matrix_from_db


def assign(centroids , observations , code_map):
    assignment = defaultdict(list)
    for ob_i , ob in enumerate(observations):
        min_dist = float("inf")
        closest_cent = None
        for cent_i , cent in enumerate(centroids):
            dist = sum(((cent - ob)) ** 2)**0.5
            if dist < min_dist:
                min_dist = dist
                closest_cent = cent_i
        #print closest_cent , code_map[ob_i]
        assignment[closest_cent].append(code_map[ob_i])
    return assignment

if __name__ == "__main__":
    obj_fp = "dist_mat_epi_166.dat"
    pdb_src = "epi_166/pdb_file/*"
    code_map = get_codes_from_file(pdb_src)
    print len(code_map.keys())
    try:
        mat = load(open(obj_fp,'r'))
    except IOError:
        mat = get_sim_matrix_from_db(code_map)
        dump(mat,open(obj_fp,'w'))

    for k in xrange(3,10):
        observations = whiten(mat)
        centroids,distortions = kmeans(observations , k )

        assignments = assign(centroids , observations , code_map)

        with open("epi_166/kmeans_result/%d_clusters.dat" %k , 'w') as f:
            for a in assignments.values():
                f.write("%s\n" %(' '.join(a)))
    mat = mat / mat.diagonal()
    write_to_file(code_map,mat)
