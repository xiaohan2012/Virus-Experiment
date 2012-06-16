import glob
import os
import pymongo 
from numpy import array,zeros
from scipy.cluster.vq import vq, kmeans, whiten
from pickle import dump,load
from collections import defaultdict


def get_codes_from_file(path):
    codes = {}
    for i,pdb_path in enumerate(glob.glob(pdb_src)):
        code= os.path.basename(pdb_path).split('.')[0] 
        codes[i] = code
    
    return codes

def get_matrix(code_map):
    conn = pymongo.Connection()
    db = pymongo.database.Database(conn , "virus_cluster")
    col = db["epi_166_sim_mat"]
    size = len(code_map.keys())
    mat = zeros((size,size))

    for i1,code1 in code_map.items():
        for i2,code2 in code_map.items():
            if i1 == i2:
                continue
            res = col.find_one({"complex1":code1,"complex2":code2})
            if res:
                mat[i1,i2] = -res["similarity"]#inverse it
                continue

            res = col.find_one({"complex1":code2,"complex2":code1})
            if res:
                mat[i1,i2] = -res["similarity"]#inverse it
                continue
    
    return mat

def write_to_file(code_map,mat):
    f = open("dist_mat.csv",'w')
    f.write(",%s\n" %",".join(code_map.values()))

    for i1,code1 in code_map.items():
        f.write("%s,%s\n" %(code1,','.join(map(str,mat[i1,:]))))

    f.close()

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
    obj_fp = "dist_mat.dat"
    pdb_src = "epi_166/pdb_file/*"
    code_map = get_codes_from_file(pdb_src)
    try:
        mat = load(open(obj_fp,'r'))
    except IOError:
        mat = get_matrix(code_map)
        dump(mat,open(obj_fp,'w'))

    for k in xrange(3,10):
        observations = whiten(mat)
        centroids,distortions = kmeans(observations , k )

        assignments = assign(centroids , observations , code_map)

        with open("kmeans_result/%d_clusters.dat" %k , 'w') as f:
            for a in assignments.values():
                f.write("%s\n" %(' '.join(a)))

    write_to_file(code_map,mat)
