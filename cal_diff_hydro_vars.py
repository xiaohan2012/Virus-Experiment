import os,  sys
import glob
import pymongo

from fp_gen2 import Complex
from fp2_clustering import gen_sim_matrix
from logger import TaskFileLogger

h_fp = "hydro_variations.dat"
def load_hydro_var(h_fp = "hydro_variations.dat"):
    with open(h_fp,'r') as f:
        reses = f.readline().split()
        d_ = {}
        for l in f.readlines():
            linkdb = l.split()[0]
            d_[linkdb] = {}
            for res,num in zip(reses,l.split()[1:]):
                d_[linkdb][res] = float(num)
    return d_

def gen_fps():
    """generate finger print files"""
    logger = TaskFileLogger("GenFP")

    h_vars = load_hydro_var()
    data_src = "402/pdb_file/*"
    output_dir = "402/fp_result"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    for h_type,var_d in h_vars.items():
        print "considering %s" %h_type

        t_output_dir = os.path.join(output_dir,h_type)
        if not os.path.exists(t_output_dir):
            print "creating path %s" %t_output_dir
            os.mkdir(t_output_dir)
        logger.log("%s started" %(h_type))

        for fname in glob.glob(data_src):
            complex_id = os.path.basename(fname).split('.')[0] 
            fp_path = os.path.join(t_output_dir,complex_id + ".fp" )
            if os.path.exists(fp_path):
                print "%s processed" %complex_id
                continue
            print "processing %s,fp saved as %s" %(fname , fp_path )
            c = Complex(fname,hydro_dict = var_d)
            c.get_fp()
            c.write_fp_to_file(fp_path)

        logger.log("%s finished" %(h_type))

def gen_sim_mats():
    """generate the distance matrix"""
    logger = TaskFileLogger("GenMat")

    h_vars = load_hydro_var()
    h_names = h_vars.keys()
    fp_root = '402/fp_result'
    
    conn = pymongo.Connection()
    dbname = "virus_cluster"
    db = pymongo.database.Database(conn , dbname)

    for h_name in h_names:

        logger.log("%s started" %(h_name))
        gen_sim_matrix(db,\
                        col_name = "%s_dist_mat" %h_name,\
                        fp_dir = os.path.join(fp_root,h_name))
        logger.log("%s finished" %(h_name))

if __name__ == "__main__":
    gen_fps()
    gen_sim_mats()
