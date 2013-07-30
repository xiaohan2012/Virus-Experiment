"""
calculating the similarity matrix from finger print files
and the save the matrix in database
"""

from similarity import *
from util.time_measure import tic,toc
from config import *

import os
import glob
import pymongo 
import bson
import time
from threading import Thread
from Queue import Queue


def gen_sim_matrix(db,col_name = "402_another" ,pdb_src = "/home/rxzhu/code/Virus-Experiment/402/pdb_file/*" , fp_dir = "/home/rxzhu/code/Virus-Experiment/402/fp_result"):
    db[col_name].ensure_index([("complex1",pymongo.ASCENDING),("complex2",pymongo.ASCENDING)],unique = True)
    db[col_name].ensure_index([("complex2",pymongo.ASCENDING),("complex1",pymongo.ASCENDING)],unique = True)
    count = -1
    for pdb_fp1 in glob.glob(pdb_src):
        count += 1
        print count

        complex1_id = os.path.basename(pdb_fp1).split('.')[0] 
        fp1_path = os.path.join(fp_dir,complex1_id + ".fp" )
        
        c1 = FPWithComplex(pdb_fp1, fp1_path)
        
        for pdb_fp2 in glob.glob(pdb_src):

            complex2_id = os.path.basename(pdb_fp2).split('.')[0] 
            fp2_path = os.path.join(fp_dir,complex2_id + ".fp" )

            if db[col_name].find_one({"complex1":complex1_id,"complex2":complex2_id}) or\
               db[col_name].find_one({"complex2":complex1_id,"complex1":complex2_id}):
                print "%s and %s already calculated" %(complex1_id,complex2_id)
            else:                

                c2 = FPWithComplex(pdb_fp2, fp2_path)
                
                val1, val2, val3 = similarity_between(c1, c2)
                
                db[col_name].save({"complex1":complex1_id,"complex2":complex2_id,"val1":val1 ,"val2":val2 ,"val3":val3 })
    return True

def generate_tasks(pdb_src):
    for pdb_fp1 in glob.glob(pdb_src):
        for pdb_fp2 in glob.glob(pdb_src):
            complex1_id = os.path.basename(pdb_fp1).split('.')[0] 
            complex2_id = os.path.basename(pdb_fp2).split('.')[0] 
            yield complex1_id , complex2_id
    
class CalPairThread(Thread):
    def __init__(self , queue , db , col_name = None , pdb_src = "/home/rxzhu/code/Virus-Experiment/402/pdb_file/*" , fp_dir = "/home/rxzhu/code/Virus-Experiment/402/fp_result"):
        Thread.__init__(self)
        self.q = queue
        self.db = db
        self.col_name = col_name
        self.pdb_src =  pdb_src 
        self.fp_dir = fp_dir

    def run(self):
        
        while self.q.qsize():
            complex1_id , complex2_id = self.q.get()#get complex id
            #reassemble pdb file path and finger print path of complex 1
            pdb1_fp = os.path.join( os.path.dirname(self.pdb_src) , complex1_id + ".pdb" )
            fp1_path = os.path.join(fp_dir,complex1_id + ".fp" )

            #reassemble pdb file path and finger print path of complex 2
            pdb2_fp = os.path.join( os.path.dirname(self.pdb_src) , complex2_id + ".pdb" )
            fp2_path = os.path.join(fp_dir,complex2_id + ".fp" )
            
            #print pdb1_fp , fp1_path , pdb2_fp , fp2_path 
            
            #prevent repeated calculations
            if not db[col_name].find_one({"complex1":complex1_id,"complex2":complex2_id}) and \
               not db[col_name].find_one({"complex2":complex1_id,"complex1":complex2_id}):
                pair = dist_mat(pdb1_fp , pdb2_fp , fp1_path , fp2_path)

                while pair.find_cluster():
                    pass
                val1 , val2 , val3 = get_similarity(pair.clusters,pair)
                db[col_name].save({"complex1":complex1_id,"complex2":complex2_id,"val1":val1 ,"val2":val2 ,"val3":val3 })
                

            #notify task is done
            print complex1_id,complex2_id,"done","%d tasks left" %self.q.qsize()
            self.q.task_done()


if __name__ == "__main__":
    """
    thread_count = 5#thread count
    col_name = "thread_level_parallel_test"#collection name 
    fp_dir = os.path.join(fp_root ,"ARGP820101")#finger print directory

    queue = Queue()

    for t in generate_tasks(data_src):
        queue.put(t)

    threads = []

    db["log"].save({"node":'bnode1',"type":"started","time":time.time(),"task_name":'thread level parallel test'})
    for i in xrange(thread_count):
        threads.append( CalPairThread(queue , db , col_name , data_src , fp_dir) )
    for t in threads:
        t.setDaemon(True)#don't forget this
        t.start()

    queue.join()        
    print "task done"
    db["log"].save({"node":'bnode1',"type":"finished","time":time.time(),"task_name":'thread level parallel test'})
    # gen_sim_matrix(db)
    """
    from ve.config import pdb_166_src, epi166_fp
    from ve.dbconfig import db
    
    gen_sim_matrix(db,col_name = "epi_166" ,pdb_src = pdb_166_src , fp_dir = epi166_fp)