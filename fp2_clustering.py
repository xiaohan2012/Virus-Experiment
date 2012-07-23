from clustering import *
import os
import glob
import pymongo 
import bson

def gen_sim_matrix(db,col_name = "402_another" ,pdb_src = "402/pdb_file/*" , fp_dir = "402/fp_result"):
    db[col_name].ensure_index([("complex1",pymongo.ASCENDING),("complex2",pymongo.ASCENDING)],unique = True)
    db[col_name].ensure_index([("complex2",pymongo.ASCENDING),("complex1",pymongo.ASCENDING)],unique = True)
    count = -1
    for pdb_fp1 in glob.glob(pdb_src):
        count += 1
        print count
        for pdb_fp2 in glob.glob(pdb_src):
            complex1_id = os.path.basename(pdb_fp1).split('.')[0] 
            fp1_path = os.path.join(fp_dir,complex1_id + ".fp" )

            complex2_id = os.path.basename(pdb_fp2).split('.')[0] 
            fp2_path = os.path.join(fp_dir,complex2_id + ".fp" )

            print complex1_id,complex2_id
            
            if db[col_name].find_one({"complex1":complex1_id,"complex2":complex2_id}) or\
               db[col_name].find_one({"complex2":complex1_id,"complex1":complex2_id}):
                print "%s and %s already calculated" %(complex1_id,complex2_id)
            else:                
                #save the similarity value
                pair = dist_mat(pdb_fp1 , pdb_fp2 , fp1_path , fp2_path)

                while pair.find_cluster():
                    pass
                val1 , val2 , val3 = get_similarity(pair.clusters,pair)
                db[col_name].save({"complex1":complex1_id,"complex2":complex2_id,"val1":val1 ,"val2":val2 ,"val3":val3 })
    return True
if __name__ == "__main__":
    conn = pymongo.Connection()
    dbname = "virus_cluster"
    db = pymongo.database.Database(conn , dbname)
    gen_dist_matrix(db)
