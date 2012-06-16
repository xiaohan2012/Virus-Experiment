from clustering import *
import os
import glob
import pymongo 

if __name__ == "__main__":
    conn = pymongo.Connection()
    db = pymongo.database.Database(conn , "virus_cluster")

    pdb_src = "epi_166/pdb_file/*"
    fp_dir = "epi_166/fp_result"
    count = -1
    for pdb_fp1 in glob.glob(pdb_src):
        count += 1
        print count
        for pdb_fp2 in glob.glob(pdb_src):
            if pdb_fp1 == pdb_fp2:
                continue
            complex1_id = os.path.basename(pdb_fp1).split('.')[0] 
            fp1_path = os.path.join(fp_dir,complex1_id + ".fp" )

            complex2_id = os.path.basename(pdb_fp2).split('.')[0] 
            fp2_path = os.path.join(fp_dir,complex2_id + ".fp" )

            print complex1_id,complex2_id

            if db["epi_166_sim_mat"].find_one({"complex1":complex1_id,"complex2":complex2_id}) or\
               db["epi_166_sim_mat"].find_one({"complex2":complex1_id,"complex1":complex2_id}):
                print "%s and %s already calculated" %(complex1_id,complex2_id)
            else:                
                #save the similarity value
                pair = dist_mat(pdb_fp1 , pdb_fp2 , fp1_path , fp2_path)

                while pair.find_cluster():
                    pass
                sim_val = get_similarity(pair.clusters,pair)
                db["epi_166_sim_mat"].save({"complex1":complex1_id,"complex2":complex2_id,"similarity":sim_val})



