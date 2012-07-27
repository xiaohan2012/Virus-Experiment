import pymongo

host = "bnode1"
dbname = "virus_cluster"

fp_root = '/home/rxzhu/code/Virus-Experiment/402/fp_result'
data_src = "/home/rxzhu/code/Virus-Experiment/402/pdb_file/*"
output_dir = "/home/rxzhu/code/Virus-Experiment/402/fp_result"

db = pymongo.database.Database(pymongo.Connection(host), dbname)

