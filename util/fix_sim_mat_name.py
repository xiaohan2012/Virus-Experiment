import os
import re
from glob import glob


from config import sim_dist_pickle_dir 
def fix_mat_name():
    for fp in glob(os.path.join(sim_dist_pickle_dir,"*")):
        mat_id = re.findall(r"[A-Z]{3,}\d{3,}",fp)[0]
        new_fp = os.path.join(os.path.dirname(fp), mat_id + ".mat")
        os.rename(fp,new_fp)

if __name__ == "__main__":
    fix_mat_name()

