import glob
import os
from pickle import load

from util.aa2code import get_codes_from_file
from config import data_src , sim_dist_pickle_dir , mat_csv_path

def write_to_file(code_map , mat , csv_path):
    f = open(csv_path,'w')
    f.write(",%s\n" %",".join(code_map.values()))

    for i1,code1 in code_map.items():
        f.write("%s,%s\n" %(code1,','.join(map(str,mat[i1,:]))))

    f.close()


if __name__ == "__main__":
    code_map = get_codes_from_file(data_src)

    for pick_path in glob.glob(os.path.join(sim_dist_pickle_dir , '*')):
        mat = load(open(pick_path,'r'))
        mat = mat / mat.diagonal()
        path = os.path.join(mat_csv_path , os.path.basename(pick_path).split('.')[0] + '.csv')
        print path
        write_to_file(code_map,mat , path)
