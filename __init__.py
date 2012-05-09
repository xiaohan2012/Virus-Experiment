__author__ = 'xiaohan'
from schrodinger import structure
import fp_gen1
import avg_sift
import schrodinger.utils.fileutils as fileutils
import os
import glob
import re
def start_with_atom(line):
    """
    tell if the line is atom line
    """
    return line.split()[0].lower() == 'atom'

def is_hetatm(line):
    return line.split()[0].lower() == 'hetatm' or re.match("hetatm",line.lower())

def belongs_to_antibody(line):
    """
    tell if the line belongs to antibody
    """
    try:
        #return line.split()[4].upper() in 'HLIMJN'
        return line.split()[4].upper() in 'ABCD'
    except:
        return False

def split_big_mol_sml_mol(complex_path='',binder_path='',receptor_path=''):
    with open(receptor_path,'w') as rec_f:
        with open(binder_path,'w') as bnd_f:
            with open(complex_path) as f:
                for line in f.readlines():
                    if start_with_atom(line):
                        rec_f.write(line[:60]+'\n')
                    elif is_hetatm(line):
                        bnd_f.write(line[:60]+'\n')


def split_protein_protein_complex_manual(complex_path='',binder_path='',receptor_path=''):
    """
    @params:complex structure file path,to be saved antigen path,to be saved antibody path
    @return antibody structure, antigen structure
    """

    with open(receptor_path,'w') as ab_f:
        with open(binder_path,'w') as ag_f:
            with open(complex_path) as f:
                for line in f.readlines():
                    if start_with_atom(line):
                        if belongs_to_antibody(line):
                            ab_f.write(line[:60]+'\n')
                        else:
                            ag_f.write(line[:60]+'\n')


def gen_protein_protein_complex_avg_sift(complex_st_path,processed_data_path = 'processed_data'):
    complex_id,  ext = fileutils.splitext(os.path.basename(complex_st_path))

    sift_path = '%s/%s/%s_pattern.dat' %(processed_data_path,complex_id,complex_id)

    if os.path.exists(sift_path):
        print complex_id,'is processed'
        return


    antibody, antigen = split_protein_protein_complex_manual(complex_st_path,processed_data_path)#load complex structure and split it

    fp_gen_path = fp_gen1.gen_fp(receptor=antibody,binder = antigen,complex_id = complex_id,root_path = processed_data_path)#get the finger print
    avg_sift.gen_avg_sift(fp_gen_path,sift_path)#generate average sift


def cal_avg_sift_from_complex(complex_path='data/HL_chain/1RD8-1918/complex.1000.pdb',output_dir="/home/xiaohan/Desktop",split_fun=globals()['split_protein_protein_complex_manual']):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir);
    binder_path = os.path.join(output_dir,'binder.pdb')
    receptor_path = os.path.join(output_dir,'receptor.pdb')
    print binder_path,receptor_path 

    split_fun(complex_path = complex_path,binder_path = binder_path ,receptor_path=receptor_path)

    fp_path = os.path.join(output_dir,'fp.out')
    fp_gen_path = fp_gen1.gen_fp(receptor_file=receptor_path,binder_file = binder_path,fp_path= fp_path)#get the finger print
    #fp_gen_path = fp_gen1.gen_fp(receptor_file = binder_path,binder_file = receptor_path,fp_path= fp_path)#get the finger print

    sift_path=os.path.join(output_dir,'avg_sift.out')
    print sift_path
    avg_sift.gen_avg_sift(fp_gen_path,sift_path)#generate average sift

if __name__ == '__main__':
    data_src = '/home/xiaohan/Desktop/data/*'
    output_base_dir = '/home/xiaohan/Desktop/result'
    #data_src = 'data/HL_chain/HL_Q464S3/*'
    for fname in glob.glob(data_src):
        print fname
        complex_id = fname.split('/')[-1].split('.')[0]
        output_dir = os.path.join(output_base_dir,complex_id)
        if os.path.exists(os.path.join(output_dir,"avg_sift.out")):
            print "%s is processed" %complex_id
            continue
        cal_avg_sift_from_complex(fname,output_dir,split_fun=globals()['split_big_mol_sml_mol'])
