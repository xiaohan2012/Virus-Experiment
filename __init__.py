__author__ = 'xiaohan'
from schrodinger import structure
import fp_gen1
import avg_sift
import schrodinger.utils.fileutils as fileutils
import os
import glob

def start_with_atom(line):
    """
    tell if the line is atom line
    """
    return line.split()[0].lower() == 'atom'

def belongs_to_antibody(line):
    """
    tell if the line belongs to antibody
    """
    try:
        return line.split()[4].upper() in 'HLIMJN'
    except:
        return False


def split_protein_protein_complex_manual(complex_path='',antigen_path='',antibody_path=''):
    """
    @params:complex structure file path,to be saved antigen path,to be saved antibody path
    @return antibody structure, antigen structure
    """

    with open(antibody_path,'w') as ab_f:
        with open(antigen_path,'w') as ag_f:
            with open(complex_path) as f:
                for line in f.readlines():
                    if start_with_atom(line):
                        if belongs_to_antibody(line):
                            ab_f.write(line[:60]+'\n')
                        else:
                            ag_f.write(line[:60]+'\n')

    antibody_id,_ = fileutils.splitext(os.path.basename(antibody_path))
    antigen_id,_ = fileutils.splitext(os.path.basename(antigen_path))

    antibody = structure.StructureReader(antibody_path).next()
    antibody.title = antibody_id

    antigen = structure.StructureReader(antigen_path).next()
    antigen.title = antigen_id

    return antibody,antigen

def gen_protein_protein_complex_avg_sift(complex_st_path,processed_data_path = 'processed_data'):
    complex_id,  ext = fileutils.splitext(os.path.basename(complex_st_path))

    sift_path = '%s/%s/%s_pattern.dat' %(processed_data_path,complex_id,complex_id)

    if os.path.exists(sift_path):
        print complex_id,'is processed'
        return


    antibody, antigen = split_protein_protein_complex_manual(complex_st_path,processed_data_path)#load complex structure and split it

    fp_gen_path = fp_gen1.gen_fp(receptor=antibody,binder = antigen,complex_id = complex_id,root_path = processed_data_path)#get the finger print
    avg_sift.gen_avg_sift(fp_gen_path,sift_path)#generate average sift

def cal_avg_sift_from_complex(complex_path='data/HL_chain/1RD8-1918/complex.1000.pdb'):
    print complex_path
    _,chain_name,complex_id,instance_name = complex_path.split('/')

    p_chain_path=os.path.join('processed_data',chain_name)
    if not os.path.exists(p_chain_path):
        os.mkdir(p_chain_path)

    p_complex_path=os.path.join(p_chain_path,complex_id)
    if not os.path.exists(p_complex_path):
        os.mkdir(p_complex_path)

    p_instance_path=os.path.join(p_complex_path,instance_name)[:-4]
    if not os.path.exists(p_instance_path):
        os.mkdir(p_instance_path)

    antigen_path = os.path.join(p_instance_path,'antigen.pdb')
    antibody_path = os.path.join(p_instance_path,'antibody.pdb')
    print antigen_path,antibody_path 
    antibody,antigen = split_protein_protein_complex_manual(complex_path = complex_path,antigen_path = antigen_path ,antibody_path=antibody_path)

    fp_path = os.path.join(p_instance_path,'fp.out')
    print antibody,antigen
    try:
        fp_gen_path = fp_gen1.gen_fp(receptor=antibody,binder = antigen,fp_path= fp_path)#get the finger print
    except:
        print 'failed\n'
        return
    else:
        print 'good\n'

    sift_path=os.path.join(p_instance_path,'avg_sift.out')
    print sift_path
    avg_sift.gen_avg_sift(fp_gen_path,sift_path)#generate average sift
    #gen_protein_protein_complex_avg_sift('1918_complex.pdb')

if __name__ == '__main__':
    for fname in glob.glob('data/HL_chain/HL_Q464S3/*'):
        cal_avg_sift_from_complex(fname)
