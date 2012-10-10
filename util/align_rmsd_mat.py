import numpy as np

from aa2code import *

def align(from_mat,to_mat_shape,from_inv_cm,to_cm):
    rc,cc = to_mat_shape
    tmp_mat = np.zeros(to_mat_shape)

    for r in xrange(rc):
        for c in xrange(cc):
            #print r,c
            r_chain_id = to_cm[r]
            c_chain_id = to_cm[c]

            #print r_chain_id,c_chain_id
            try:
                r_ind = from_inv_cm[r_chain_id]
                c_ind = from_inv_cm[c_chain_id]
            except KeyError:
                tmp_mat[r,c] = -1
                continue
            
            #print r_ind,c_ind
            tmp_mat[r,c] = from_mat[r_ind,c_ind]
            #print
    return tmp_mat            

if __name__ == "__main__":
    from_inv_cm = get_rmsd_inv_codemap()
    to_cm = gen_166_codemap()
    mat = np.loadtxt(rmsd_matrix_path)
    mat1 = align(mat, (166,166) , from_inv_cm, to_cm)
   
    print mat1
    #test
    import random

    r = random.choice(range(0,166))
    c = random.choice(range(0,166))
    print mat1[r,c]
    
    r_chain = to_cm[r]
    c_chain = to_cm[c]
    
    r = from_inv_cm[r_chain]
    c = from_inv_cm[c_chain]
    print mat[r,c]
    print 

