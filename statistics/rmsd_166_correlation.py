import numpy as np
import pylab

from util.aa2code import *
from util.align_rmsd_mat import align
from config import *

def get_aligned_rmsd_mat():
    from_inv_cm = get_rmsd_inv_codemap()
    to_cm = gen_166_codemap()
    mat = np.loadtxt(rmsd_matrix_path)
    return align(mat, (166,166) , from_inv_cm, to_cm)

def load_epi_mat():
    return np.loadtxt(dist_mat_path )

def gen_column_data(chain_id):
    inv_map = gen_166_inv_codemap()
    #print inv_map.keys()                
    c_ind = inv_map[chain_id]
    
    rmsd_mat = get_aligned_rmsd_mat()
    epi_mat = load_epi_mat()
    
    rmsd_col = rmsd_mat[:,c_ind]
    epi_col = epi_mat[:,c_ind]
    
    #filter out the -1 rows(marked in the alignment step)
    non_zero_inds = np.where(rmsd_col>=0)
    rmsd_col = rmsd_col[non_zero_inds]#/ np.max(rmsd_col)
    epi_col = epi_col[non_zero_inds]
    
    if len(rmsd_col) and len(epi_col):
        kb = np.polyfit(rmsd_col,epi_col,1)
        k,b = kb

        cf_x = np.arange(np.floor(np.min(rmsd_col)),np.ceil(np.max(rmsd_col)),\
                (np.max(rmsd_col) - np.min(rmsd_col))/10)
        cf_y = k * cf_x + b
        return  chain_id, rmsd_col, epi_col, cf_x, cf_y 

def compare(chain_id):
    result = gen_column_data(chain_id)
    if result:
        rmsd_col, epi_col, cf_x, cf_y = result

        pylab.scatter(rmsd_col , epi_col)
        pylab.hold(True)
        pylab.plot(cf_x,cf_y,"-r")

        pylab.xlabel("RMSD vector")
        pylab.ylabel("Our similarity vector")
        pylab.title("Curve fitting result for %s column of RMSD matrix and similarity matrix" %chain_id)

        pylab.xlim(np.floor(np.min(rmsd_col)),np.ceil(np.max(rmsd_col)))
        pylab.ylim(0,np.max(epi_col))

        pylab.savefig("/home/xiaohan/Desktop/rmsd_166_comp/%s.png" %chain_id)
        pylab.hold(False)

def gen_all_data():
    rmsd_vec = get_aligned_rmsd_mat().flatten()
    epi_vec = load_epi_mat().flatten()

    #filter out the -1 rows(marked in the alignment step)
    non_zero_inds = np.where(rmsd_vec>=0)
    rmsd_vec = rmsd_vec[non_zero_inds] #/ np.max(rmsd_vec)
    epi_vec = epi_vec[non_zero_inds]

    kb = np.polyfit(rmsd_vec,epi_vec,1)
    k,b = kb

    cf_x = np.arange(0,np.ceil(np.max(rmsd_vec)),0.5)
    cf_y = k * cf_x + b
    
    return rmsd_vec, epi_vec, cf_x, cf_y 

def compare_all():
    rmsd_vec, epi_vec, cf_x, cf_y = gen_all_data()
    pylab.scatter(rmsd_vec , epi_vec)
    pylab.hold(True)
    pylab.plot(cf_x,cf_y,"-r")
    pylab.xlabel("RMSD matrix")
    pylab.ylabel("Our similarity matrix")
    pylab.title("Curve fitting result for all items in RMSD and our similarity matrix")

    pylab.savefig("/home/xiaohan/Desktop/rmsd_166_comp/all.png")
    pylab.hold(False)

def gen_above_line_data(start,end):
    x,y = start
    diff_x , diff_y = start - end
    k = diff_y / diff_x
    b = y - k*x
    
    print k,b
    rmsd_vec = get_aligned_rmsd_mat().flatten()
    epi_vec = load_epi_mat().flatten()

    ind = np.where(rmsd_vec*k  + b - epi_vec <= 0)

    rmsd_vec = rmsd_vec[ind]
    epi_vec = epi_vec[ind]

    kb = np.polyfit(rmsd_vec,epi_vec,1)
    k,b = kb
    cf_x = np.arange(0,np.ceil(np.max(rmsd_vec)),0.5)
    cf_y = k * cf_x + b

    return rmsd_vec , epi_vec , cf_x , cf_y

def compare_above_line_kb(start,end):
    rmsd_vec , epi_vec , cf_x , cf_y = gen_above_line_data(start,end)

    pylab.scatter(rmsd_vec , epi_vec)
    pylab.hold(True)

    pylab.plot(cf_x,cf_y,"-r")
    
    pylab.ylim(0,1)
    pylab.xlim(0,16)
    pylab.ylabel("Our similarity matrix")
    pylab.xlabel("RMSD matrix")
    pylab.title("Curving fitting for points above line y = -0.02*x + 0.6")
    pylab.savefig("/home/xiaohan/Desktop/rmsd_166_comp/above_line.png")

    pylab.hold(False)

def gen_js_data(columns, start, end, path):
    column_data = [gen_column_data(column) for column in columns]
    all_data = gen_all_data()
    above_line_data = gen_above_line_data(start, end)
    
    with open(path, "w") as f:
        ###column data####
        f.write("var cols = {\n");
        for name, rmsd, epi, cf_x, cf_y in column_data:
            f.write("\t'%s':{\n" %name);
            f.write("\t\t'xy':%s,\n" %repr(map(list,zip(rmsd.tolist(),epi.tolist()))))
            #f.write("\t\t'y':%s,\n" %repr())
            f.write("\t\t'cf_xy':%s\n\t},\n" %repr(map(list,zip(cf_x.tolist(),cf_y.tolist()))))
            #f.write("\t\t'cf_y':%s};\n" %repr(cf_y.tolist()))
        f.write("};\n")
        
        ###all data###
        rmsd, epi, cf_x, cf_y = all_data
        f.write("var all = {\n");
        f.write("\t'xy':%s,\n" %(repr(map(list,zip(rmsd.tolist(),epi.tolist())))))
        f.write("\t'cf_xy':%s\n" %(repr(map(list,zip(cf_x.tolist(),cf_y.tolist())))))
        f.write("};\n")

        ###above line data###
        rmsd, epi, cf_x, cf_y = above_line_data 
        f.write("var above_line= {\n");
        f.write("\t'xy':%s,\n" %(repr(map(list,zip(rmsd.tolist(),epi.tolist())))))
        f.write("\t'cf_xy':%s\n" %(repr(map(list,zip(cf_x.tolist(),cf_y.tolist())))))
        f.write("};\n")

if __name__ == "__main__":
    """
    start = np.array((0,0.6))
    end = np.array((20,0.2))
    compare_above_line_kb(start,end)
    compare_all()
    chain_ids = gen_166_codemap().values()
    for cid in ["1BGX_T","1E6J_P","1FBI_X","2JEL_P","3LH2_S"] +["2P45_C", "1J1O_Y", "1G7J_C", "2P48_A", "1FE8_A"]:
        try:
            compare(cid)
        except KeyError:    
            pass
    """
    cols = ["1BGX_T","1E6J_P","1FBI_X","2JEL_P","3LH2_S"] +["1J1O_Y", "1G7J_C", "2P48_A", "1FE8_A"] + ["1MHH_E","1JPS_T"]#"2P45_C", 
    start = np.array((0,0.6))
    end = np.array((20,0.2))
    path = os.path.join(base,"js_plot/rmsd/data.js")

    gen_js_data(cols, start, end, path)
