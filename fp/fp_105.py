import sys
import os

from itertools import chain
from glob import glob

from fp_75_gen import get_15bits

from residue_util import *
from fp import BaseResidueFingerprint,BaseComplexFingerprint
from complex import TriangleComplex

from ve.util.load_pdb import load_pdb_struct
from ve.util.file import fpstr2file
from ve.util.residue import BaseResidue
from ve.util.dist import FP105DistCache

from ve.config import *

class ElectricResidue(BaseResidue):
    def __init__(self,residue):
        BaseResidue.__init__(self,residue)
        self._set_ca_atom()

        self.dist_min = 0

        self.fp = None
    
        #adding the required methods
        init_resdist_util(self)
        init_neighbour_util(self)
        
        #add electric finger print calculation function
        init_electric_fp_util(self)

        self.bitlength = 30
        
    def dist2residue(self,other):
        return self.distcache.get(self,other)

    def find_neighbours(self,others,dist=4):
        """find neighbours"""
        for other in others:
            if  self.dist2residue(other) <= dist and self != other:
                self.nbs.append(other)
        return self.nbs

    def _set_ca_atom(self):
        self.ca = None
        for a in self.atom:
            if a.pdbname.strip() == "CA":
                self.ca = a
                break
        if self.ca is None:
            sys.stderr.write("resnum %d\n" %self.resnum)
            raise ValueError("self.ca should not be none")

    def gen_fp(self,others):
        self.gen_electric_fp(others)
        return self.fp

    def gen_fp_triangle(self,others):
        self.gen_electric_fp_triangle(others)
        return self.fp



        


class Complex(TriangleComplex):
    def __init__(self,antigen,antibody):
        TriangleComplex.__init__(self,antigen,antibody) 
        self.bitlength = 105
        self.fp = BaseComplexFingerprint()

        self.distcache = FP105DistCache()

    def gen_15_bits(self):
        self.atg_fp = get_15bits(self.atg, self.atb)
        self.atb_fp = get_15bits(self.atb, self.atg)

    def gen_fps(self):
        #gen the first 15 bits
        self.gen_15_bits()

        for residue in self.atg.residues:
            others = self.atg.residues#1~30, iterate antigen side
            fp1 = residue.gen_fp(others)

            others = self.atb.residues#31~60, iterate antibody side
            fp2 = residue.gen_fp(others)

            others = self.triangles#61~90. iterate triangles
            fp3 = residue.gen_fp_triangle(others)
            
            fp1.append(fp2)
            fp1.append(fp3)
            
            #append the following bits to the first 15 bits
            self.atg_fp[residue].append(fp1)

        for residue in self.atb.residues:
            others = self.atg.residues#1~30, iterate antigen side
            fp1 = residue.gen_fp(others)

            others = self.atb.residues#31~60, iterate antibody side
            fp2 = residue.gen_fp(others)

            others = self.triangles#61~90. iterate triangles
            fp3 = residue.gen_fp_triangle(others)
            
            fp1.append(fp2)
            fp1.append(fp3)
            
            #append the following bits to the first 15 bits
            self.atb_fp[residue].append(fp1)

        return self.atg_fp,self.atb_fp

    def write_structure(self,i_path,o_path,resnums):
        with open(i_path,"r") as input, open(o_path,"w") as output:
            for l in input.readlines():
                resnum = ''.join(l[22:26]).strip()
                if resnum in resnums:output.write(l)

    def write_epitope(self,source_fp,output_fp):
        resnums = map(lambda a:str(a.resnum),self.epitope)
        i_path = os.path.join(source_fp,"antigen.pdb")
        o_path = os.path.join(output_fp,"epitope.pdb")
        
        self.write_structure(i_path,o_path,resnums)

    def write_paratope(self,source_fp,output_fp):
        resnums = map(lambda a:str(a.resnum),self.paratope)
        i_path = os.path.join(source_fp,"antibody.pdb")
        o_path = os.path.join(output_fp,"paratope.pdb")
        
        self.write_structure(i_path,o_path,resnums)


####Real life date processing part#####
def data237_fp_gen(refresh = True):
    path = os.path.join(data237_complex_root,"*")
    print path
    for fp in glob(path):

        complex_id = os.path.basename(fp)
        #lst = ["1SLG_D","1T6V_L","1N4X_L","1JV5_A"]
        #lst = ["1SLG_D"]
        #if complex_id not in lst:
            #continue
        
        #prepare neccessary directory
        if not os.path.exists(data237_fp105_root): os.makedirs(data237_fp105_root)

        fp_path = os.path.join(data237_fp105_root,complex_id)
        
        if not refresh and os.path.exists(fp_path):
            print "%s processed" %complex_id
        else:
            print "processing %s" %complex_id
            atb_path = os.path.join(fp,"antibody.pdb")
            atg_path = os.path.join(fp,"antigen.pdb")
            try:
                antibody = load_pdb_struct(atb_path,residue_cls = ElectricResidue)
                antigen = load_pdb_struct(atg_path,residue_cls = ElectricResidue)
                complex = Complex(antigen,antibody)

                out_path =os.path.join(data237_fp105_root,complex_id) 

                if not os.path.exists(out_path):os.makedirs(out_path)

                atg_fp,atb_fp = complex.gen_fps()
                atg_fp.tofile(os.path.join(out_path,"antigen.csv"))
                atb_fp.tofile(os.path.join(out_path,"antibody.csv"))


            except:
                sys.stderr.write("%s encountered error\n" %complex_id)
                continue
    
def fp_gen_test(complex_id):
    print "procssing %s" %complex_id

    path = os.path.join(data237_complex_root,complex_id)
    print path

    atb_path = os.path.join(path,"antibody.pdb")
    atg_path = os.path.join(path,"antigen.pdb")
    antibody = load_pdb_struct(atb_path,residue_cls = ElectricResidue)
    antigen = load_pdb_struct(atg_path,residue_cls = ElectricResidue)

    complex = Complex(antigen,antibody)

    fp_path = os.path.join(data237_fp105_root,"%s.csv" %complex_id)
    fp = complex.gen_fp()
    fp.tofile(fp_path)

def data237_epitope_paratope(refresh = False):
    source_path = os.path.join(data237_complex_root,"*")
    print source_path


    for fp in glob(source_path):
        complex_id = os.path.basename(fp)

        if not refresh and os.path.exists(os.path.join(data237_epitope_root,complex_id,"paratope.pdb")):
            print "%s processed" %complex_id
        else:
            print "processing %s" %complex_id
            try:
                atb_path = os.path.join(fp,"antibody.pdb")
                atg_path = os.path.join(fp,"antigen.pdb")
                antibody = load_pdb_struct(atb_path,residue_cls = ElectricResidue)
                antigen = load_pdb_struct(atg_path,residue_cls = ElectricResidue)

                complex = Complex(antigen,antibody)
                output_path = os.path.join(data237_epitope_root, complex_id)

                if not os.path.exists(output_path):
                    os.makedirs(output_path)

                complex.write_epitope(fp,output_path)
                complex.write_paratope(fp,output_path)
                print "epitope and paratope saved"
            except:
                sys.stderr.write("%s encountered error\n" %complex_id)
                continue



def epi_gen_test(complex_id):
    print "procssing %s" %complex_id

    path = os.path.join(data237_complex_root,complex_id)
    print path
    output_path = os.path.join(data237_epitope_root, complex_id)

    atb_path = os.path.join(path,"antibody.pdb")
    atg_path = os.path.join(path,"antigen.pdb")
    antibody = load_pdb_struct(atb_path,residue_cls = ElectricResidue)
    antigen = load_pdb_struct(atg_path,residue_cls = ElectricResidue)

    complex = Complex(antigen,antibody)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    complex.write_epitope(path,output_path)
    complex.write_paratope(path,output_path)

    

if __name__ == "__main__":
    #data237_fp_gen(refresh = False)
    data237_fp_gen(False)
    #data237_epitope_paratope(False)
    #fp_gen_test("1SLG_D")
