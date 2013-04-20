from ve.util.residue import BaseResidue

import sys
import logging

logging.basicConfig( stream=sys.stderr )
logging.getLogger("FP370").setLevel( logging.DEBUG )
logger = logging.getLogger("FP370")

class MyResidue(BaseResidue):
    def __init__(self, res):
        BaseResidue.__init__(self, res)
        
        from ve.fp.residue_util.res_geom import init_geom_trait
        init_geom_trait(self)
        
from ve.util.complex import BaseComplex        

class MyComplex(BaseComplex):
    def __init__(self,c_id, atg, atb):
        self.c_id = c_id
        self.atg = atg
        self.atb = atb
        
        #Split cylinde part finger print
        from ve.fp.complex_util.split_cylinder import init_gcb_split_cylinder_trait
        init_gcb_split_cylinder_trait(self)
        
        #physical chemical finger print
        from ve.fp.complex_util.physi_chemi import init_physi_chemi_fp_trait
        init_physi_chemi_fp_trait(self)
        
        #interactive force finger print
        from ve.fp.complex_util.interactive_force_fp import init_interative_force_fp_trait
        init_interative_force_fp_trait(self)
        
        
    def get_fp(self):
        logger.info("gen fp for %s", self.c_id)

        #antigen side split cylinder fingerprint
        fp1 = self.get_atg_fp_by_split_cylinder()

        #antigen side physical chemical fingerprint
        fp2 = self.get_physi_chemi_atg_fp()

        #antibody side split cylinder fingerprint
        fp3 = self.get_atb_fp_by_split_cylinder()

        #antigen side physical chemical fingerprint
        fp4 = self.get_physi_chemi_atb_fp()
        
        #interactive force finger print 
        #antigen as the receptor
        fp5 = self.gen_if_complex_fp_atg()
        
        return fp1.append(fp2).append(fp3).append(fp4).append(fp5)

def load_atg_atb(c_id):
    """
    (str) => (antigen, antibody)

    shortcut for loading the antigen and antibody
    
    >>> atg, atb = load_atg_atb("1SLG_D")
    >>> type(atg)
    <class 've.util.structure.mystructure'>
    """
    import os
    from ve.config import data237_complex_root as complex_dir
    
    from ve.util.load_pdb import load_pdb_struct
    
    return load_pdb_struct(os.path.join(complex_dir, c_id, "antigen.pdb"), MyResidue), load_pdb_struct(os.path.join(complex_dir, c_id, "antibody.pdb"), MyResidue)

def test():
    import doctest
    doctest.testmod()

def gen_fp_for_complex(cid):
    """
    (str) => None

    generate fp for complex cid

    >>> gen_fp_for_complex("1SLG_D")
    """
    import os

    from ve.util.load_pdb import load_pdb_struct

    #the fingerprint path
    from ve.config import data237_root
    fp_dir = os.path.join(data237_root, "fp_370")

    #load atg and atb
    atg, atb = load_atg_atb(cid)

    try:
        #init the complex
        c = MyComplex(cid, atg, atb)
    
        #gen the finger print
        fp = c.get_fp()
    except:
        import sys
        sys.stderr.write("%s error occured" %cid)
        
        return

    #the fingerprint path
    fp_path = os.path.join(fp_dir, "%s.fp" %cid)
        
    #output to file
    fp.tofile(fp_path)
        
def main():
    """
    main function
    generate the 370-bit finger prints for all complexes in 237 dataset
    """

    from ve.util.load_pdb import complex_ids
    
    c_ids = complex_ids()
    
    for cid in c_ids:
        print cid
        gen_fp_for_complex(cid)

if __name__ == "__main__":
    import sys
    if len(sys.argv) == 2 and sys.argv[1] == "test":
        test()
    else:
        main()
