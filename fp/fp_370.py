from ve.util.residue import BaseResidue

from ve.util.logger import make_logger
logger = make_logger("FP 370")

from ve.fp.residue_util.geom import GeometryTrait
from ve.fp.residue_util.dist import ResDistTrait

class MyResidue(BaseResidue, GeometryTrait, ResDistTrait):
    def __init__(self, res, **kwargs):
        super(MyResidue,self).__init__(res, **kwargs)
        
        
from ve.util.complex import BaseComplex        
from ve.fp.complex_util.split_cylinder import GBCSplitCylinderTrait
from ve.fp.complex_util.physi_chemi import PhysiChemiTrait
from ve.fp.complex_util.interactive_force_fp import InteractiveForceTrait

class MyComplex(BaseComplex, GBCSplitCylinderTrait, PhysiChemiTrait, InteractiveForceTrait):
    def __init__(self, **kwargs):
        super(MyComplex,self).__init__(radius = 20, radius_step = 2, height = 40, height_step = 5, **kwargs)
        
    def get_fp(self, which_as_rec):
        logger.info("gen fp for %s", self.c_id)

        print "antigen side split cylinder fingerprint"
        #antigen side split cylinder fingerprint
        fp1 = self.get_atg_fp_by_split_cylinder()

        print "antigen side physical chemical fingerprint"
        #antigen side physical chemical fingerprint
        fp2 = self.get_physi_chemi_atg_fp()

        print "antibody side split cylinder fingerprint"
        #antibody side split cylinder fingerprint
        fp3 = self.get_atb_fp_by_split_cylinder()

        print "antigen side physical chemical fingerprint"
        #antigen side physical chemical fingerprint
        fp4 = self.get_physi_chemi_atb_fp()

        print "interactive force finger print "
        #interactive force finger print 
        #antigen as the receptor
        if which_as_rec == "atg":
            fp5 = self.gen_if_complex_fp_atg()
        elif which_as_rec == "atb":
            fp5 = self.gen_if_complex_fp_atb()
                        
        
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
        
def main1():
    """
    main function
    generate the 370-bit finger prints for all complexes in 237 dataset
    """

    from ve.util.load_pdb import complex_ids
    
    c_ids = complex_ids()
    
    for cid in c_ids:
        print cid
        gen_fp_for_complex(cid)

def main2(which_as_rec):
    import os
    from ve.util.load_pdb import complex_ids, load_complexes
    from ve.config import data480_root, data480_complex_root, data237_complex_root, data237_root
    ids = ["3DVN_XY"]
    #complex_ids(data480_complex_root)
    cs = load_complexes(ids, directory = data480_complex_root,complex_cls = MyComplex, residue_cls = MyResidue)
    
    preexisted_fp_dir = os.path.join(data237_root, "fp_370_%s" %which_as_rec)
    fp_dir = os.path.join(data480_root, "fp_370_%s" %which_as_rec)

    for c in cs:
        print c.c_id
        fp_path = os.path.join(fp_dir, "%s.fp" %c.c_id)
        if os.path.exists(os.path.join(preexisted_fp_dir, "%s.fp" %c.c_id)) or os.path.exists(fp_path):
            print "%s preexists" %c.c_id
        else:

            try:
                fp = c.get_fp(which_as_rec)
                fp.tofile(fp_path)
            except Exception as e:
                from ve.util.error import get_error_info
                print c.c_id, get_error_info(e)


def usage():
    return """
python fp_370.py test | atg | atb
    """
    
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 2 and sys.argv[1] == "test":
        test()
    elif len(sys.argv) == 2 and sys.argv[1] == "atg":
        """
        3B2U_E
        1A2Y_C
        3DVN_XY
        """
        main2("atg")
    elif len(sys.argv) == 2 and sys.argv[1] == "atb":
        main2("atb")
    else:
        print usage()
