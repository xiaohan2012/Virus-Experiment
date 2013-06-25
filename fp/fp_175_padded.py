from __future__ import division

"""
Padded version of 175-bit finger print(80 + 80 + 15)
"""

from ve.fp.fp_75_gen import get_15bits
from ve.fp.complex import TriangleComplex

from ve.fp.complex_util.padding import PaddedComplexFingerPrint, OverallSpatialDistribution
from ve.fp.complex_util.split_cylinder import SplitCylinderViaComplexPlaneTrait, SplitCylinderViaResiduePlaneTrait
from ve.fp.complex_util.res_spat_dist import ResidueSpatialDistributionTrait
from ve.fp.complex_util.geom import GeometryTrait
from ve.fp.complex_util.cache import ComplexFingerPrintCache as C

overall_atg_dist,overall_atb_dist =  OverallSpatialDistribution.from_cache()

class AxialPlaneBasedComplex(TriangleComplex, #triangle genenration
                             ResidueSpatialDistributionTrait,#padding
                             GeometryTrait):#for complex center calculation
    def gen_fp_str(self, use_tri = True, atg_as_receptor = True, use_cache = True):
        #antigen side
        if use_tri:
            if use_cache:
                C.set_signature("atg_tri_%s_plane" %self.plane_type)
                fps1 = C.load(self.c_id, self,  complex_fp_cls =  PaddedComplexFingerPrint)
            else:
                fps1 = self.gen_fp_by_splitting_cylinder(bases=self.atg.residues,
                                                         targets=[(0,self.get_triangles())],
                                                         fps = PaddedComplexFingerPrint())
        else:
            if use_cache:
                C.set_signature("atg_res_%s_plane" %self.plane_type)
                fps1 = C.load(self.c_id, self,  complex_fp_cls =  PaddedComplexFingerPrint)
            else:
                fps1 = self.gen_fp_by_splitting_cylinder(bases=self.atg.residues,
                                                         targets=[(0,self.atg.residues)], 
                                                         fps = PaddedComplexFingerPrint())
        
        #padded fp str
        _, atg_res_dist = self.get_atg_res_spat_dist()
        str1 = fps1.fp_str(overall_atg_dist, atg_res_dist, number_type=int)
        
        print "str1", len(str1.split(","))
        print overall_atg_dist, atg_res_dist
        
        #antibody side
        if use_cache:
            C.set_signature("atb_res_%s_plane" %self.plane_type)
            fps2 = C.load(self.c_id, self,  complex_fp_cls =  PaddedComplexFingerPrint)
        else:
            fps2 = self.gen_fp_by_splitting_cylinder(bases=self.atb.residues,
                                                     targets=[(0,self.atb.residues)],
                                                     fps = PaddedComplexFingerPrint())
        
        #padded fp str
        _, atb_res_dist = self.get_atb_res_spat_dist()
        str2 = fps2.fp_str(overall_atb_dist, atb_res_dist, number_type=int)

        print "str2", len(str2.split(","))
        print overall_atb_dist, atb_res_dist
        
        #interactive force
        if atg_as_receptor:
            #atg as receptor
            if use_cache:
                C.set_signature("interactive_force_atG_as_rec")
                fps3 = C.load(self.c_id, self, complex_fp_cls =  PaddedComplexFingerPrint)
            else:
                fps3 = get_15bits(receptor = self.atg, binder = self.atb,
                                  fp = PaddedComplexFingerPrint())
            str3 = fps3.fp_str(overall_atg_dist, atg_res_dist, number_type = int)
        else:
            #atb as receptor
            if use_cache:
                C.set_signature("interactive_force_atB_as_rec")
                fps3 = C.load(self.c_id, self, complex_fp_cls =  PaddedComplexFingerPrint)
            else:
                fps3 = get_15bits(binder = self.atg, receptor = self.atb,
                                  fp = PaddedComplexFingerPrint())
            str3 = fps3.fp_str(overall_atb_dist, atb_res_dist, number_type = int)
            
        print "str3", len(str3.split(","))
            
        return ",".join([str1, str2, str3])

class ComplexPlaneBasedComplex(AxialPlaneBasedComplex, SplitCylinderViaComplexPlaneTrait):
    def __init__(self, *args, **kwargs):
        self.plane_type = "complex"
        super(ComplexPlaneBasedComplex,self).__init__(*args, **kwargs)
    
class ResiduePlaneBasedComplex(AxialPlaneBasedComplex, SplitCylinderViaResiduePlaneTrait):
    def __init__(self, *args, **kwargs):
        self.plane_type = "residue"
        super(ResiduePlaneBasedComplex,self).__init__(*args, **kwargs)
    
def main(fp_dir, use_complex_plane = True, atg_as_rec = True, use_tri = True, use_cache = True):
    from ve.fp.fp_80 import Residue
    from ve.util.load_pdb import complex_ids, load_complexes
    from ve.config import data237_fp175_padded_root
    
    complex_cls = ComplexPlaneBasedComplex if use_complex_plane else ResiduePlaneBasedComplex
    cids = ["1FJ1_F"]
    #cids = complex_ids()
    cs = load_complexes(cids, complex_cls = complex_cls, residue_cls = Residue)
    for c in cs:
        fp_str = c.gen_fp_str(atg_as_receptor = atg_as_rec, use_cache = use_cache, use_tri = use_tri)
        print c.c_id
        with open(fp_dir + "/" + c.c_id, "w") as f:
            f.write(fp_str)
            
        try:
            pass
        except:
            print "%s encountered error" %c.c_id

def usage():
    print """
Usage:
    
    python fp_175_padded.py  complex|residue atg|atb tri|res
    """

def bug_list(plane_type, rec_target, tri_or_res):
    tpl = (plane_type, rec_target, tri_or_res)
    if tpl == [("complex", "atg", "res"), ("complex", "atb", "res")]:
        return ["1SLG_D", "2XTJ_A", "1FJ1_E"]
    elif "tri" in tpl:
        """if tri is involved, then the fp lengths are very diversified, so it may needs separate investigation"""
        return None
    elif tpl in [("residue", "atg", "res"), ("residue", "atb", "res")]:
        return ["1SLG_D"]
        
    
if __name__ == '__main__':
    import sys, os
    from ve.config import data237_fp175_padded_root

    try:
        plane_type, rec_target, tri_or_res = sys.argv[1:]
    except ValueError:
        usage()
        sys.exit(0)
        
    fp_dir = os.path.join(data237_fp175_padded_root, "%s-%s-%s" %(plane_type, rec_target, tri_or_res))
    if not os.path.exists(fp_dir):
        os.makedirs(fp_dir)
        
    if plane_type == "complex":
        use_complex_plane = True
    elif plane_type == "residue":
        use_complex_plane = False
    else:
        usage()
        sys.exit(0)

    if rec_target == "atg":
        atg_as_rec = True
    elif rec_target == "atb":
        atg_as_rec = False
    else:
        usage()
        sys.exit(0)

    if tri_or_res == "tri":
        use_tri = True
    elif tri_or_res == "res":
        use_tri = False
    else:
        usage()
        sys.exit(0)        
        
    main(fp_dir, use_complex_plane = use_complex_plane, use_tri = use_tri, atg_as_rec = atg_as_rec, use_cache = True)