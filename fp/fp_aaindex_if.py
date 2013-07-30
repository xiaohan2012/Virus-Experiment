from __future__ import division

"""
Padded version of 175-bit finger print(80 + 80 + 15)
"""
from ve.util.complex import BaseComplex
from ve.fp.complex_util.interactive_force_fp import InteractiveForceTrait
from ve.fp.complex_util.res_spat_dist import ResidueSpatialDistributionTrait
from ve.fp.complex_util.aaindex import AAIndexTrait

from ve.fp.complex_util.padding import PaddedComplexFingerPrint, OverallSpatialDistribution
from ve.fp.complex_util.cache import ComplexFingerPrintCache as C

overall_atg_dist,overall_atb_dist, overall_tri_dist =  OverallSpatialDistribution.from_cache()

class Complex(BaseComplex,
              InteractiveForceTrait,
              ResidueSpatialDistributionTrait,
              AAIndexTrait):
    def gen_fp_str(self, tri_or_res, atg_as_receptor = True, use_cache = True):
        _, atg_res_dist = self.get_atg_res_spat_dist()
        
        if tri_or_res is "res":
            fps1 = self.get_atg_aaindex_fp(fps = PaddedComplexFingerPrint())
        
            str1 = fps1.fp_str(overall_atg_dist, atg_res_dist, number_type=float)
        else:#tri
            fps1 = self.get_tri_aaindex_fp(fps = PaddedComplexFingerPrint())
        
            #padded fp str
            _, tri_dist = self.get_tri_spat_dist()
            str1 = fps1.fp_str(overall_tri_dist, tri_dist, number_type=float)
            
        #antibody side
        fps2 = self.get_atb_aaindex_fp(fps = PaddedComplexFingerPrint())
            
        #padded fp str
        _, atb_res_dist = self.get_atb_res_spat_dist()
        str2 = fps2.fp_str(overall_atb_dist, atb_res_dist, number_type=float)
        
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
            
        return ",".join([str1, str2, str3])

def main(fp_dir, res_or_tri= "tri", atg_as_rec = True, use_cache = True):
    from ve.fp.fp_80 import Residue
    from ve.util.load_pdb import complex_ids, load_complexes
    from ve.config import data237_fp175_padded_root, data480_complex_root
    
    cids = complex_ids()
    cids = ["1FJ1_F", "3BN9_B","3B9K_EF"]
    
    cs = load_complexes(cids, complex_cls = Complex, residue_cls = Residue)
    for c in cs:
        print c.c_id
        fp_str = c.gen_fp_str(res_or_tri, atg_as_receptor = atg_as_rec, use_cache = use_cache)

        try:
            #fp_str = c.gen_fp_str(res_or_tri, atg_as_receptor = atg_as_rec, use_cache = use_cache)
            with open(fp_dir + "/" + c.c_id + ".fp", "w") as f:
                f.write(fp_str)
        except:
            print "%s encountered error" %c.c_id

        
if __name__ == '__main__':
    import os
    from ve.config import data480_root
    fp_dir = os.path.join(data480_root, "fp_aaindex_if_padded_tri")
    main(fp_dir)