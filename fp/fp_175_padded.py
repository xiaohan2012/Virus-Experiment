from __future__ import division

"""
Padded version of 175-bit finger print
"""

from ve.fp.fp_75_gen import get_15bits
from ve.fp.complex import TriangleComplex

from ve.fp.complex_util.padding import PaddedComplexFingerPrint, OverallSpatialDistribution
from ve.fp.complex_util.split_cylinder import SplitCylinderViaComplexPlaneTrait, SplitCylinderViaResiduePlaneTrait
from ve.fp.complex_util.res_spat_dist import ResidueSpatialDistributionTrait

overall_atg_dist,overall_atb_dist =  OverallSpatialDistribution.from_cache()

class AxialPlaneBasedComplex(TriangleComplex, #triangle genenration
                             ResidueSpatialDistributionTrait):#padding
    def gen_fp_str(self, atg_as_receptor = True):
        #gen fp
        fps1 = self.gen_fp_by_splitting_cylinder(bases=self.atg.residues,
                                     targets=[(0,self.get_triangles())],
                                     fps = PaddedComplexFingerPrint(self.ring_count))
        #padded fp str
        _, atg_res_dist = self.get_atg_res_spat_dist()
        str1 = fps1.fp_str(overall_atg_dist, atg_res_dist, number_type=int)

        #gen fp
        fps2 = self.gen_fp_by_splitting_cylinder(bases=self.atb.residues,
                                                 targets=[(0,self.atb.residues)],
                                                 fps = PaddedComplexFingerPrint(self.ring_count))

        #padded fp str
        _, atb_res_dist = self.get_atb_res_spat_dist()
        str2 = fps2.fp_str(overall_atb_dist, atb_res_dist, number_type=int)

        #gen fp and padded fp str
        if atg_as_receptor:
            fps3 = get_15bits(receptor = self.atg, binder = self.atb, fp = PaddedComplexFingerPrint(self.ring_count))
            str3 = fps3.fp_str(overall_atg_dist, atg_res_dist, number_type = int)
        else:
            fps3 = get_15bits(binder = self.atg, receptor = self.atb, fp = PaddedComplexFingerPrint(self.ring_count))
            str3 = fps3.fp_str(overall_atb_dist, atg_res_dist, number_type = int)

        return ",".join([str1, str2, str3])

class ComplexPlaneBasedComplex(AxialPlaneBasedComplex, SplitCylinderViaComplexPlaneTrait):
    pass

class ResiduePlaneBasedComplex(AxialPlaneBasedComplex, SplitCylinderViaResiduePlaneTrait):
    pass
    
    
def main():
    from ve.fp.fp_80 import Residue
    from ve.util.load_pdb import complex_ids, load_complexes

    cs = load_complexes(complex_ids(), complex_cls = ComplexPlaneBasedComplex, residue_cls = Residue)
    for c in cs:
        try:
            print "%s:%s" %(c.c_id, c.gen_fp_str())
        except:
            print "%s encountered error" %c.c_id
    
if __name__ == '__main__':
    main()