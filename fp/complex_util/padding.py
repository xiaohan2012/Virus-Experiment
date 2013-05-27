from collections import Counter
from ve.fp.complex_util.res_spat_dist import ResiduePositionDistribution

import os
from ve.config import data237_cache as cache_dir
dist_cahce_path = os.path.join(cache_dir, "overall_residue_distribution")

from cPickle import dump, load

class OverallSpatialDistribution(object):
    """
    Pad fingerprint based on residue partition within the splitted cylinder
    """
        
    @classmethod
    def overall_dist(cls, complexes, to_cache=True, path = dist_cahce_path):
        """
        (OverallSpatialDistribution, list of Complex) -> ResiduePositionDistribution, ResiduePositionDistribution

        return the Counter recording the maximum residue count in each ring of the splitted cylinder
        """
        atg_dist, atb_dist = ResiduePositionDistribution(), ResiduePositionDistribution()
        
        for c in complexes:
            atg_dist += c.get_atg_res_spat_dist()[0]
            atb_dist += c.get_atb_res_spat_dist()[0]
            
        #to cache
        if to_cache:
            dump((atg_dist, atb_dist), open(path, "w"))
            
        return atg_dist, atb_dist
        
    @classmethod    
    def from_cache(cls,  path = dist_cahce_path):
        return load(open(path, "r"))
    
from ve.fp.fp import BaseComplexFingerprint, BaseResidueFingerprint

from collections import namedtuple
ComplexPaddedFinperPrintPickable = namedtuple("ComplexPaddedFinperPrintPickable", "ring_count mapping")

class PaddedComplexFingerPrint(BaseComplexFingerprint):
    """
    Padded version of complex fingerprint based on fingerprint template
    """
    def __init__(self, ring_count = 0, residue_fp_cls = BaseResidueFingerprint):
        """
        (PaddedComplexFingerPrint, int, int, type) -> NoneType
        """
        self.ring_count = ring_count
        super(PaddedComplexFingerPrint,self).__init__(residue_fp_cls = residue_fp_cls)

    @classmethod
    def make_instance_from_pickable(cls, p, res_fp_cls):
        return cls(p.ring_count, res_fp_cls)
        
    def to_pickable(self):
        """(BaseComplexFingerprint) -> ComplexFinperPrintPickable"""
        return ComplexPaddedFinperPrintPickable(self.ring_count, self.get_mapping())
        
    def fake_fp_str(self, length, value='0',separator=','):
        """
        a string of `length`  separated by `separator` with `value` in between `separator`
        """
        return separator.join([value] * length)
        
    def get_res_fp_length(self):
        """Residue finger print length"""
        return self[self.keys()[0]].bitlength
        
    def fp_str(self, cnt_dist, res_dist, number_type= float):
        """
        (PaddedComplexFingerPrint, ResiduePositionDistribution, ResidueDistribution, type) -> str

        rearange the residue fingerprint and perform some necessary padding
        
        #residue distribution equivalent to dict(int-> list of Residue)
        cnt_dist: the overall residue count distribution
        res_idst: the residue distribution of the complex
        """
        fp_strs = []
        res_fp_len = self.get_res_fp_length()
        #for the ith ring
        for i in xrange(self.ring_count):
            #there is `cnt` residues in this ring
            cnt = cnt_dist[i]
            #for each residue position
            for j in xrange(cnt):
                #haven't reach the end && fps has record for the residue
                if j < len(res_dist[i]) and self.has_key(res_dist[i][j]):
                    fp = self[res_dist[i][j]]
                    fp_strs.append(fp.fp_body(number_type))
                else:
                    fp_strs.append(self.fake_fp_str(res_fp_len))
                    
        return ",".join(fp_strs)
        
    def __eq__(self, other):
        return super(PaddedComplexFingerPrint, self).__eq__(other) and self.ring_count == other.ring_count

    def __str__(self):
        return "Padded Complex Fingerprint with %d rings" %self.ring_count

        
from ve.fp.complex_util.res_spat_dist import ResidueSpatialDistributionTrait
from ve.util.complex import BaseComplex

class ComplexWithResidueSpatialDistribution(BaseComplex, ResidueSpatialDistributionTrait):
    pass
        
def get_overall_residue_distribution():
    from ve.util.load_pdb import load_complexes, complex_ids
    from ve.fp.test.common import GeometryResidue
    
    cs = load_complexes(complex_ids(), complex_cls = ComplexWithResidueSpatialDistribution, residue_cls = GeometryResidue)
    overall_atg_dist, overall_atb_dist = OverallSpatialDistribution.overall_dist(cs)
    
    return overall_atg_dist, overall_atb_dist


if __name__ == '__main__':
    get_overall_residue_distribution()