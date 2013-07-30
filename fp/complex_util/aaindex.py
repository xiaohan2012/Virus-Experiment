from __future__ import division
from ve.fp.fp import BaseComplexFingerprint
from ve.fp.complex_util.triangle import ResidueTriangleTrait

hydro = {'A': 0.61, 'C': 1.07, 'E': 0.47, 'D': 0.46, 'G': 0.07, 'F': 2.02, 'I': 2.22, 'H': 0.61, 'K': 1.15, 'M': 1.18, 'L': 1.53, 'N': 0.06, 'Q': 0, 'P': 1.95, 'S': 0.05, 'R': 0.6, 'T': 0.05, 'W': 2.65, 'V': 1.32, 'Y': 1.88}
hbond = {'A': 0, 'C': 0, 'E': 1, 'D': 1, 'G': 0, 'F': 0, 'I': 0, 'H': 1, 'K': 2, 'M': 0, 'L': 0, 'N': 2, 'Q': 2, 'P': 0, 'S': 1, 'R': 4, 'T': 1, 'W': 1, 'V': 0, 'Y': 1}
electro = {'A': -0.01, 'C': 0.12, 'E': 0.07, 'D': 0.15, 'G': 0, 'F': 0.03, 'I': -0.01, 'H': 0.08, 'K': 0, 'M': 0.04, 'L': -0.01, 'N': 0.06, 'Q': 0.05, 'P': 0, 'S': 0.11, 'R': 0.04, 'T': 0.04, 'W': 0, 'V': 0.01, 'Y': 0.03}

class AAIndexTrait(ResidueTriangleTrait):
    def get_aaindex_fp(self, iterables, fps):
        fp_len = 3
        
        for r in iterables:
            res_code = r.getCode()
            hy_val = hydro[res_code]
            hb_val = hbond[res_code]
            e_val = electro[res_code]
            
            if not fps.has_res(r):
                fps.add_res(r, fp_len)
            
            fps[r][0] = hy_val
            fps[r][1] = hb_val
            fps[r][2] = e_val
        return fps

    def get_tri_aaindex_fp(self, fps = BaseComplexFingerprint()):
        from collections import Counter
        res_freq = Counter((r  for t in self.get_triangles() for r in t.get_residues()))
        fp_len = 3
        
        for t in self.get_triangles():
            hy_val, hb_val, e_val = 0.0, 0.0, 0.0
            for r in t.get_residues():
                res_code = r.getCode()
                hy_val += hydro[res_code] / res_freq[r]
                hb_val += hbond[res_code] / res_freq[r]
                e_val += electro[res_code] / res_freq[r]
            
            if not fps.has_res(t):
                fps.add_res(t, fp_len)
                
            fps[t][0] = hy_val
            fps[t][1] = hb_val
            fps[t][2] = e_val
        return fps
        
        
    def get_atg_aaindex_fp(self, fps = BaseComplexFingerprint()):
        return self.get_aaindex_fp(self.atg.residues, fps)
        
    def get_atb_aaindex_fp(self, fps = BaseComplexFingerprint()):
        return self.get_aaindex_fp(self.atb.residues, fps)