"""
Physical chemical feature based finger print for complex
"""

from collections import defaultdict
from types import MethodType

hydro_dict = {'A':0.61,'C':1.07,'D':0.46,'E':0.47,'F':2.02,'G':0.07,'H':0.61,'I':2.22,'K':1.15,'L':1.53,'M':1.18,'N':0.06,'P':1.95,'Q':0.0,'R':0.6,'S':0.05,'T':0.05,'V':1.32,'W':2.65,'Y':1.88}

charged_dict={'A':-0.01,'C':0.12,'D':0.15,'E':0.07,'F':0.03,'G':0.0,'H':0.08,'I':-0.01,'K':0.0,'L':-0.01,'M':0.04,'N':0.06,'P':0.0,'Q':0.05,'R':0.04,'S':0.11,'T':0.04,'V':0.01,'W':0.0,'Y':0.03}

h_bond_dict={'A':0,'C':0,'D':1,'E':1,'F':0,'G':0,'H':1,'I':0,'K':2,'L':0,'M':0,'N':2,'P':0,'Q':2,'R':4,'S':1,'T':1,'V':0,'W':1,'Y':1}


def physi_chemi_fp(self,others):
    """
    (list of Residue) => None

    Generate the physical-chemical feature based finger print
    
    iterate through a list of residues, group them based on their relaitve location to the complex center.
    calculate scores for each layer in three aspects
    and generate a length-30 finger print
    """
    bitlength=30

    #group surrounding residues according to distance
    layer_res_list = defaultdict(list)
    for other in others:
        #distance to residue 
        dist = self.get_geom_center().dist2point(other.center)

        #within range
        if dist >= self.dist_min and dist <= self.dist_max:
            dist_ind =  int(dist / self.dist_step)
            layer_res_list[dist_ind].append(other)
    
    from ve.fp.fp import HeadlessFingerprint
    #init the fingerprint
    fp = HeadlessFingerprint(bitlength)
    
    #calculate each layer's score in 3 aspects
    for i in xrange(self.layer_count):# at layer i
        #init the scores for this layer
        h_bond,charged,hydro = 0,0,0
        
        #for each residue in this layer, consider its score
        for res in layer_res_list[i]:
            code = res.getCode()

            hydro += hydro_dict[code]
            charged += charged_dict[code]
            h_bond += h_bond_dict[code]
        
        #three scores for this layer
        fp[i], fp[10+i], fp[20+i] = hydro, charged, h_bond

    return fp

def get_physi_chemi_atg_fp(self):
    """finger print for antigen side"""
    return self.physi_chemi_fp(self.atg.residues)

def get_physi_chemi_atb_fp(self):
    """fingerprint for antibody side"""
    return self.physi_chemi_fp(self.atb.residues)


def init_physi_chemi_fp_trait(self,dist_max=20, dist_step=2):
    
    self.physi_chemi_fp = MethodType(physi_chemi_fp,self)
    
    self.dist_min = 0
    self.dist_max = dist_max
    self.dist_step = dist_step
    self.layer_count = self.dist_max / self.dist_step
    
    from ve.fp.complex_util.geom import init_geom_trait
    init_geom_trait(self)
    
    
    self.get_physi_chemi_atg_fp = MethodType(get_physi_chemi_atg_fp, self)
    self.get_physi_chemi_atb_fp = MethodType(get_physi_chemi_atb_fp, self)
