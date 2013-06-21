from customcollections import OrderedDefaultDict
from collections import defaultdict,OrderedDict
from itertools import izip

import numpy as np

from collections import namedtuple

ResidueFingerPrintPickable = namedtuple("ResidueFingerPrintPickable", "res_id bitlength mapping")

class BaseResidueFingerprint(OrderedDefaultDict):
    """
    Fingerprint class for Residue/Triangle
    """
    
    def __init__(self,res, bitlength, values = None):
        """(Residue, int, dict or list) => BaseResidueFingerprint"""
        OrderedDefaultDict.__init__(self,float)
        self.bitlength = bitlength
        self.min_idx = 0
        self.max_idx = self.min_idx + bitlength
        self.res = res
        
        if values is not None:
            #if it's dict
            if isinstance(values, dict):
                self.set_val(values)
            #if it's list
            elif isinstance(values, list):
                self.set_val(OrderedDict(enumerate(values)))
            else:
                raise ValueError("invalid values type, either dict or list")

    def fp_body(self, number_type=float):
        """
        (type) => str
        specify the numeric display type(float or int)
        return the finger print body
        """
        mapping = {int: "%d",float: "%.2f"}

        return ",".join([mapping[number_type] %self[i] for i in xrange(self.min_idx,self.max_idx)])
    
    def fp_array(self):
        return np.array([self[i] for i in xrange(self.min_idx,self.max_idx)])

    def fp_header(self):
        """
        () => str
        return the finger print header
        """
        from res_triangle import ResTriangle
        
        #the header string(residue id or triangle id)
        return  isinstance(self.res, ResTriangle) and self.res.fingerprint or str(self.res.resnum)

    def fp_str(self, number_type=float):
        """(type) => str
        output the fingerprint string, number display format specified by `number_type`
        """
        body_str = self.fp_body(number_type)
        header_str = self.fp_header()

        #return the concatenated finger print
        return "%s,%s" %(header_str, body_str)
    
    def set_val(self,val_dict):
        """
        (dict of int->int) => None
        copy the keys and values from val_dict to it
        """
        for k,v in val_dict.items():
            self[k] = v
            
    def __add__(self, other):
        """essentially, vector addition"""
        
        fp_list = (self.fp_array() + other.fp_array()).tolist()
        if isinstance(self, HeadlessFingerprint):
            return HeadlessFingerprint(len(fp_list), values=fp_list)
        else:
            return self.__class__(self.res, len(fp_list), values=fp_list)
        
    def __copy__(self):
        """maunally copying"""
        fp = BaseResidueFingerprint(self.res, self.bitlength, self)
        return fp

    def append(self,other):
        """
        (residue) => None
        Concatenate it with another residue fingerprint, assuming both of their index are zero-based
        """
        from copy import copy
        fp = copy(self)
        
        for idx2,val in other.items():
            fp[fp.max_idx + idx2] = val
        fp.max_idx += other.bitlength
        fp.bitlength += other.bitlength
        
        #be functional!
        return fp

    def tofile(self,fp):
        f = open(fp,"w")
        f.write(self.fp_str())
        f.close()
        
    def to_pickable(self):
        """
        (BaseResidueFingerprint) -> ResidueFingerPrintPickable
        """
        d = OrderedDefaultDict(float)
        for k,v in self.items():
            d[k] = v
            
        return ResidueFingerPrintPickable(self.res.res_id, self.bitlength, d)

    @classmethod
    def from_pickable(cls, p, c):
        """(BaseResidueFingerprint) -> (p: Pickable, c: Complex)"""
        res = c.get_res_from_resid(p.res_id)
        res.c = c#set the host complex
        
        return cls(res, p.bitlength, p.mapping)
        
    def __eq__(self, other):
        return sorted(self.items()) == sorted(other.items())
        
ComplexFinperPrintPickable = namedtuple("ComplexFinperPrintPickable", "mapping")

#complex fingerprint
class BaseComplexFingerprint(OrderedDict):
    
    def __init__(self, residue_fp_cls = BaseResidueFingerprint):
        OrderedDict.__init__(self)
        self.residue_fp_cls = residue_fp_cls

    def has_res(self,res):
        return self.has_key(res)

    def add_res(self,*args):
        res = args[0]
        self[res] = self.residue_fp_cls(*args)

    def fp_str(self, number_type= float):
        return "\n".join([fp.fp_str(number_type) for fp in self.values()])

    def tofile(self,fp):
        f = open(fp,"w")
        f.write(self.fp_str())
        f.close()

    def residues(self):
        return self.keys()
    
    def get_bitlength(self):
        return self.values()[0].bitlength
        
    def get_complex(self):
        return self.values()[0].res.c

    def get_mapping(self):
        return OrderedDict(map(lambda (r,fp): (r.res_id, fp.to_pickable()),
                               self.items()))

    def to_pickable(self):
        """(BaseComplexFingerprint) -> ComplexFinperPrintPickable"""
        return ComplexFinperPrintPickable(self.get_mapping())

    @classmethod
    def make_instance_from_pickable(cls, p, res_fp_cls):
        return cls(res_fp_cls)
        
    @classmethod
    def from_pickable(cls, p, c, res_fp_cls = BaseResidueFingerprint):
        """
        (type , ComplexFinperPrintPickable, Complex, type) -> BaseComplexFingerprint
        """
        #construct the complex finger print
        fp = cls.make_instance_from_pickable(p, res_fp_cls)

        #set the residue finger print
        for res_id, res_p in p.mapping.items():
            res_fp = res_fp_cls.from_pickable(res_p, c)
            fp[res_fp.res] = res_fp

        return fp    

    def __eq__(self, other):
        for (res1,fp1),(res2, fp2) in zip(self.items(),other.items()):
            if res1 != res2 or fp1 != fp2:
                print res1,res2, res1 == res2
                print fp1, fp2, fp1 == fp2
                return False
        return True
        
    def get_res_fp_by_resid(self, res_id):
        """
        (ComplexFingerPrint, str) -> ResidueFingerPrint

        Get residue fingerprint by residue id
        """
        sel = filter(lambda r: r.res_id ==res_id, self.keys())
        if sel:
            return self[sel[0]]

    def __str__(self):
        return "%s: %s with %d res fps" %(self.__class__, self.get_complex().c_id, len(self.values()))

    def __repr__(self):
        return str(self)
class HeadlessFingerprint(BaseResidueFingerprint):
    """Geometric-center-based complex finger print"""
    def __init__(self, bitlength, values = None):
        BaseResidueFingerprint.__init__(self, None, bitlength, values=values)

    def fp_str(self, number_type=float):
        """() => str
        the fingerprint string
        """
        return self.fp_body(number_type)

    def __copy__(self):
        """maunally copying"""
        fp = HeadlessFingerprint(self.bitlength, values = self)
        return fp

    
class FP0112(BaseComplexFingerprint):
    def __init__(self,targets):
        if len(targets) == 1:
            bitlength=50
        else:
            step = targets[1][0] - targets[0][0]
            bitlength = step * len(targets)
        BaseComplexFingerprint.__init__(self,bitlength)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
