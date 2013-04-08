from customcollections import OrderedDefaultDict
from collections import defaultdict,OrderedDict
from itertools import izip

class BaseResidueFingerprint(OrderedDefaultDict):
    """Fingerprint class for Residue/Triangle
    """
    
    def __init__(self,res,bitlength):
        """(Residue, int) => BaseResidueFingerprint"""

        OrderedDefaultDict.__init__(self,float)
        self.bitlength = bitlength
        self.min_idx = 0
        self.max_idx = self.min_idx + bitlength
        self.res = res

    def fp_body(self):
        """
        () => str
        return the finger print body
        """
        return ",".join(["%.2f" %self[i] for i in xrange(self.min_idx,self.max_idx)])
    
    def fp_header(self):
        """
        () => str
        return the finger print header
        """
        from res_triangle import ResTriangle
        
        #the header string(residue id or triangle id)
        return  isinstance(self.res, ResTriangle) and self.res.fingerprint or str(self.res.resnum)

    def fp_str(self):
        """() => str
        output the fingerprint string
        """
        body_str = self.fp_body()
        header_str = self.fp_header()

        #return the concatenated finger print
        return "%s,%s" %(header_str, body_str)
    
    def append(self,other):
        """
        (residue) => None
        Concatenate it with another residue fingerprint, assuming both of their index are zero-based
        """
        for idx2,val in other.items():
            self[self.max_idx + idx2] = val
        self.max_idx += other.bitlength
        self.bitlength += other.bitlength

#complex fingerprint
class BaseComplexFingerprint(OrderedDict):
    def __init__(self,bitlength=None, residue_fp_cls = BaseResidueFingerprint):
        OrderedDict.__init__(self)
        self.residue_fp_cls = residue_fp_cls
                    
    def has_res(self,res):
        return self.has_key(res)

    def add_res(self,*args):
        res = args[0]
        self[res] = self.residue_fp_cls(*args)

    def fp_str(self):
        return "\n".join([fp.fp_str() for fp in self.values()])

    def tofile(self,fp):
        f = open(fp,"w")
        f.write(self.fp_str())
        f.close()

    def residues(self):
        return self.keys()

class GCBCFingerprint(BaseResidueFingerprint):
    """Geometric-center-based complex finger print"""
    def __init__(self, bitlength):
        BaseResidueFingerprint.__init__(self, None, bitlength)

    def fp_str(self):
        """() => str
        the fingerprint string
        """
        return self.fp_body()
    
class FP0112(BaseComplexFingerprint):
    def __init__(self,targets):
        if len(targets) == 1:
            bitlength=50
        else:
            step = targets[1][0] - targets[0][0]
            bitlength = step * len(targets)
        print bitlength
        BaseComplexFingerprint.__init__(self,bitlength)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
