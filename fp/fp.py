from customcollections import OrderedDefaultDict
from collections import defaultdict,OrderedDict
from itertools import izip

class BaseResidueFingerprint(OrderedDefaultDict):
    def __init__(self,res,bitlength):
        OrderedDefaultDict.__init__(self,float)
        self.bitlength = bitlength
        self.min_idx = 0
        self.max_idx = self.min_idx + bitlength
        self.res = res

    def fp_str(self):
        string = ",".join(["%.2f" %self[i] for i in xrange(self.min_idx,self.max_idx)])
        from res_triangle import ResTriangle
        res_id =  isinstance(self.res, ResTriangle) and self.res.fingerprint or str(self.res.resnum)
        return "%s,%s" %(res_id, string)
    
    def append(self,other):
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

class FP0112(BaseComplexFingerprint):
    def __init__(self,targets):
        if len(targets) == 1:
            bitlength=50
        else:
            step = targets[1][0] - targets[0][0]
            bitlength = step * len(targets)
        print bitlength
        BaseComplexFingerprint.__init__(self,bitlength)
