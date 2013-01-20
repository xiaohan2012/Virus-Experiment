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
        return "%d,%s" %(self.res.resnum, string)
    
    def append(self,other):
        for idx2,val in other.items():
            self[self.max_idx + idx2] = val
        self.max_idx += other.bitlength
        self.bitlength += other.bitlength

#complex fingerprint
class BaseComplexFingerprint(OrderedDefaultDict):
    def __init__(self,bitlength=None):
        if bitlength:
            self.min_idx = 0
            self.max_idx = self.min_idx+bitlength
            self.has_bitlength = True
            OrderedDefaultDict.__init__(self, lambda :defaultdict(int))
        else:            
            OrderedDefaultDict.__init__(self, BaseResidueFingerprint)
            self.has_bitlength = False
                    

    def fp_str(self):
        if self.has_bitlength:
            str_segs = []
            for r,fp in self.items():
                str_seg = ','.join(["%d" %self[r][i] for i in xrange(self.min_idx,self.max_idx)])
                str_seg = "%s,%s" %(r,str_seg)
                str_segs.append(str_seg)
            return '\n'.join(str_segs)
        else:    
            return "\n".join([fp.fp_str() for fp in self.values()])
    """
    def fp_str(self):
        str_segs = []
        for r,fp in self.items():
            str_seg = ','.join(["%d" %self[r][i] for i in xrange(self.min_idx,self.max_idx)])
            str_seg = "%s,%s" %(r,str_seg)
            str_segs.append(str_seg)
        return '\n'.join(str_segs)
    """
    def tofile(self,fp):
        f = open(fp,"w")
        f.write(self.fp_str())
        f.close()

class FP0112(BaseComplexFingerprint):
    def __init__(self,targets):
        step = targets[1][0] - targets[0][0]
        bitlength = step * len(targets)
        print bitlength
        BaseComplexFingerprint.__init__(self,bitlength)
