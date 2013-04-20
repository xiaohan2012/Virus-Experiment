"""
labeled matrix 
"""
import numpy as np
from numpy import ndarray

np.set_printoptions(precision=3, suppress=True)

class lmatrix(ndarray):
    """Labeled Matrix class"""

    def __new__(cls, labels=[]):
        """labels: list of str"""
        obj = ndarray.__new__(cls, (len(labels), len(labels)))

        obj.labels = labels
        
        #label to index mapping
        obj.label2index_mapping = dict((l,i) for i,l in enumerate(obj.labels))
        
        return obj

    def __array_finalize__(self, obj):    
        if obj is None: return
        
        self.labels = getattr(obj, "labels", None)
        self.label2index_mapping = getattr(obj, "label2index_mapping", None)

    def __getitem__(self, key):
        if isinstance(key, tuple):
            s1,s2 = key
            if isinstance(s1,str) and isinstance(s2,str):#s1 and s2 are both labels
                i1 = self.label2index_mapping[s1]
                i2 = self.label2index_mapping[s2]
                return self.view(ndarray)[i1,i2]
            elif isinstance(s1,str):#s1 is string
                i1 = self.label2index_mapping[s1]
                return self.view(ndarray)[i1,s2]
            elif isinstance(s2,str):#s2 is string
                i2 = self.label2index_mapping[s2]
                return self.view(ndarray)[s1,i2]

                
        return (self.view(ndarray)[key]).view(self.__class__)
    
    def __setitem__(self, key, item):
        if isinstance(key, tuple):
            s1, s2 = key
            if isinstance(s1,str) and isinstance(s2,str):#s1 and s2 are both labels
                i1 = self.label2index_mapping[s1]
                i2 = self.label2index_mapping[s2]
                super(lmatrix, self).__setitem__((i1, i2), item)
                return
            elif isinstance(s1,str):#s1 is label
                i1 = self.label2index_mapping[s1]
                super(lmatrix, self).__setitem__((i1, s2), item)
                return
            elif isinstance(s2,str):#s2 is label
                i2 = self.label2index_mapping[s2]
                super(lmatrix, self).__setitem__((s1, i2), item)
                return
        #otherwise
        super(lmatrix, self).__setitem__(key, item)

#    def __str__(self):
#        return "   %s\n%s" %('  '.join(self.labels), super(lmatrix, self).__str__())

def main():
    labels = ["l1", "l2", "l3"]
    m = lmatrix(labels)

    m[0,:] = [1,2,3]
    m["l2",:] = [4,5,6]
    
    m["l3","l1"] = 7
    m[:,"l3"] = [3,6,9]
    m["l3",:2] = [7,8]

    print m["l1",:]
    print m[:,'l2']
    print m
if __name__ == "__main__":
    main()
