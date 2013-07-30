"""
labeled matrix 
"""
import numpy as np
from numpy import ndarray

#np.set_printoptions(precision=3, suppress=True)

class LMatrix(ndarray):
    """Labeled Matrix class"""

    def __new__(cls, rlabels, clabels = None, data = None):
        """
        rlabels: list of hashable obj, like str, for the rows,
        
        clabels: list of hashable obj, like str, for the cols.
        if clabels not presented, it is the same as rlabels
        
        matrix: None by default,
        if presented, be `np.array` like object
        
        labels_synonyms: iother names for attr labels
        """

        #if clabels not presented, then copy rlabels to it
        if clabels is None: 
            clabels = rlabels
        
        #if matrix is presented, pass it to the new function
        if data is not None:

            r_cnt,c_cnt = data.shape
            
            #the row count and col count should equal
            if r_cnt != len(rlabels) or c_cnt != len(clabels):
                raise ValueError("label size and matrix dimension not match ( %dx%d required, %dx%d given)" %(len(rlabels),
                                                                                                              len(clabels),
                                                                                                              r_cnt,
                                                                                                              c_cnt))
            obj = np.asarray(data).view(cls)
        else:
            rc, cc = len(rlabels), len(clabels)
            obj = ndarray.__new__(cls, (rc, cc), buffer = np.zeros( (rc,cc) ))

            
        obj.rlabels = rlabels
        obj.clabels = clabels
            
        #label to index mapping
        obj.rlabel2index_mapping = dict(map(lambda (i,l): (l,i), enumerate(obj.rlabels)))
        obj.clabel2index_mapping = dict(map(lambda (i,l): (l,i), enumerate(obj.clabels)))
        
        return obj
    
    def __array_finalize__(self, obj):
        if obj is None: return
        
        self.rlabels = getattr(obj, "rlabels", None)
        self.clabels = getattr(obj, "clabels", None)
        
        self.rlabel2index_mapping = getattr(obj, "rlabel2index_mapping", None)
        self.clabel2index_mapping = getattr(obj, "clabel2index_mapping", None)

    def __getitem__(self, key):
        if isinstance(key, tuple):
            ridx, cidx  = key
            if not self.is_traditional_indexing(ridx) and not self.is_traditional_indexing(cidx):
                ridx = self.rlabel2index_mapping[ridx]
                cidx = self.clabel2index_mapping[cidx]
                return self.view(ndarray)[ridx,cidx]
            elif not self.is_traditional_indexing(ridx):
                ridx = self.rlabel2index_mapping[ridx]
                return self.view(ndarray)[ridx,cidx]
            elif not self.is_traditional_indexing(cidx):
                cidx = self.clabel2index_mapping[cidx]
                return self.view(ndarray)[ridx,cidx]
               
        return (self.view(ndarray)[key]).view(self.__class__)
        
    def is_traditional_indexing(self, idx):
        """whether use the traditional indexing method, integer or slice of integer"""
        return isinstance(idx,slice) or isinstance(idx,int)
        
    def __setitem__(self, key, item):
        
        if isinstance(key, tuple):
            ridx, cidx  = key
            if not self.is_traditional_indexing(ridx) and not self.is_traditional_indexing(cidx):
                ridx = self.rlabel2index_mapping[ridx]
                cidx = self.clabel2index_mapping[cidx]
                super(LMatrix, self).__setitem__((ridx, cidx), item)
                return
            elif not self.is_traditional_indexing(ridx):
                ridx = self.rlabel2index_mapping[ridx]
                super(LMatrix, self).__setitem__((ridx, cidx), item)
                return
            elif not self.is_traditional_indexing(cidx):
                cidx = self.clabel2index_mapping[cidx]
                super(LMatrix, self).__setitem__((ridx, cidx), item)
                return
        #otherwise
        super(LMatrix, self).__setitem__(key, item)
