import os
import numpy as np

def make_fp_str_loader(directory):
    #load_fp_string closure
    def load_fp_string(cid):
        return open(os.path.join(directory, "%s.fp" %cid)).read()

    return load_fp_string

def make_fp_str_saver(directory):
    #load_fp_string closure
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    def save_fp_string(cid, string):
        return open(os.path.join(directory, "%s.csv" %cid), "w").write(string)

    return save_fp_string

def make_single_line_converter(slice_obj = None):
    """
    (slice) => (str => np.array)

    make a single line fp string converter specified by the slice object
    
    """
    def converter(string):
        """(str) => np.array
        
        given a fp string, return the np.array
        """
        if slice_obj is None:
            return np.array(map(float, string.strip().split(",")))
        else:
            return np.array(map(float, string.strip().split(",")))[slice_obj]
            
    return converter

def make_dataloader(directory, converter):
    """
    (str, (str:fp string => np.array)) => (str: cid => np.array)

    given the data source directory and string->array converter function
    return the dataloader function

    """
    #get the str_loader function
    str_loader = make_fp_str_loader(directory)
    
    #define the dataloader closure
    def dataloader(cid):
        string = str_loader(cid)
        return converter(string)

    return dataloader

def load_cids(path):
    """() => set of str
    
    return a set of the complex ids
    
    >>> len(load_cids("fp_370_atg.txt"))
    236
    """
    return sorted(list(set([c for l in open(path).readlines() for c in l.split()[:2]])))
    

def load_simmat(path):
    """(str) => lmatrix

    load the lmatrix from file
    """
    
    from lmatrix import lmatrix

    #init matrixc
    m = lmatrix(load_cids(path))

    #set values
    for l in open(path).readlines():
        c1,c2,value = l.split()
        value = float(value)
        m[c1,c2] = value
        m[c2,c1] = value

    return m

def print_simmat_in_csv(path):
    mat = load_simmat(path)
    #print mat[:10,:10]
    print mat.to_csv_str()

        
def test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
#    test()
    import sys
    path = sys.argv[1]
    print_simmat_in_csv(path)

