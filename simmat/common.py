import os
import numpy as np

def make_fp_str_loader(directory):
    #load_fp_string closure
    def load_fp_string(cid):
        return open(os.path.join(directory, "%s.fp" %cid)).read()

    return load_fp_string

def single_line_headless_converter(string):
    return np.array(map(float, string.strip().split(",")))

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

from ve.config import data237_root

fp370_atg_dataloader = make_dataloader(os.path.join(data237_root, "fp_370_atg"), single_line_headless_converter)
fp370_atb_dataloader = make_dataloader(os.path.join(data237_root, "fp_370_atb"), single_line_headless_converter)




def load_simmat(path):
    """(str) => lmatrix

    load the lmatrix from file
    """
    def load_cids():
        """() => set of str
        
        return a set of the complex ids
        
        >>> len(load_cids("fp_370_atg.txt"))
        236
        """
        return set([c for l in open(path).readlines() for c in l.split()[:2]])
    
    from lmatrix import lmatrix

    #init matrixc
    m = lmatrix(load_cids())

    #set values
    for l in open(path).readlines():
        c1,c2,value = l.split()
        value = float(value)
        m[c1,c2] = value
        m[c2,c1] = value

    return m

def main():
    mat = load_simmat("fp_370_atg.txt")
    print mat[:10,:10]

def test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
#    test()
    main()
