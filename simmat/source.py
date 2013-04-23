import os
import numpy as np

def make_fp_str_loader(directory):
    #load_fp_string closure
    def load_fp_string(cid):
        return open(os.path.join(directory, "%s.fp" %cid)).read()

    return load_fp_string

def make_single_line_converter(slice_obj):
    """
    (slice) => (str => np.array)

    make a single line fp string converter specified by the slice object
    
    """
    def converter(string):
        """(str) => np.array
        
        given a fp string, return the np.array
        """
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

from ve.config import data237_root

fp370_atg_dataloader = make_dataloader(os.path.join(data237_root, "fp_370_atg"), make_single_line_converter(slice(0, 370)))
fp370_atb_dataloader = make_dataloader(os.path.join(data237_root, "fp_370_atb"), make_single_line_converter(slice(0, 370)))

first_110_atg_dataloader = make_dataloader(os.path.join(data237_root, "fp_370_atg"), make_single_line_converter(slice(0, 110)))
second_110_atg_dataloader = make_dataloader(os.path.join(data237_root, "fp_370_atg"), make_single_line_converter(slice(110, 220)))
last_150_atg_dataloader = make_dataloader(os.path.join(data237_root, "fp_370_atg"), make_single_line_converter(slice(220, 370)))

first_110_atb_dataloader = make_dataloader(os.path.join(data237_root, "fp_370_atb"), make_single_line_converter(slice(0, 110)))
second_110_atb_dataloader = make_dataloader(os.path.join(data237_root, "fp_370_atb"), make_single_line_converter(slice(110, 220)))
last_150_atb_dataloader = make_dataloader(os.path.join(data237_root, "fp_370_atb"), make_single_line_converter(slice(220, 370)))

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
