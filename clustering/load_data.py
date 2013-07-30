from lmatrix import LMatrix
import numpy as np

def load_mat(path):
    """
    >>> mat = load_mat("data/Multiprot.csv")
    >>> mat.shape
    (166, 166)
    >>> mat["1A2Y_C","3MAC_A"]
    0.41176499999999999
    """
    data = np.loadtxt(open(path,"rb"),delimiter=";",skiprows=1)
    labels = map(lambda s: s.strip(), open(path,"rb").readline().split(";"))
    
    return LMatrix(labels, data = data)

if __name__ == '__main__':
    import doctest
    doctest.testmod()