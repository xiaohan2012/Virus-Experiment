"""
metrics for similarity calculation
"""
import numpy as np

def make_single_line_fp_similatirty_calculator(metric_func):
    """
    given the metric function, 
    return a single line finger print similarity calculator
    """
    #create the closure
    def calc(fp1, fp2):
        return metric_func(fp1, fp2)
    
    #return the closure
    return calc

def corr_coef(vec1, vec2):
    """
    (np.array, np.array) => double

    return the correlation coeffient of vec1 and vec2

    >>> corr_coef([1,2,3], [3,2,1])
    -1.0
    >>> corr_coef([1,2,3], [1,2,3])
    1.0
    >>> corr_coef([1,2,3], [1,2,3])
    1.0
    """
    import numpy as np
    if np.all(vec1 == 0) and np.all(vec2 == 0):
        return 1.0
    elif np.all(vec1 == 0) or np.all(vec2 == 0):
        return 0.0
    else:
        return np.corrcoef(vec1, vec2)[0,1]

def test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    test()
