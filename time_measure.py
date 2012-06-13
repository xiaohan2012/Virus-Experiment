from time import time
def measure_time(fun):
    def wrapper(*args , **kwargs):
        start = time()
        result = fun(*args , **kwargs)
        print "executing %s:%f" %(fun.__name__ , time() - start)
        return result
    return wrapper
