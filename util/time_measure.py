"""
program running time meansurement utility 
"""
from time import time
def measure_time(fun):
    def wrapper(*args , **kwargs):
        start = time()
        result = fun(*args , **kwargs)
        print "executing %s:%f" %(fun.__name__ , time() - start)
        return result
    return wrapper

start_time = None

def tic(name = "something"):
    global start_time
    start_time = time()
    #print "%s started" %name
    return start_time

def toc(name = "something"):
    global start_time
    elapse = time() - start_time
    print "%f elapsed when %s finished" %(elapse,name)
    return elapse

