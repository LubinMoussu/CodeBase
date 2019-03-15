
from time import time
import functools

def function_timing(f):
    def wrap(*args):
        time1 = time()
        ret = f(*args)
        time2 = time()
        print('{:s} function took {:.3f} ms'.format(f.__name__, (time2-time1)*1000.0))

        return ret
    return wrap

def decorator_passing_arguments(f):
    def wrap(*args):
        print("{0}: {1}".format(f.__name__, str(*args) ))
        ret = f(*args)
        return ret
    return wrap

def synchronized(lock):
    """ Synchronization decorator """
    def wrap(f):
        @functools.wraps(f)
        def newFunction(*args, **kw):
            with lock:
                return f(*args, **kw)
        return newFunction
    return wrap

if __name__ == "__main__":

    print("Start")
