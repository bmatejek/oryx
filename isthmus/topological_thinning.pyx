import os
import time


cimport cython
cimport numpy as np
import ctypes
import numpy as np


cdef extern from 'cpp-topological-thinning.h':
    void CppTopologicalThinning(const char *prefix, const char *lookup_table_directory, long label)


# generate skeletons for this volume
def TopologicalThinning(prefix, label):
    start_time = time.time()
    
    # convert the numpy arrays to c++
    lut_directory = os.path.dirname(__file__)

    # call the topological skeleton algorithm
    CppTopologicalThinning(prefix, lut_directory, label)
    
    print 'Topological thinning time in {:0.2f} seconds'.format(time.time() - start_time)
