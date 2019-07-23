cimport cython
cimport numpy as np

import ctypes
from libcpp cimport bool
import numpy as np



from oryx.utilities import dataIO


cdef extern from 'cpp-seg2seg.h':
    void CppConnectBackground(int *segmentation, long grid_size[3], int max_label)



def ConnectBackground(segmentation, max_label):
    # everything needs to be long ints to work with c++
    assert (segmentation.dtype == np.int32)

    # transform into c array
    cdef np.ndarray[int, ndim=3, mode='c'] cpp_segmentation = np.ascontiguousarray(segmentation, dtype=ctypes.c_int32)
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_grid_size = np.ascontiguousarray(segmentation.shape, dtype=ctypes.c_int64)

    # call the c++ function
    CppConnectBackground(&(cpp_segmentation[0,0,0]), &(cpp_grid_size[0]), max_label)
    
    del cpp_grid_size