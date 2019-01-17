import math
import os



cimport cython
cimport numpy as np
import ctypes
import numpy as np



cdef extern from 'cpp-wiring.h':
    void CppExtractWiringDiagram(long *segmentation, unsigned char *synapses, long grid_size[3])



# extract the wiring diagram for this segmentation/synapse pair
def ExtractWiringDiagram(segmentation, synapses):
    # everything needs to be long ints and unsigned chars to work with c++
    if not segmentation.dtype == np.int64: segmentation = segmentation.astype(np.int64)
    if not synapses.dtype == np.uint8: synapses = synapses.astype(np.uint8)

    # get the size of the array
    assert (segmentation.shape == synapses.shape)
    grid_size = segmentation.shape
    
    # convert the numpy arrays to c++
    cdef np.ndarray[long, ndim=3, mode='c'] cpp_segmentation = np.ascontiguousarray(segmentation, dtype=ctypes.c_int64)
    cdef np.ndarray[unsigned char, ndim=3, mode='c'] cpp_synaspses = np.ascontiguousarray(synapses, dtype=ctypes.c_uint8)
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_grid_size = np.ascontiguousarray(segmentation.shape, dtype=ctypes.c_int64)
    #lut_directory = os.path.dirname(__file__)

    # call the topological skeleton algorithm
    CppExtractWiringDiagram(&(cpp_segmentation[0,0,0]), &(cpp_synaspses[0,0,0]), &(cpp_grid_size[0]))
