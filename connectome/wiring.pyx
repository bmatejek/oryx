import math
import os
import time


cimport cython
cimport numpy as np
import ctypes
import numpy as np



cdef extern from 'cpp-wiring.h':
    void CppExtractWiringDiagram(const char *prefix, const char *lookup_table_directory, long *segmentation, long *synapses, long grid_size[3])



# extract the wiring diagram for this segmentation/synapse pair
def ExtractWiringDiagram(prefix, segmentation, synapses):
    # start running time statistics
    start_time = time.time()

    # everything needs to be long ints and unsigned chars to work with c++
    if not segmentation.dtype == np.int64: segmentation = segmentation.astype(np.int64)
    if not synapses.dtype == np.int64: synapses = synapses.astype(np.int64)

    # get the size of the array
    assert (segmentation.shape == synapses.shape)
    grid_size = segmentation.shape
    
    # convert the numpy arrays to c++
    cdef np.ndarray[long, ndim=3, mode='c'] cpp_segmentation = np.ascontiguousarray(segmentation, dtype=ctypes.c_int64)
    cdef np.ndarray[long, ndim=3, mode='c'] cpp_synaspses = np.ascontiguousarray(synapses, dtype=ctypes.c_int64)
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_grid_size = np.ascontiguousarray(grid_size, dtype=ctypes.c_int64)
    lut_directory = os.path.dirname(__file__)

    # call the topological skeleton algorithm
    CppExtractWiringDiagram(prefix, lut_directory, &(cpp_segmentation[0,0,0]), &(cpp_synaspses[0,0,0]), &(cpp_grid_size[0]))

    # print out statistics for wiring extraction
    print 'Extracted wiring diagram in {:0.2f} seconds'.format(time.time() - start_time)