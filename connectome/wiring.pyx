import os
import time


cimport cython
cimport numpy as np
import ctypes
import numpy as np



cdef extern from 'cpp-wiring.h':
    void CppExtractWiringDiagram(const char *prefix, long label, const char *lookup_table_directory)

cdef extern from 'cpp-dijkstra.h':
    void CppPostProcess(const char *prefix, long label)



# extract the wiring diagram for this prefix and label
def ExtractWiringDiagram(prefix, label):
    # start running time statistics
    start_time = time.time()

    # the look up table is in the oryx/connectome folder
    lut_directory = os.path.dirname(__file__)

    # call the topological skeleton algorithm
    CppExtractWiringDiagram(prefix, label, lut_directory)

    # print out statistics for wiring extraction
    print 'Extracted wiring diagram in {:0.2f} seconds'.format(time.time() - start_time)



# post process the volume to correct segment errors
def PostProcess(prefix, label):
    # start running time statistics
    start_time = time.time()

    # call the post processing algorithm
    CppPostProcess(prefix, label)

    # print out statistics 
    print 'Ran Dijkstra postprocessing in {:0.2f} seconds'.format(time.time() - start_time)