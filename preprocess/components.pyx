import os


cimport cython
cimport numpy as np




cdef extern from 'cpp-components.h':
     void CppForceConnectivity(char *prefix, long label)


     
def Preprocess(prefix, label):
    # return if somae missing
    if not os.path.exists('somae/{}/{:06d}.pts'.format(prefix, label)): return
    if not os.path.exists('original_data/segmentations/{}/{:06d}.pts'.format(prefix, label)): return

    # call the c++ function
    CppForceConnectivity(prefix.encode('utf-8'), label)
