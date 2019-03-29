import os


cimport cython
cimport numpy as np




cdef extern from 'cpp-components.h':
     void CppForceConnectivity(char *prefix, long label)


     
def Preprocess(prefix, label):
    # return if somae missing
    if not os.path.exists('somae/{}/{:06d}.pts'.format(prefix, label)): return

    # call the c++ function
    CppForceConnectivity(prefix, label)
