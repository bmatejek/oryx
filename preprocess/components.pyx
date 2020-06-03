import os
import time


cimport cython
cimport numpy as np




cdef extern from 'cpp-components.h':
     void CppForceConnectivity(char *prefix, long label)



def Preprocess(prefix, label):
    # return if somae missing
    if not os.path.exists('original_data/segmentations/{}/{:06d}.pts'.format(prefix, label)): return
    start_time = time.time()

    # call the c++ function
    CppForceConnectivity(prefix.encode('utf-8'), label)

    print ('Completed {} in {:0.2f} seconds'.format(label, time.time() - start_time))
