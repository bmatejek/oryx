import os
import h5py
import struct


cimport cython
cimport numpy as np
from libcpp cimport bool
import ctypes
import numpy as np



cdef extern from 'cpp-components.h':
     void CppForceConnectivity(char *data, long grid_size[3], long soma_index)


     
def JWRPreprocess(label):
    # read the raw data
    filename = 'raw_data/segmentations/JWR/cell{:03d}_d.h5'.format(label)
    with h5py.File(filename, 'r') as hf:
         data = np.array(hf[hf.keys()[0]]).astype(np.int8)

    # convert the data to a c++ array
    cdef np.ndarray[char, ndim=3, mode='c'] cpp_data = np.ascontiguousarray(data, dtype=ctypes.c_int8)
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_grid_size = np.ascontiguousarray(data.shape, dtype=ctypes.c_int64)
    nentries = data.size
    
    # get the soma location
    soma_filename = 'somae/JWR/{:06d}.pts'.format(label)
    with open(soma_filename, 'rb') as fd:
        zres, yres, xres, npoints, = struct.unpack('qqqq', fd.read(32))
        soma_index, = struct.unpack('q', fd.read(8))

    # call the c++ function
    CppForceConnectivity(&(cpp_data[0,0,0]), &(cpp_grid_size[0]), soma_index)
    
    # save this file
    output_filename = 'raw_data/segmentations/JWR/cell{:03d}_connected_d.h5'.format(label)
    with h5py.File(output_filename, 'w') as hf:
        hf.create_dataset('main', data=cpp_data, compression='gzip')

    
    
