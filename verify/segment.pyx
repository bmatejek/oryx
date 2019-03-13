import os
import h5py


cimport cython
cimport numpy as np
from libcpp cimport bool
import ctypes
import numpy as np



from oryx.utilities import dataIO



cdef extern from 'cpp-segment.h':
     void CppVerifySNEMISegment(long *data, long nentries, long *point_cloud, long npoints, long label, bool reverse_only)
     void CppVerifyJWRSegment(char *data, long nentries, long *point_cloud, long npoints, bool reverse_only)



# make sure that SNEMI segments match
def SNEMIPointCloud(prefix):
    # read the raw data
    segment_filename = 'raw_data/segmentations/{}/seg.h5'.format(prefix)
    with h5py.File(segment_filename, 'r') as hf:
        seg_data = np.array(hf[hf.keys()[0]]).astype(np.int64)

    labels = [int(label[:-4]) for label in sorted(os.listdir('segmentations/{}'.format(prefix)))]

    # transform raw data to numpy array
    cdef np.ndarray[long, ndim=3, mode='c'] cpp_seg_data = np.ascontiguousarray(seg_data, dtype=ctypes.c_int64)
    nentries = seg_data.size
    # need to define here
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_point_cloud
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_surface_point_cloud

    for label in labels:
        point_cloud = np.array(dataIO.ReadPoints(prefix, label, 'segmentations'))
        npoints = point_cloud.size

        # transform point cloud data to numpy array
        cpp_point_cloud = np.ascontiguousarray(point_cloud, dtype=ctypes.c_int64)

        # call verification function
        CppVerifySNEMISegment(&(cpp_seg_data[0,0,0]), nentries, &(cpp_point_cloud[0]), npoints, label, False)

        surface_point_cloud = np.array(dataIO.ReadPoints(prefix, label, 'surfaces'))
        nsurface_points = surface_point_cloud.size

        # transform point cloud data to numpy array
        cpp_surface_point_cloud = np.ascontiguousarray(surface_point_cloud, dtype=ctypes.c_int64)

        # call verification function
        CppVerifySNEMISegment(&(cpp_seg_data[0,0,0]), nentries, &(cpp_surface_point_cloud[0]), nsurface_points, label, True)

        del cpp_point_cloud
        del cpp_surface_point_cloud



# make sure that SNEMI segments match
def JWRPointCloud(label):

    # need to define here
    cdef np.ndarray[char, ndim=3, mode='c'] cpp_seg_data
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_point_cloud
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_surface_point_cloud
    
    # read the raw data
    segment_filename = 'raw_data/segmentations/JWR/cell{:03d}_connected_d.h5'.format(label)
    with h5py.File(segment_filename, 'r') as hf:
        seg_data = np.array(hf[hf.keys()[0]]).astype(np.int8)

    nentries = seg_data.size

    point_cloud = np.array(dataIO.ReadPoints('JWR', label, 'segmentations'))
    npoints = point_cloud.size

    # transform raw data to numpy array
    cpp_seg_data = np.ascontiguousarray(seg_data, dtype=ctypes.c_int8)

    # transform point cloud data to numpy array
    cpp_point_cloud = np.ascontiguousarray(point_cloud, dtype=ctypes.c_int64)
	
    # call verification function
    CppVerifyJWRSegment(&(cpp_seg_data[0,0,0]), nentries, &(cpp_point_cloud[0]), npoints, False)
    del cpp_point_cloud

    surface_point_cloud = np.array(dataIO.ReadPoints('JWR', label, 'surfaces'))
    nsurface_points = surface_point_cloud.size

    # transform point cloud data to numpy array
    cpp_surface_point_cloud = np.ascontiguousarray(surface_point_cloud, dtype=ctypes.c_int64)

    # call verification function
    CppVerifyJWRSegment(&(cpp_seg_data[0,0,0]), nentries, &(cpp_surface_point_cloud[0]), nsurface_points, True)	
    del cpp_surface_point_cloud 
    
    del cpp_seg_data	
