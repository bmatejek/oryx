import h5py


cimport cython
cimport numpy as np
from libcpp cimport bool
import ctypes
import numpy as np



from oryx.utilities import dataIO



cdef extern from 'cpp-segment.h':
     void CppVerifySNEMI3DSegment(long *data, long nentries, long *point_cloud, long npoints, long label, bool reverse_only)


# make sure that SNEMI segments match
def VerifySNEMI3D(prefix):
    # read the raw data
    segment_filename = 'raw_data/segmentations/{}/seg.h5'.format(prefix)
    with h5py.File(segment_filename, 'r') as hf:
        seg_data = np.array(hf[hf.keys()[0]])

    labels = dataIO.Labels(prefix)

    # transform raw data to numpy array
    cdef np.ndarray[long, ndim=3, mode='c'] cpp_seg_data = np.ascontiguousarray(seg_data, dtype=ctypes.c_int64)
    nentries = seg_data.size
    # need to define here
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_point_cloud
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_surface_point_cloud

    for label in labels:
        point_cloud = np.array(dataIO.ReadSegmentationPoints(prefix, label))
        npoints = point_cloud.size

        # transform point cloud data to numpy array
        cpp_point_cloud = np.ascontiguousarray(point_cloud, dtype=ctypes.c_int64)

        # call verification function
        CppVerifySNEMI3DSegment(&(cpp_seg_data[0,0,0]), nentries, &(cpp_point_cloud[0]), npoints, label, False)

        surface_point_cloud = np.array(dataIO.ReadSurfacePoints(prefix, label))
        nsurface_points = surface_point_cloud.size

        # transform point cloud data to numpy array
        cpp_surface_point_cloud = np.ascontiguousarray(surface_point_cloud, dtype=ctypes.c_int64)

        # call verification function
        CppVerifySNEMI3DSegment(&(cpp_seg_data[0,0,0]), nentries, &(cpp_surface_point_cloud[0]), nsurface_points, label, True)