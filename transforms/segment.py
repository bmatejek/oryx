import struct
import numpy as np
import h5py
import os
import glob


from numba import jit


from oryx.utilities import dataIO
from oryx.utilities.constants import *



#######################################
### ZEBRAFINCH/FIB25 TRANSFORM CODE ###
#######################################

@jit(nopython=True)
def SectionExtractPointCloud(data, z_start):
    zres, yres, xres = data.shape

    max_label = np.amax(data) + 1

    # start with a default element to be removed
    # needed for numba to get types
    point_clouds = [[-1] for iv in range(max_label)]

    for iz in range(zres):
        for iy in range(yres):
            for ix in range(xres):
                if not data[iz,iy,ix]: continue
                
                # need to add z_start since this section is not necessarilly at the bottom
                iv = (z_start + iz) * yres * xres + iy * xres + ix
                point_clouds[data[iz,iy,ix]].append(iv)

    # pop off the default element
    for iv in range(max_label):
        point_clouds[iv] = point_clouds[iv][1:]

    return point_clouds



def H5Section2PointCloud(prefix, filename, section_index, section_width):
    # get grid size as filler for saved files
    zres, yres, xres = dataIO.GridSize(prefix)

    # open this file and read data
    with h5py.File(filename, 'r') as hf:
        data = np.array(hf['main'])

    # z index should not start at zero most times
    z_start = section_width * section_index
    point_clouds = SectionExtractPointCloud(data, z_start)

    for label, point_cloud in enumerate(point_clouds):
        # skip over missing elements from these slices
        if not len(point_cloud): continue

        # write the point cloud to file
        output_filename = 'original_data/segmentations/{}/sections/section-{:03d}-label-{:06d}.pts'.format(prefix, section_index, label)
        with open(output_filename, 'w') as fd:
            npoints = len(point_cloud)
            fd.write(struct.pack('qqqq', zres, yres, xres, npoints))
            fd.write(struct.pack('%sq' % npoints, *point_cloud))



def CombineSectionPointClouds(prefix):
    # get the grid size for this prefix
    zres, yres, xres = dataIO.GridSize(prefix)

    sub_directory = 'original_data/segmentations/{}/sections'.format(prefix)
    
    # get the maximum label from the filenames
    max_label = max(int(filename.split('-')[-1][:-4]) for filename in os.listdir(sub_directory)) + 1

    # go through every label 
    for label in range(max_label):
        filenames = sorted(glob.glob('{}/*-label-{:06d}.pts'.format(sub_directory, label)))
        if not len(filenames): continue
        
        # write the point cloud to this final file
        output_filename = 'original_data/segmentations/{}/{:06d}.pts'.format(prefix, label)
        with open(output_filename, 'wb') as wfd:
            # write nonsense number of points first and overwrite later
            npoints = 0
            wfd.write(struct.pack('qqqq', zres, yres, xres, npoints))
            for filename in filenames:
                with open(filename, 'rb') as rfd:
                    zres, yres, xres, section_npoints = struct.unpack('qqqq', rfd.read(32))
                    point_cloud = struct.unpack('%sq' % section_npoints, rfd.read(8 * section_npoints))
                    wfd.write(struct.pack('%sq' % section_npoints, *point_cloud))
                    npoints += section_npoints

            # reset now that we know the number of true points
            wfd.seek(0)
            wfd.write(struct.pack('qqqq', zres, yres, xres, npoints))



##########################
### JWR TRANSFORM CODE ###
##########################

@jit(nopython=True)
def ExtractJWRPointCloud(data):
    zres, yres, xres = data.shape

    point_cloud = []

    # all non-zero points belong to the segment
    for iz in range(zres):
        for iy in range(yres):
            for ix in range(xres):
                if not data[iz,iy,ix]: continue
                
                iv = iz * yres * xres + iy * xres + ix
                point_cloud.append(iv)

    return point_cloud



def JWRPointCloud(label):
    filename = 'raw_data/segmentations/JWR/cell{:03d}_d.h5'.format(label)

    # open this binary file
    with h5py.File(filename, 'r') as hf:
        # use np.array to decompress
        data = np.array(hf[hf.keys()[0]])
    
    # verify the resolutions match
    zres, yres, xres = dataIO.GridSize('JWR')
    assert (zres == data.shape[OR_Z] and yres == data.shape[OR_Y] and xres == data.shape[OR_X])

    # get all of the non-zero points in a list
    point_cloud = ExtractJWRPointCloud(data)

    # write the point cloud to file
    output_filename = 'original_data/segmentations/JWR/{:06d}.pts'.format(label)
    with open(output_filename, 'wb') as fd:
        npoints = len(point_cloud)
        fd.write(struct.pack('qqqq', zres, yres, xres, npoints))
        fd.write(struct.pack('%sq' % npoints, *point_cloud))



############################
### SNEMI TRANSFORM CODE ###
############################

@jit(nopython=True)
def ExtractSNEMIPointClouds(data):
    zres, yres, xres = data.shape

    max_label = np.amax(data) + 1

    # start with a default element to be removed
    # needed for numba to get types
    point_clouds = [[-1] for iv in range(max_label)]

    for iz in range(zres):
        for iy in range(yres):
            for ix in range(xres):
                if not data[iz,iy,ix]: continue
                
                # need to add z_start since this section is not necessarilly at the bottom
                iv = iz * yres * xres + iy * xres + ix
                point_clouds[data[iz,iy,ix]].append(iv)

    # pop off the default element
    for iv in range(max_label):
        point_clouds[iv] = point_clouds[iv][1:]

    return point_clouds



def SNEMIPointCloud(prefix):
    filename = 'raw_data/segmentations/{}/seg.h5'.format(prefix)

    # open this h5 file
    with h5py.File(filename, 'r') as hf:
        # use np.array to decompress
        data = np.array(hf[hf.keys()[0]])

    # verify the resolutions match
    zres, yres, xres = dataIO.GridSize(prefix)
    assert (zres == data.shape[OR_Z] and yres == data.shape[OR_Y] and xres == data.shape[OR_X])

    # get all of the non-zero points in a list
    point_clouds = ExtractSNEMIPointClouds(data)

    for label, point_cloud in enumerate(point_clouds):
        # skip over missing elements from these slices
        if not len(point_cloud): continue

        # write the point cloud to file
        output_filename = 'segmentations/{}/{:06d}.pts'.format(prefix, label)
        with open(output_filename, 'wb') as fd:
            npoints = len(point_cloud)
            fd.write(struct.pack('qqqq', zres, yres, xres, npoints))
            fd.write(struct.pack('%sq' % npoints, *point_cloud))
