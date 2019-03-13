import os
import h5py
import struct


import numpy as np
from PIL import Image


from oryx.data_structures import meta_data


def GridSize(prefix):
    # return the size of the grid for this prefix
    return meta_data.MetaData(prefix).GridSize()



def Resolution(prefix):
    # return the resolution for this prefix
    return meta_data.MetaData(prefix).Resolution()



def ReadImage(filename):
    # return the image corresponding to this file
    im = np.array(Image.open(filename))

    return im



def ReadH5File(filename):
    # return the first h5 dataset from this file
    with h5py.File(filename, 'r') as hf:
        data = np.array(hf[hf.keys()[0]])

    return data



def ReadPoints(prefix, label, dataset):
    # get the filename for the segmentation
    point_cloud_filename = '{}/{}/{:06d}.pts'.format(dataset, prefix, label)

    prefix_zres, prefix_yres, prefix_xres = GridSize(prefix)

    with open(point_cloud_filename, 'rb') as fd:
        zres, yres, xres, npoints, = struct.unpack('qqqq', fd.read(32))
        assert (zres == prefix_zres)
        assert (yres == prefix_yres)
        assert (xres == prefix_xres)
        point_cloud = struct.unpack('%sq' % npoints, fd.read(8 * npoints))

    return point_cloud



def ReadAllPoints(prefix, dataset):
    labels = [int(label[:-4]) for label in sorted(os.listdir('{}/{}'.format(dataset, prefix)))]
    
    point_clouds = {}

    # read all individual point clouds
    for label in labels:
        point_clouds[label] = ReadPoints(prefix, label, dataset)

    return point_clouds



def ReadRadii(prefix, label):
    # get the filename with all of the widths
    radius_filename = 'radii/{}/{:06d}.pts'.format(prefix, label)

    prefix_zres, prefix_yres, prefix_xres = GridSize(prefix)

    radii = {}

    with open(radius_filename, 'rb') as fd:
        zres, yres, xres, nelements, = struct.unpack('qqqq', fd.read(32))
        assert (zres == prefix_zres)
        assert (yres == prefix_yres)
        assert (xres == prefix_xres)

        for _ in range(nelements):
            index, neighbor_index, radius, = struct.unpack('qqd', fd.read(24))
            radii[index] = (neighbor_index, radius)

    # return the dictionary of radii for each skeleton point
    return radii