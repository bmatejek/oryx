import os
import h5py
import struct


import numpy as np
from PIL import Image


from oryx.data_structures import meta_data


def GridSize(prefix):
    # return the size of the grid for this prefix
    return meta_data.MetaData(prefix).GridSize()



def ReadImage(filename):
    # return the image corresponding to this file
    im = np.array(Image.open(filename))

    return im



def ReadH5File(filename):
    # return the first h5 dataset from this file
    with h5py.File(filename, 'r') as hf:
        data = np.array(hf[hf.keys()[0]])

    return data



def ReadSegmentationPoints(prefix, label):
    # get the filename for the segmentation
    point_cloud_filename = 'segmentations/{}/{:06d}.pts'.format(prefix, label)

    with open(point_cloud_filename, 'rb') as fd:
        zres, yres, xres, npoints, = struct.unpack('qqqq', fd.read(32))
        point_cloud = struct.unpack('%sq' % npoints, fd.read(8 * npoints))

    return point_cloud



def ReadAllSegmentationPoints(prefix):
    labels = [int(label[:-4]) for label in sorted(os.listdir('segmentations/{}'.format(prefix)))]
    
    point_clouds = {}

    # read all individual point clouds
    for label in labels:
        point_clouds[label] = ReadSegmentationPoints(prefix, label)

    return point_clouds



def ReadSynapsePoints(prefix, label):
    # get the filename for the synapses
    synapse_filename = 'synapses/{}/{:06d}.pts'.format(prefix, label)

    with open(synapse_filename, 'rb') as fd:
        zres, yres, xres, nsynapses, = struct.unpack('qqqq', fd.read(32))
        synapses = struct.unpack('%sq' % nsynapses, fd.read(8 * nsynapses))

    return synapses



def ReadAllSynapsePoints(prefix):
    labels = [int(label[:-4]) for label in sorted(os.listdir('synapses/{}'.format(prefix)))]

    synapses = {}

    # read all synapse points
    for label in labels:
        synapses[label] = ReadSynapsePoints(prefix, label)

    return synapses