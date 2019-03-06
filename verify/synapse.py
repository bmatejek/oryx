import struct
import os
import h5py


from numba import jit
import numpy as np


from oryx.utilities import dataIO
from oryx.utilities.constants import *



########################
### JWR VERIFICATION ###
########################

def JWRSynapse():
    labels = [int(label[:-4]) for label in sorted(os.listdir('segmentations/JWR'))]
    
    prefix = 'JWR'

    # downsample synapse
    downsample = (1, 8, 8)

    # get the grid size
    zres, yres, xres = dataIO.GridSize(prefix)

    # make sure each label has a one to one mapping
    for label in labels:
        synapses = dataIO.ReadSynapsePoints('JWR', label)

        synapses_found = set()

        # read in all of the synapses from the raw data
        with open('raw_data/synapses/JWR/cell{:03d}_d.txt'.format(label), 'r') as fd:
            for synapse in fd:
                attributes = synapse.strip().split()

                ix = int(attributes[0]) / downsample[OR_X]
                iy = int(attributes[1]) / downsample[OR_Y]
                iz = int(attributes[2]) / downsample[OR_Z]

                iv = iz * yres * xres + iy * xres + ix

                assert (iv in synapses)

                synapses_found.add(iv)

        # make sure that every synapse is in raw data
        for synapse in synapses:
            assert (synapse in synapses_found)





############################
### SNEMI3D VERIFICATION ###
############################

@jit(nopython=True)
def SegmentInPointCloud(seg_data, point_clouds):
    # get the grid size
    zres, yres, xres = seg_data.shape

    # make sure every voxel is in the point cloud
    for iz in range(zres):
        print iz
        for iy in range(yres):
            for ix in range(xres):
                if not seg_data[iz,iy,ix]: continue
                
                iv = iz * yres * xres + iy * xres + ix

                assert (iv in point_clouds[seg_data[iz,iy,ix]])



def SNEMISegment(prefix):
    # read in the segment h5 file
    with h5py.File('raw_data/segmentations/{}/seg.h5'.format(prefix), 'r') as hf:
        seg_data = np.array(hf[hf.keys()[0]])

    # read in all of the segmentations
    segmentations = dataIO.ReadAllSegmentationPoints(prefix)

    # get the grid size
    zres, yres, xres = dataIO.GridSize(prefix)

    # make sure that each voxel is in the point cloud
    #SegmentInPointCloud(seg_data, segmentations)

    # make sure every voxel in point clouds is correct
    for segment in segmentations:
        #print segment
        for iv in segmentations[segment]:
            iz = iv / (yres * xres)
            iy = (iv - iz * yres * xres) / xres
            ix = iv % xres

            assert (seg_data[iz,iy,ix] == segment)



def SNEMISynapse(prefix):
    # read in the synapse h5 file
    with h5py.File('raw_data/synapses/{}/synapses.h5'.format(prefix), 'r') as hf:
        syn_data = np.array(hf[hf.keys()[0]])

    # read in all of the synapses 
    syn_per_seg = dataIO.ReadAllSynapsePoints(prefix)

    # get the grid size
    zres, yres, xres = dataIO.GridSize(prefix)

    # make sure that every location actually falls on a synapse
    for segment in syn_per_seg:
        synapses = syn_per_seg[segment]

        for synapse in synapses:
            # make sure this location is a synapse
            iz = synapse / (yres * xres)
            iy = (synapse - iz * yres * xres) / xres
            ix = synapse % xres

            # make sure this is a non-zero location
            assert (syn_data[iz,iy,ix])

    # read in the segmentation h5 file
    with h5py.File('raw_data/segmentations/{}/seg.h5'.format(prefix), 'r') as hf:
        seg_data = np.array(hf[hf.keys()[0]])

    syn_seg_pairs = set()
    syn_seg_found = set()

    # make sure every (synapse, segment) pair occurs once and only once
    for iz in range(zres):
        for iy in range(yres):
            for ix in range(xres):
                synapse = syn_data[iz,iy,ix]
                segment = seg_data[iz,iy,ix]

                # skip background data
                if not synapse or not segment: continue

                syn_seg_pairs.add((synapse, segment))

                # see if this point is in the list of synapses
                iv = iz * yres * xres + iy * xres + ix
                if iv in syn_per_seg[segment]:
                    # make sure there is only one element for this pair
                    assert (not (synapse, segment) in syn_seg_found)

                    # add this pair to the list of found locations
                    syn_seg_found.add((synapse, segment))

    for (synapse, segment) in syn_seg_pairs:
        assert ((synapse, segment) in syn_seg_found)


