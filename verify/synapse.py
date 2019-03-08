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

def JWRSynapses(label):
    # make sure each label has a one to one mapping
    synapses = dataIO.ReadSynapsePoints('JWR', label)
    point_cloud = set(dataIO.ReadSegmentationPoints('JWR', label))

    for synapse in synapses:
        assert (synapse in point_cloud)



############################
### SNEMI3D VERIFICATION ###
############################

def SNEMISynapses(prefix):
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


