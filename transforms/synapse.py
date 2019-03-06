import os
import struct
import h5py


import numpy as np
from numba import jit


from oryx.utilities import dataIO
from oryx.utilities.constants import *



############################
### FIB25 TRANSFORM CODE ###
############################

def Fib25Synapses():
    prefix = 'Fib25'

    # get the grid size to convert to linear coordinates
    zres, yres, xres = dataIO.GridSize(prefix)

    # get the labels for this dataset
    labels = []
    for filename in os.listdir('segmentations/{}'.format(prefix)):
        if not os.path.isfile('segmentations/{}/{}'.format(prefix, filename)): continue        
        labels.append(int(filename[:-4]))

    synapses_per_segment = {}
    for label in labels:
        synapses_per_segment[label] = []

    pre_filename = 'raw_data/synapses/Fib25/synapse_gid_ffn_pre_v1.txt'
    with open(pre_filename, 'r') as fd:
        for line in fd:
            # remove the new line and separate parts
            line = line.strip().split()
            
            # get the id and location
            segment = int(line[0])
            iz = int(line[1])
            iy = int(line[2])
            ix = int(line[3])

            # verify input
            assert (0 <= ix and ix < xres)
            assert (0 <= iy and iy < yres)
            assert (0 <= iz and iz < zres)
            
            # oddly this occurs 
            if not segment in synapses_per_segment: continue

            iv = iz * yres * xres + iy * xres + ix
            synapses_per_segment[segment].append(iv)

    post_filename = 'raw_data/synapses/Fib25/synapse_gid_ffn_post_v1.txt'
    with open(post_filename, 'r') as fd:
        for line in fd:
            # remove the new line and separate parts
            line = line.strip().split()

            # get the id and location
            segment = int(line[1])
            iz = int(line[2])
            iy = int(line[3])
            ix = int(line[4])

            # verify input
            assert (0 <= ix and ix < xres)
            assert (0 <= iy and iy < yres)
            assert (0 <= iz and iz < zres)

            # oddly this occurs
            if not segment in synapses_per_segment: continue

            iv = iz * yres * xres + iy * xres + ix
            synapses_per_segment[segment].append(iv)

    # save all synapses for each label
    for segment in synapses_per_segment:
        filename = 'synapses/Fib25/{:06d}.pts'.format(segment)
        
        with open(filename, 'wb') as fd:
            nsynapses = len(synapses_per_segment[segment])
        
            fd.write(struct.pack('qqqq', zres, yres, xres, nsynapses))
            fd.write(struct.pack('%sq' % nsynapses, *synapses_per_segment[segment]))



##########################
### JWR TRANSFORM CODE ###
##########################

def JWRSynapses():
    prefix = 'JWR'

    # get the grid size to convert to linear coordinates
    zres, yres, xres = dataIO.GridSize(prefix)
    # JWR segmentation is downsampled by factor of 8 in x and y
    downsample_factor = (1, 8, 8)
    

    # get the labels for this dataset
    labels = []
    for filename in os.listdir('segmentations/{}'.format(prefix)):
        if not os.path.isfile('segmentations/{}/{}'.format(prefix, filename)): continue        
        labels.append(int(filename[:-4]))

    for label in labels:
        filename = 'raw_data/synapses/JWR/cell{:03d}_d.txt'.format(label)

        synapses = []

        with open(filename, 'r') as fd:
            for line in fd:
                # remove the new line and separate parts
                line = line.strip().split()

                ix = int(line[0]) / downsample_factor[OR_X]
                iy = int(line[1]) / downsample_factor[OR_Y]
                iz = int(line[2]) / downsample_factor[OR_Z]
                
                # verify input
                assert (0 <= ix and ix < xres)
                assert (0 <= iy and iy < yres)
                assert (0 <= iz and iz < zres)

                iv = iz * yres * xres + iy * xres + ix
                synapses.append(iv)
        
        # save the synapses in the correct form
        filename = 'synapses/{}/{:06d}.pts'.format(prefix, label)

        with open(filename, 'wb') as fd:
            nsynapses = len(synapses)

            fd.write(struct.pack('qqqq', zres, yres, xres, nsynapses))
            fd.write(struct.pack('%sq' % nsynapses, *synapses))



############################
### SNEMI TRANSFORM CODE ###
############################

def FindSynapseSegmentPairs(segmentations, data):
    # get the grid size for the data
    zres, yres, xres = data.shape

    syn_seg_pairs = {}

    for segment in segmentations:
        # go through each point in each segment
        for iv in segmentations[segment]:
            iz = iv / (yres * xres)
            iy = (iv - iz * yres * xres) / xres
            ix = iv % xres
            
            # not a synapse location
            if not data[iz,iy,ix]: continue

            synapse = data[iz,iy,ix]
            if not (synapse, segment) in syn_seg_pairs: syn_seg_pairs[(synapse, segment)] = []
                
            # add this point to the list of overlapping locations
            syn_seg_pairs[(synapse, segment)].append((iz, iy, ix))

    return syn_seg_pairs



@jit(nopython=True)
def MedianCoordinate(coordinates):
    # find the average in the coordinates
    avg_ix = 0
    avg_iy = 0
    avg_iz = 0
    ncoordinates = len(coordinates)

    for (iz, iy, ix) in coordinates:
        avg_ix += float(ix) / ncoordinates
        avg_iy += float(iy) / ncoordinates
        avg_iz += float(iz) / ncoordinates

    # find the coordinate closest to this location
    minimum_distance = float(10e6)
    median_point = (-1, -1, -1)
    for (iz, iy, ix) in coordinates:
        distance = (avg_ix - ix) * (avg_ix - ix) + (avg_iy - iy) * (avg_iy - iy) + (avg_iz - iz) * (avg_iz - iz)
        if distance < minimum_distance:
            minimum_distance = distance
            median_point = (iz, iy, ix)

    return median_point

        

def SNEMISynapses(prefix):
    # read in the segmentation points
    segmentations = dataIO.ReadAllSegmentationPoints(prefix)

    # read in the synapse h5 file
    with h5py.File('raw_data/synapses/{}/synapses.h5'.format(prefix), 'r') as hf:
        data = np.array(hf[hf.keys()[0]])

    # get the grid size for this data
    zres, yres, xres = data.shape

    syn_seg_pairs = FindSynapseSegmentPairs(segmentations, data)
    syn_per_seg = {}

    for (synapse, segment) in syn_seg_pairs:
        coordinates = syn_seg_pairs[(synapse, segment)]

        (iz, iy, ix) = MedianCoordinate(coordinates)

        # convert to linear coordinates
        iv = iz * yres * xres + iy * xres + ix

        if not segment in syn_per_seg:
            syn_per_seg[segment] = []
        
        # add in the median point for this synapse
        syn_per_seg[segment].append(iv)

    # save all of the files
    for segment in syn_per_seg:
        filename = 'synapses/{}/{:06d}.pts'.format(prefix, segment)

        with open(filename, 'wb') as fd:
            nsynapses = len(syn_per_seg[segment])
            fd.write(struct.pack('qqqq', zres, yres, xres, nsynapses))
            fd.write(struct.pack('%sq' % nsynapses, *syn_per_seg[segment]))