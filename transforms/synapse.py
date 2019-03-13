import os
import struct
import h5py
import math
import scipy.spatial


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
    labels = [int(label[:-4]) for label in sorted(os.listdir('segmentations/{}'.format(prefix)))]

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
    
    # JWR has a downsampled segmentation
    downsample_rate = (1, 8, 8)

    # get the labels for this dataset
    labels = [int(label[:-4]) for label in sorted(os.listdir('segmentations/{}'.format(prefix)))]

    for label in labels:
        # get the original filename
        filename = 'raw_data/synapses/JWR/cell{:03d}_d.txt'.format(label)

        # read the segmentation points for this label and convert to numpy array
        surface_point_cloud = dataIO.ReadPoints(prefix, label, 'surfaces')
        segment_point_cloud = set(dataIO.ReadPoints(prefix, label, 'segmentations'))
        npoints = len(surface_point_cloud)

        np_point_cloud = np.zeros((npoints, 3), dtype=np.int32)
        for index, iv in enumerate(surface_point_cloud):
            iz = iv / (yres * xres)
            iy = (iv - iz * yres * xres) / xres
            ix = iv % xres

            np_point_cloud[index,:] = (ix, iy, iz)
            index += 1

        synapses = []

        mse = 0.0
        with open(filename, 'r') as fd:
            for line in fd:
                # remove the new line and separate parts
                line = line.strip().split()

                ix = int(line[0]) / downsample_rate[OR_X]
                iy = int(line[1]) / downsample_rate[OR_Y]
                iz = int(line[2]) / downsample_rate[OR_Z]

                # if already in segment there are no problems
                iv = iz * yres * xres + iy * xres + ix
                if iv in segment_point_cloud:
                    synapses.append(iv)
                    continue
                
                # create a 2D vector for this point
                vec = np.zeros((1, 3), dtype=np.int32)
                vec[0,:] = (ix, iy, iz)

                closest_point = surface_point_cloud[scipy.spatial.distance.cdist(np_point_cloud, vec).argmin()]

                point_iz = closest_point / (yres * xres)
                point_iy = (closest_point - point_iz * yres * xres) / xres
                point_ix = closest_point % xres
                
                distance = math.sqrt((ix - point_ix) * (ix - point_ix) + (iy - point_iy) * (iy - point_iy) + (iz - point_iz) * (iz - point_iz))
                # skip over clearly wrong synapses
                if distance > 100: continue
                mse += distance
                
                synapses.append(closest_point)
                
        print 'Mean Squared Error: {}'.format(mse / len(synapses))

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
    segmentations = dataIO.ReadAllPoints(prefix, 'segmentations')

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
