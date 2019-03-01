import os
import struct



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

            fd.write(struct.pack('qqqq', nsynapses, zres, yres, xres))
            fd.write(struct.pack('%sq' % nsynapses, *synapses))
