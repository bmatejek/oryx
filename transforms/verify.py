import struct
import os
import h5py


import numpy as np
from oryx.utilities import dataIO



def VerifySynapseLocations(prefix):
    # read in all segment locations
    filenames = sorted(os.listdir('segmentations/{}'.format(prefix)))
    
    zres, yres, xres = dataIO.GridSize(prefix)

    for filename in filenames:
        print '{} {}'.format(prefix, filename)
        segment_filename = 'segmentations/{}/{}'.format(prefix, filename)
        synapse_filename = 'synapses/{}/{}'.format(prefix, filename)

        with open(segment_filename, 'rb') as fd:
            segment_zres, segment_yres, segment_xres, segment_npoints, = struct.unpack('qqqq', fd.read(32))

            segment_points = struct.unpack('%sq' % segment_npoints, fd.read(8 * segment_npoints))
        
        # create a set for easy search
        segment_points = set(segment_points)

        with open(synapse_filename, 'rb') as fd:
            synapse_zres, synapse_yres, synapse_xres, synapse_npoints, = struct.unpack('qqqq', fd.read(32))
            assert (segment_zres == synapse_zres)
            assert (segment_yres == synapse_yres)
            assert (segment_xres == synapse_xres)
            
            synapse_points = struct.unpack('%sq' % synapse_npoints, fd.read(8 * synapse_npoints))
        
        for synapse in synapse_points:
            if not (synapse in segment_points):
                nincorrect += 1
            else:
                ncorrect += 1
        
        print '{} {} {}'.format(max_x, max_y, max_z)

        print '{} {}'.format(ncorrect, nincorrect)
        print '  Point Cloud Size: {}'.format(segment_npoints)
        print '  Synapses Size: {}'.format(synapse_npoints)



def VerifyFib25Segmentations():

    continue
