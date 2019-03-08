import os
import struct


from oryx.utilities import dataIO



##########################
### JWR TRANSFORM CODE ###
##########################

def JWRSomae(label):
    prefix = 'JWR'

    # get the grid size to convert to linear coordinates
    zres, yres, xres = dataIO.GridSize(prefix)

    # get the original filename
    filename = 'raw_data/somae/JWR/cell{:03d}_d.txt'.format(label)
    
    with open(filename, 'r') as fd:
        # remove the new line
        line = fd.readline().strip().split(',')
        
        ix = int(line[0])
        iy = int(line[1])
        iz = int(line[2])

        # convert to linear index
        iv = iz * yres * xres + iy * xres + ix

        output_filename = 'somae/JWR/{:06d}.pts'.format(label)
        with open(output_filename, 'wb') as fd:
            npoints = 1
            fd.write(struct.pack('qqqq', zres, yres, xres, npoints))
            fd.write(struct.pack('q', iv))

        
        
