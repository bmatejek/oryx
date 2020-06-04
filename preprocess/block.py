import os
import math
import time
import h5py
import struct



import numpy as np



from oryx.utilities.dataIO import GridSize, ReadPoints
from oryx.utilities.constants import *



def GenerateBlocks(prefix, label):
    # get the number of blocks
    zres, yres, xres = GridSize(prefix)
    block_sizes = (1024, 1024, 1024)

    nzblocks = int(math.ceil(float(zres) / block_sizes[OR_Z]))
    nyblocks = int(math.ceil(float(yres) / block_sizes[OR_Y]))
    nxblocks = int(math.ceil(float(xres) / block_sizes[OR_X]))

    start_time = time.time()

    # get the filename for the segmentation
    point_cloud_filename = 'segmentations/{}/{:06d}.pts'.format(prefix, label)
    if not os.path.exists(point_cloud_filename): return

    points = ReadPoints(prefix, label, 'segmentations')

    points_per_block = {}

    for iz in range(nzblocks):
        for iy in range(nyblocks):
            for ix in range(nxblocks):
                points_per_block[(iz, iy, ix)] = []

    for point in points:
        iz = point // (yres * xres)
        iy = (point - iz * yres * xres) // xres
        ix = point % xres

        zblock = iz // block_sizes[OR_Z]
        yblock = iy // block_sizes[OR_Y]
        xblock = ix // block_sizes[OR_X]

        iz = iz - zblock * block_sizes[OR_Z]
        iy = iy - yblock * block_sizes[OR_Y]
        ix = ix - xblock * block_sizes[OR_X]

        iv = iz * block_sizes[OR_Y] * block_sizes[OR_X] + iy * block_sizes[OR_X] + ix

        points_per_block[(zblock, yblock, xblock)].append(iv)

    for iz in range(nzblocks):
        for iy in range(nyblocks):
            for ix in range(nxblocks):
                output_filename = 'blocks/sections/{:04d}z-{:04d}y-{:04d}x-{:016d}.pts'.format(iz, iy, ix, label)

                with open(output_filename, 'wb') as fd:
                    npoints = len(points_per_block[(iz, iy, ix)])
                    fd.write(struct.pack('qqqq', zres, yres, xres, npoints))
                    fd.write(struct.pack('%sq' % npoints, *points_per_block[(iz, iy, ix)]))

    print ('Completed label {} in {:0.2f} seconds.'.format(label, time.time() - start_time))



def CombineBlocks(prefix):
    # get the number of blocks
    zres, yres, xres = GridSize(prefix)
    block_sizes = (1024, 1024, 1024)

    nzblocks = int(math.ceil(float(zres) / block_sizes[OR_Z]))
    nyblocks = int(math.ceil(float(yres) / block_sizes[OR_Y]))
    nxblocks = int(math.ceil(float(xres) / block_sizes[OR_X]))

    # get the labels for this dataset
    labels = [int(label[:-4]) for label in sorted(os.listdir('segmentations/{}'.format(prefix)))]

    for iz in range(nzblocks):
        for iy in range(nyblocks):
            for ix in range(nxblocks):

                start_time = time.time()
                
                seg = np.zeros((block_sizes[OR_Z], block_sizes[OR_Y], block_sizes[OR_X]), dtype=np.int64)

                for label in labels:
                    input_filename = 'blocks/sections/{:04d}z-{:04d}y-{:04d}x-{:016d}.pts'.format(iz, iy, ix, label)
                    if not os.path.exists(input_filename): continue

                    # read in the points for this label
                    with open(input_filename, 'rb') as fd:
                        input_zres, input_yres, input_xres, input_nelements, = struct.unpack('qqqq', fd.read(32))
                        assert (input_zres == zres)
                        assert (input_yres == yres)
                        assert (input_xres == xres)

                        for _ in range(input_nelements):
                            index, = struct.unpack('q', fd.read(8))

                            iw = index // (block_sizes[OR_Y] * block_sizes[OR_X])
                            iv = (index - iw * block_sizes[OR_Y] * block_sizes[OR_X]) // block_sizes[OR_X]
                            iu = index % block_sizes[OR_X]

                            seg[iw,iv,iu] = label

                output_filename = 'raw_segmentations/Fib25/{:04d}x{:04d}x{:04d}/{:04d}z-{:04d}y-{:04d}x.h5'.format(block_sizes[OR_X], block_sizes[OR_Y], block_sizes[OR_Z], iz, iy, ix)

                with h5py.File(output_filename, 'w') as hf:
                    hf.create_dataset('main', data=seg, compression='gzip')

                print ('Finished block {:04d}z-{:04d}y-{:04d}x in {:0.2f} seconds'.format(iz, iy, ix, time.time() - start_time))
