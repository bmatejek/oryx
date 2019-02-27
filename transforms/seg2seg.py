import h5py
import os
import glob
import re
import h5py



import numpy as np
from PIL import Image



from oryx.utilities import dataIO
from oryx.utilities.constants import *



# quick hack to sort filenames naturally
def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_naturally(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)



def H52ImageStack(prefix, input_directory):
    # create the output directory if it does not exist
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    # get all of the files in the input directory
    filenames = glob.glob('{}/*h5'.format(input_directory))
    sort_naturally(filenames)

    index = 0
    # open each h5 file and write png files
    for filename in filenames:
        with h5py.File(filename, 'r') as hf:
            data = np.array(hf[hf.keys()[0]])
            assert (np.amax(data) < 2**32)
            data = data.astype(np.uint32)

        zres = data.shape[0]
        for iz in range(zres):
            output_filename = 'images/{}/{:06d}.png'.format(prefix, index)
            
            # extract this slice and save the image
            data_slice = data[iz,:,:]
            im = Image.fromarray(data_slice)
            im.save(output_filename)
            
            index += 1



def ImageStack2Sections(prefix):
    zres, yres, xres = dataIO.GridSize(prefix)

    # how many blocks are needed at the given size
    nyblocks = yres / block_size[OR_Y] + (not (yres % block_size[OR_Y] == 0))
    nxblocks = xres / block_size[OR_X] + (not (xres % block_size[OR_X] == 0))

    # begin creating the blocks
    for iz in range(zres):
        # read in the png file
        data_slice = np.array(Image.open('images/{}/{:06d}.png'.format(prefix, iz)))

        # go through each block and create tiles
        for iy in range(nyblocks):
            for ix in range(nxblocks):
                start_x = ix * block_size[OR_X]
                end_x = min((ix + 1) * block_size[OR_X], xres)
                
                start_y = iy * block_size[OR_Y]
                end_y = min((iy + 1) * block_size[OR_Y], yres)

                # write this small tile to file
                output_filename = 'sections/{}/{:06d}-{:06d}x{:06d}-{:06d}x{:06d}.png'.format(prefix, iz, start_y, end_y, start_x, end_x)
                data_segment = data_slice[start_y:end_y,start_x:end_x]
                im = Image.fromarray(data_segment)
                im.save(output_filename)

    

def Sections2Blocks(prefix):
    zres, yres, xres = dataIO.GridSize(prefix)

    # how many blocks are needed at the given size
    nzblocks = zres / block_size[OR_Z] + (not (zres % block_size[OR_Z] == 0))
    nyblocks = yres / block_size[OR_Y] + (not (yres % block_size[OR_Y] == 0))
    nxblocks = xres / block_size[OR_X] + (not (xres % block_size[OR_X] == 0))

    for iyb in range(nyblocks):
        for ixb in range(nxblocks):
            start_x = ixb * block_size[OR_X]
            end_x = min((ixb + 1) * block_size[OR_X], xres)
                
            start_y = iyb * block_size[OR_Y]
            end_y = min((iyb + 1) * block_size[OR_Y], yres)

            for izb in range(nzblocks):
                start_z = izb * block_size[OR_Z]
                end_z = min((izb + 1) * block_size[OR_Z], zres)
                
                # create an empty data array
                data = np.zeros((end_z - start_z, end_y - start_y, end_x - start_x), dtype=np.uint32)

                # read in all of the section images
                index = 0
                for iz in range(start_z, end_z):
                    # get the input filename
                    input_filename = 'sections/{}/{:06d}-{:06d}x{:06d}-{:06d}x{:06d}.png'.format(prefix, iz, start_y, end_y, start_x, end_x)
                    
                    data[index,:,:] = np.array(Image.open(input_filename))
                    index += 1
                
                # save the h5 block
                output_filename = 'blocks/{}/{:06d}x{:06d}-{:06d}x{:06d}-{:06d}x{:06d}.h5'.format(prefix, start_z, end_z, start_y, end_y, start_x, end_x)
                with h5py.File(output_filename, 'w') as hf:
                    hf.create_dataset('main', data=data, compression='gzip')
