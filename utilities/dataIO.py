import h5py


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