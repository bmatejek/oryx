import h5py

import numpy as np



def ReadH5File(filename, dataset=None):
    # read the h5py file
    with h5py.File(filename, 'r') as hf:
        # read the first dataset if none given
        if dataset == None: data = np.array(hf[hf.keys()[0]])
        else: data = np.array(hf[dataset])