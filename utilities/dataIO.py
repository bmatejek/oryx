import h5py

import numpy as np

from oryx.data_structures import meta_data



def ReadMetaData(prefix):
    # return the meta data for this prefix
    return meta_data.MetaData(prefix)



def Resolution(prefix):
    # return the resolution for this prefix
    return meta_data.MetaData(prefix).Resolution()



def ReadH5File(filename, dataset=None):
    # read the h5py file
    with h5py.File(filename, 'r') as hf:
        # read the first dataset if none given
        if dataset == None: data = np.array(hf[hf.keys()[0]])
        else: data = np.array(hf[dataset])

        # allow affinities and images to not be int64, everything else gets converted
        if data.dtype == np.float32 or data.dtype == np.uint8 or data.dtype == np.int64: return data
        else: return data.astype(np.int64)



def ReadSegmentationData(prefix):
    filename, dataset = meta_data.MetaData(prefix).SegmentationFilename()

    return ReadH5File(filename, dataset).astype(np.int64)



def ReadSynapseData(prefix):
    filename, dataset = meta_data.MetaData(prefix).SynapseFilename()

    return ReadH5File(filename, dataset).astype(np.int64)