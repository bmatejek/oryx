import h5py


from oryx.data_structures import meta_data


def GridSize(prefix):
    # return the size of the grid for this prefix
    return meta_data.MetaData(prefix).GridSize()
