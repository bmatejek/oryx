import h5py
import numpy as np
import struct


from oryx.utilities import dataIO



def CreateVolumetricSomae(prefix, label, subset):
    if 'JWR':
        somae_filename = 'raw_data/somae/{}/cell{:03d}_d.h5'.format(prefix, label)
        downz, downy, downx = 4, 4, 4

    # get the grid size
    zres, yres, xres = dataIO.GridSize(prefix)
    
    with h5py.File(somae_filename, 'r') as hf:
        somae = np.array(hf['main'])

    # read all segment points
    point_cloud = dataIO.ReadPoints(prefix, label, subset)
    somae_point_cloud = set()

    for index in point_cloud:
        iz = index / (yres * xres)
        iy = (index - iz * yres * xres) / xres
        ix = index % xres

        zdown = iz / downz
        ydown = iy / downy
        xdown = ix / downx

        if somae[zdown, ydown, xdown]:
            somae_point_cloud.add(index)

    somae_point_cloud = list(somae_point_cloud)

    output_filename = 'volumetric_somae/{}/{}/{:06d}.pts'.format(subset, prefix, label)
    with open(output_filename, 'wb') as fd:
        nsegment_points = len(somae_point_cloud)
        fd.write(struct.pack('qqqq', zres, yres, xres, nsegment_points))
        fd.write(struct.pack('%sq' % nsegment_points, *somae_point_cloud))
