import math
import random
import time
import scipy.spatial



import numpy as np



from oryx.utilities import dataIO
from oryx.utilities.constants import *



def EvaluateRadii(prefix, label):
    start_time = time.time()
    
    # get the resolution, surface voxels, and radii for this prefix label pair
    resolution = dataIO.Resolution(prefix)
    zres, yres, xres = dataIO.GridSize(prefix)
    
    # surface information
    surface_point_cloud = np.array(dataIO.ReadPoints(prefix, label, 'surfaces'), dtype=np.int64)    
    npoints = len(surface_point_cloud)

    np_point_cloud = np.zeros((npoints, 3), dtype=np.int32)
    for index, iv in enumerate(surface_point_cloud):
        iz = iv // (yres * xres)
        iy = (iv - iz * yres * xres) // xres
        ix = iv % xres

        np_point_cloud[index,:] = (ix, iy, iz)
        index += 1
    
    radii = dataIO.ReadRadii(prefix, label)
    
    # keep track of the error over time
    mean_absolute_error = 0.0
    epsilon = 10e-6

    count = 0
    nradii = len(radii)
    for iv, index in enumerate(radii):
        iz = index // (yres * xres)
        iy = (index - iz * yres * xres) // xres
        ix = index % xres

        # create a 2D vector for this point
        vec = np.zeros((1, 3), dtype=np.int32)
        vec[0,:] = (ix, iy, iz)

        # get the radius at this index
        radius = radii[index]

        closest_point = surface_point_cloud[scipy.spatial.distance.cdist(np_point_cloud, vec).argmin()]
        
        point_iz = closest_point // (yres * xres)
        point_iy = (closest_point - point_iz * yres * xres) // xres
        point_ix = closest_point % xres

        # find the nearest point on the surface
        zdiff = resolution[OR_Z] * (point_iz - iz)
        ydiff = resolution[OR_Y] * (point_iy - iy)
        xdiff = resolution[OR_X] * (point_ix - ix)

        minimum_distance = math.sqrt(zdiff * zdiff + ydiff * ydiff + xdiff * xdiff)

        error = abs(radius - minimum_distance)

        mean_absolute_error += error
        count += 1

    print 'Mean Absolute Error: {:0.2f} nanometers'.format(mean_absolute_error / count)

    print time.time() - start_time


