import math
import random



import numpy as np
from numba import jit



from oryx.utilities import dataIO
from oryx.utilities.constants import *



@jit(nopython=True)
def FindNearestSurfacePoint(index, surface_point_cloud, grid_size, resolution):
    # get variables in more convenient form
    zres, yres, xres = grid_size

    iz = index / (yres * xres)
    iy = (index - iz * yres * xres) / xres
    ix = index % xres

    minimum_distance = 2 ** 32
    for point in surface_point_cloud:
        # get this point in linear coordinates
        iw = point / (yres * xres)
        iv = (point - iw * yres * xres) / xres
        iu = point % xres

        # get the distance between these points
        zdiff = resolution[OR_Z] * (iz - iw)
        ydiff = resolution[OR_Y] * (iy - iv)
        xdiff = resolution[OR_X] * (ix - iu)

        distance = zdiff * zdiff + ydiff * ydiff + xdiff * xdiff
        if distance < minimum_distance:
            minimum_distance = distance

    return math.sqrt(minimum_distance)



def EvaluateRadii(prefix, label):
    # get the resolution, surface voxels, and radii for this prefix label pair
    resolution = dataIO.Resolution(prefix)
    grid_size = dataIO.GridSize(prefix)
    
    # surface information
    surface_point_cloud = np.array(dataIO.ReadPoints(prefix, label, 'surfaces'), dtype=np.int64)    
    radii = dataIO.ReadRadii(prefix, label)
    
    mean_absolute_error = 0.0

    count = 0
    for iv, index in enumerate(radii):
        if random.random() < 0.7: continue

        if not (iv % 1000): print '{}/{}'.format(iv + 1, len(radii))
        # get the radius at this index
        radius = radii[index]
        
        # find the nearest point on the surface
        minimum_distance = FindNearestSurfacePoint(index, surface_point_cloud, grid_size, resolution)
        mean_absolute_error += abs(radius - minimum_distance)
        count += 1

    print 'Mean Absolute Error: {:0.2f} nanometers'.format(mean_absolute_error / count)
