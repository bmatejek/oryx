/* c++ file to preprocess JWR segments */

#include <stdlib.h>
#include <stdio.h>
#include <stack>
#include <unordered_set>
#include "cpp-components.h"



static long nentries;
static long row_size;
static long sheet_size;



static void IndexToIndices(long iv, long &ix, long &iy, long &iz)
{
  iz = iv / sheet_size;
  iy = (iv - iz * sheet_size) / row_size;
  ix = iv % row_size;
}



static long IndicesToIndex(long ix, long iy, long iz)
{
  return iz * sheet_size + iy * row_size + ix;
}



void CppForceConnectivity(char *segmentation, long grid_size[3], long soma_index)
{
  // create the new components array
  nentries = grid_size[OR_Z] * grid_size[OR_Y] * grid_size[OR_X];
  sheet_size = grid_size[OR_Y] * grid_size[OR_X];
  row_size = grid_size[OR_X];

  if (!segmentation[soma_index]) { fprintf(stderr, "Error in the soma index\n"); exit(-1); }
  
  long nonzero = 0;
  for (long iv = 0; iv < nentries; ++iv)
    if (segmentation[iv]) nonzero++;

  std::stack<long> voxels = std::stack<long>();
  voxels.push(soma_index);

  // create a set of vertices that are connected
  std::unordered_set<long> visited = std::unordered_set<long>();

  // perform depth first search
  while (voxels.size()) {
    // remove the pixel from the queue
    long voxel = voxels.top();
    voxels.pop();
    
    // if already visited skip
    if (visited.find(voxel) != visited.end()) continue;

    // label this voxel as visited
    visited.insert(voxel);

    // add the six neighbors to the queue
    long ix, iy, iz;
    IndexToIndices(voxel, ix, iy, iz);

    for (long iw = iz - 1; iw <= iz + 1; ++iw) {
      if (iw < 0 or iw >= grid_size[OR_Z]) continue;
      for (long iv = iy - 1; iv <= iy + 1; ++iv) {
        if (iv < 0 or iv >= grid_size[OR_Y]) continue;
        for (long iu = ix - 1; iu <= ix + 1; ++iu) {
          if (iu < 0 or iu >= grid_size[OR_X]) continue;
          long neighbor = IndicesToIndex(iu, iv, iw);
          if (neighbor == voxel) continue;

          // skip background voxels
          if (!segmentation[neighbor]) continue;

          // add this neighbor
          voxels.push(neighbor);
        }
      }
    }
  }

  // update the segmentation data
  for (long iv = 0; iv < nentries; ++iv)
    segmentation[iv] = 0;

  for (std::unordered_set<long>::iterator it = visited.begin(); it != visited.end(); ++it) 
    segmentation[*it] = 1;
}
