/* c++ file to preprocess JWR segments */

#include <stdlib.h>
#include <stdio.h>
#include <stack>
#include <unordered_set>
#include "cpp-components.h"
#include <string.h>



static long grid_size[3];
static long nentries;
static long row_size;
static long sheet_size;
static std::unordered_set<long> segment;
static long soma_index;




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


///////////////////////////////////////
//// POINT CLOUD UTILITY FUNCTIONS ////
///////////////////////////////////////

/* conventient I/O function */
void CppPopulatePointCloud(const char *prefix, const char *dataset, long label) {
    // read in the point cloud for this label
    char filename[4096];
    sprintf(filename, "%s/%s/%06ld.pts", dataset, prefix, label);

    FILE *fp = fopen(filename, "rb");
    if (!fp) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }

    long npoints;
    if (fread(&(grid_size[OR_Z]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }
    if (fread(&(grid_size[OR_Y]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }
    if (fread(&(grid_size[OR_X]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }
    if (fread(&npoints, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }
    
    // set global indexing parameters (do here since for loop calls IndicesToIndex)
    nentries = grid_size[OR_Z] * grid_size[OR_Y] * grid_size[OR_X];
    sheet_size = grid_size[OR_Y] * grid_size[OR_X];
    row_size = grid_size[OR_X];

    for (long ip = 0; ip < npoints; ++ip) {
        long voxel_index;
        if (fread(&voxel_index, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }

        long iz = voxel_index / (grid_size[OR_Y] * grid_size[OR_X]);
        long iy = (voxel_index - iz * (grid_size[OR_Y] * grid_size[OR_X])) / grid_size[OR_X];
        long ix = voxel_index % grid_size[OR_X];

        // find the new voxel index
        long iv = IndicesToIndex(ix, iy, iz);

        if (!strcmp(dataset, "original_data/segmentations")) {
          segment.insert(iv);
        }
        else if (!strcmp(dataset, "somae")) {
          soma_index = iv;
        }
        else { fprintf(stderr, "Unrecognized point cloud: %s.\n", dataset); exit(-1); }
    }

    // close file
    fclose(fp);
}



void CppForceConnectivity(char *prefix, long label)
{
  // create new segment set
  segment = std::unordered_set<long>();

  CppPopulatePointCloud(prefix, "original_data/segmentations", label);
  CppPopulatePointCloud(prefix, "somae", label);

  // potential error in soma location 
  if (segment.find(soma_index) == segment.end()) return; 
  
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

    // add the twenty six neighbors to the queue
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
          if (segment.find(neighbor) == segment.end()) continue;

          // add this neighbor
          voxels.push(neighbor);
        }
      }
    }
  }

  char output_filename[4096];
  sprintf(output_filename, "segmentations/%s/%06ld.pts", prefix, label);

  FILE *fp = fopen(output_filename, "wb");
  if (!fp) { fprintf(stderr, "Failed to write to %s\n", output_filename); return ; }

  long npoints = visited.size();
  if (fwrite(&(grid_size[OR_Z]), sizeof(long), 1, fp ) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); return; }
  if (fwrite(&(grid_size[OR_Y]), sizeof(long), 1, fp ) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); return; }
  if (fwrite(&(grid_size[OR_X]), sizeof(long), 1, fp ) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); return; }
  if (fwrite(&npoints, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); return; }

  for (std::unordered_set<long>::iterator it = visited.begin(); it != visited.end(); ++it) {
    long voxel_index = *it;
    if (fwrite(&voxel_index, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); return; }
  }
}
