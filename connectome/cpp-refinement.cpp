/* c++ file for running skeleton refinement algorithm */

#include "cpp-MinBinaryHeap.h"
#include "cpp-wiring.h"



struct DijkstraData {
    long iv;
    DijkstraData *prev;
    double distance;
    bool visited;
};



void CppSkeletonRefinement(const char *prefix, long label)
{
    // clear the global variables
    segment = std::unordered_map<long, char>();
    synapses = std::unordered_set<long>();
    static std::unordered_map<long, long> dijkstra_map = std::unordered_map<long, long>();
    
    // populate the point clouds with segment voxels and anchor points
    CppPopulatePointCloud(prefix, "skeletons", label);
    CppPopulatePointCloud(prefix, "synapses", label);
    CppPopulatePointCloud(prefix, "somae", label);
    
    // get the number of elements in the skeleton
    long nelements = segment.size();
    
    DijkstraData *voxel_data = new DijkstraData[nelements];
    if (!voxel_data) exit(-1);
    
    // initialize all data
    long index = 0;
    long source_index = -1;
    for (std::unordered_map<long, char>::iterator it = segment.begin(); it != segment.end(); ++it, ++index) {
        voxel_data[index].iv = it->first;
        voxel_data[index].prev = NULL;
        voxel_data[index].distance = infinity;
        voxel_data[index].visited = false;
        dijkstra_map[it->first] = index;

        // this is the soma
        if (it->second == 4) source_index = index;
    }

    // initialize the priority queue
    DijkstraData tmp;
    MinBinaryHeap<DijkstraData *> voxel_heap(&tmp, (&tmp.distance), segment.size());
    
    // insert the source into the heap
    voxel_data[source_index].distance = 0.0;
    voxel_data[source_index].visited = true;
    voxel_heap.Insert(source_index, &(voxel_data[source_index]));
    
    // visit all vertices
    long voxel_index;
    while (!voxel_heap.IsEmpty()) {
        DijkstraData *current = voxel_heap.DeleteMin();
        voxel_index = current->iv;

        // visit all 26 neighbors of this index
        long ix, iy, iz;
        IndexToIndices(voxel_index, ix, iy, iz);

        for (long iw = iz - 1; iw <= iz + 1; ++iw) {
            for (long iv = iy - 1; iv <= iy + 1; ++iv) {
                for (long iu = ix - 1; iu <= ix + 1; ++iu) {
                    // get the linear index for this voxel
                    long neighbor_index = IndicesToIndex(iu, iv, iw);

                    // skip if background
                    if (!segment[neighbor_index]) continue;

                    // get the corresponding neighbor data
                    long dijkstra_index = dijkstra_map[neighbor_index];
                    DijkstraData *neighbor_data = &(voxel_data[dijkstra_index]);

                    // find the distance between these voxels
                    long deltaz = (iw - iz);
                    long deltay = (iv - iy);
                    long deltax = (iu - ix);

                    // get the distance between (ix, iy, iz) and (iu, iv, iw)
                    double distance = sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);

                    // get the distance to get to this voxel through the current voxel (requires a penalty for visiting this voxel)
                    double distance_through_current = current->distance + distance;
                    double distance_without_current = neighbor_data->distance;

                    if (!neighbor_data->visited) {
                        neighbor_data->prev = current;
                        neighbor_data->distance = distance_through_current;
                        neighbor_data->visited = true;
                        voxel_heap.Insert(dijkstra_index, neighbor_data);
                    }
                    else if (distance_through_current < distance_without_current) {
                        neighbor_data->prev = current;
                        neighbor_data->distance = distance_through_current;
                        voxel_heap.DecreaseKey(dijkstra_index, neighbor_data);
                    }
                }
            }
        }
    }

    std::unordered_set<long> wiring_diagram = std::unordered_set<long>();

    // go through all of the synapses and add all of the skeleton points to the source
    for (std::unordered_set<long>::iterator it = synapses.begin(); it != synapses.end(); ++it) {
        // get the voxel and corresponding entry in the dijkstra data frame
        long voxel_index = *it;
        long dijkstra_index = dijkstra_map[voxel_index];

        DijkstraData *data = &(voxel_data[dijkstra_index]);

        while (data != NULL) {
            // add to the list of skeleton points
            long iv = data->iv;

            // convert to unpadded coordinates
            long ix, iy, iz;
            IndexToIndices(iv, ix, iy, iz);
            // unpad x, y, and z
            ix -= 1; iy -= 1; iz -= 1;
            // reconvert to linear coordinates
            iv = iz * (grid_size[OR_Y] - 2) * (grid_size[OR_X] - 2) + iy * (grid_size[OR_X] - 2) + ix;

            wiring_diagram.insert(iv);

            data = data->prev;
        }
    }
    
    long nskeleton_points = wiring_diagram.size();
    char wiring_filename[4096];
    sprintf(wiring_filename, "connectomes/%s/%06ld.pts", prefix, label);

    FILE *fp = fopen(wiring_filename, "wb"); 
    if (!fp) { fprintf(stderr, "Failed to write to %s.\n", wiring_filename); exit(-1); }
    
    // remove padding for file write
    grid_size[OR_Z] -= 2;
    grid_size[OR_Y] -= 2;
    grid_size[OR_X] -= 2;
    
    if (fwrite(&(grid_size[OR_Z]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to write to %s.\n", wiring_filename); exit(-1); }
    if (fwrite(&(grid_size[OR_Y]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to write to %s.\n", wiring_filename); exit(-1); }
    if (fwrite(&(grid_size[OR_X]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to write to %s.\n", wiring_filename); exit(-1); }
    if (fwrite(&nskeleton_points, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to write to %s.\n", wiring_filename); exit(-1); }

    for (std::unordered_set<long>::iterator it = wiring_diagram.begin(); it != wiring_diagram.end(); ++it) {
        long voxel_index = *it;
        if (fwrite(&voxel_index, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to write to %s.\n", wiring_filename); exit(-1); }
    }
    fclose(fp);
    delete[] voxel_data;
}