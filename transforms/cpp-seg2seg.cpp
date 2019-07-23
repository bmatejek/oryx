#include <stdio.h>
#include <stdlib.h>
#include <queue>
#include <unordered_set>
#include <ctime>


#define OR_Z 0
#define OR_Y 1
#define OR_X 2


static long nentries;
static long row_size;
static long sheet_size;



static void IndexToIndicies(long iv, long &ix, long &iy, long &iz)
{
    iz = iv / sheet_size;
    iy = (iv - iz * sheet_size) / row_size;
    ix = iv % row_size;
}



static long IndicesToIndex(long ix, long iy, long iz)
{
    return iz * sheet_size + iy * row_size + ix;
}




void CppConnectBackground(int *segmentation, long grid_size[3], int max_label)
{
    // create the new components array
    nentries = grid_size[OR_Z] * grid_size[OR_Y] * grid_size[OR_X];
    sheet_size = grid_size[OR_Y] * grid_size[OR_X];
    row_size = grid_size[OR_X];

    int *components = new int[nentries];
    for (long iv = 0; iv < nentries; ++iv)
        components[iv] = 0;

    // create the queue of labels
    std::queue<long> pixels = std::queue<long>();

    // add one to max label to avoid incorrect merging
    max_label++;

    long current_index = 0;
    int current_label = max_label;

    while (current_index < nentries) {
        // if the current index is already labeled or already seen, continue
        while (current_index < nentries && (components[current_index] || segmentation[current_index])) current_index++;

        // set this component and add to the queue
        components[current_index] = current_label;
        pixels.push(current_index);

        // iterate over all pixels in the queue
        while (pixels.size()) {
            // remove this pixel from the queue
            long pixel = pixels.front();
            pixels.pop();

            // add the six neighbors to the queue
            long iz, iy, ix;
            IndexToIndicies(pixel, ix, iy, iz);

            for (long iw = iz - 1; iw <= iz + 1; ++iw) {
                if (iw < 0 or iw >= grid_size[OR_Z]) continue;
                for (long iv = iy - 1; iv <= iy + 1; ++iv) {
                    if (iv < 0 or iv >= grid_size[OR_Y]) continue;
                    for (long iu = ix - 1; iu <= ix + 1; ++iu) {
                        if (iu < 0 or iu >= grid_size[OR_X]) continue;
                        long neighbor = IndicesToIndex(iu, iv, iw);
           
                        if (!segmentation[neighbor] && !components[neighbor]) {
                            components[neighbor] = current_label;
                            pixels.push(neighbor);
                        }
                    }        
                }
            }
        }
        current_label++;
    }

    // update the segmentation
    for (long iv = 0; iv < nentries; ++iv) {
        if (components[iv]) segmentation[iv] = components[iv];
    }

    delete[] components;
}




void CppFindBackgroundNeighbors(int *segmentation, long grid_size[3], int max_label) 
{
    // create the new components array
    nentries = grid_size[OR_Z] * grid_size[OR_Y] * grid_size[OR_X];
    sheet_size = grid_size[OR_Y] * grid_size[OR_X];
    row_size = grid_size[OR_X];

    // get the max segment including background
    int max_segment = 0;
    for (long iv = 0; iv < nentries; ++iv) {
        if (segmentation[iv] > max_segment) max_segment = segmentation[iv];
    }
    max_segment++;

    // create set of neighbors
    std::unordered_set<int> *neighbors = new std::unordered_set<int>[max_segment];
    for (int iv = 0; iv  < max_segment; ++iv)
        neighbors[iv] = std::unordered_set<int>();

    for (long ii = 0; ii < nentries; ++ii) {
        if ((ii % (grid_size[OR_Y] * grid_size[OR_X])) == 0) printf("%d\n", ii / (grid_size[OR_Y] * grid_size[OR_X]));
        int label = segmentation[ii];
        
        // skip over non-background elements
        if (label <= max_label) continue;

        // go through all neighbors of this background segment
        long iz, iy, ix;
        IndexToIndicies(ii, ix, iy, iz);

        for (long iw = iz - 1; iw <= iz + 1; ++iw) {
            if (iw < 0 or iw >= grid_size[OR_Z]) continue;
            for (long iv = iy - 1; iv <= iy + 1; ++iv) {
                if (iv < 0 or iv >= grid_size[OR_Y]) continue;
                for (long iu = ix - 1; iu <= ix + 1; ++iu) {
                    if (iu < 0 or iu >= grid_size[OR_X]) continue;
                    int neighbor = segmentation[IndicesToIndex(iu, iv, iw)];

                    // skip if same background label
                    if (neighbor == label) continue;

                    // make sure that the neighbor is not background
                    if (neighbor > max_label) {
                        fprintf(stderr, "Error, neighbor is not true label\n"); 
                        exit(-1);
                    }

                    neighbors[label].insert(neighbor);
                }        
            }
        }
    }

    for (long iv = 0; iv < nentries; ++iv) {
        long segment = segmentation[iv];
        if (segment <= max_label) continue;
        else if (neighbors[segment].size() == 1) {
            segmentation[iv] = *(neighbors[segment].begin());
        }
        else segmentation[iv] = 0;
    }

    delete[] neighbors;
}