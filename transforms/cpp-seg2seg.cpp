#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <queue>
#include <unordered_set>
#include <map>
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