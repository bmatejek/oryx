/* c++ file to verify segmentations */

#include <stdio.h>
#include <unordered_set>
#include "cpp-segment.h"



void CppVerifySNEMISegment(long *data, long nentries, long *point_cloud, long npoints, long label, bool reverse_only)
{
    std::unordered_set<long> point_cloud_set = std::unordered_set<long>(2 * npoints);
    for (long ip = 0; ip < npoints; ++ip)
        point_cloud_set.insert(point_cloud[ip]);

    // make sure every element in the point cloud is correct
    for (long ip = 0; ip < npoints; ++ip) {
        long iv = point_cloud[ip];
        if (data[iv] != label) { fprintf(stderr, "Error in SNEMI3D\n"); exit(-1); }
    }

    // skip forward look up for surfaces
    if (reverse_only) return;

    // make sure every entry in data is in point cloud
    for (long iv = 0; iv < nentries; ++iv) {
        // point not found in set is an error
        if ((data[iv] == label) and (point_cloud_set.find(iv) == point_cloud_set.end())) { fprintf(stderr, "Error in SNEMI3D\n"); exit(-1); }
    }
}



void CppVerifyJWRSegment(char *data, long nentries, long *point_cloud, long npoints, bool reverse_only)
{
    std::unordered_set<long> point_cloud_set = std::unordered_set<long>(2 * npoints);
    for (long ip = 0; ip < npoints; ++ip)
        point_cloud_set.insert(point_cloud[ip]);

    // make sure every element in the point cloud is correct
    for (long ip = 0; ip < npoints; ++ip) {
        long iv = point_cloud[ip];
        if (not data[iv]) { fprintf(stderr, "Error in SNEMI3D\n"); exit(-1); }
    }

    // skip forward look up for surfaces
    if (reverse_only) return;

    // make sure every entry in data is in point cloud
    for (long iv = 0; iv < nentries; ++iv) {
        // point not found in set is an error
        if (data[iv] and (point_cloud_set.find(iv) == point_cloud_set.end())) { fprintf(stderr, "Error in SNEMI3D\n"); exit(-1); }
    }
}
