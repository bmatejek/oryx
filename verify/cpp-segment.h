#ifndef __CPP_SEGMENT__
#define __CPP_SEGMENT__


// function calls across cpp files
void CppVerifySNEMISegment(long *data, long nentries, long *point_cloud, long npoints, long label, bool reverse_only);
void CppVerifyJWRSegment(char *data, long nentries, long *point_cloud, long npoints, bool reverse_only);


// universal variables and functions

static const int OR_Z = 0;
static const int OR_Y = 0;
static const int OR_X = 0;

#endif
