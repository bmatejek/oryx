#ifndef __CPP_WIRING__
#define __CPP_WIRING__


// function calls across cpp files
void CppExtractWiringDiagram(long *segmentation, long *synapses, long grid_size[3]);


// universal variables and functions

static const int OR_Z = 0;
static const int OR_Y = 1;
static const int OR_X = 2;

#endif