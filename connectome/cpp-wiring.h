#ifndef __CPP_WIRING__
#define __CPP_WIRING__


// function calls across cpp files
void CppExtractWiringDiagram(long *segmentation, unsigned char *synapses, long grid_size[3]);


// universal variables and functions

static const int IB_Z = 0;
static const int IB_Y = 1;
static const int IB_X = 2;

#endif