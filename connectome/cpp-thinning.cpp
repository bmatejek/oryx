/* c++ file to extract wiring diagram */

#include "cpp-wiring.h"



// constant variables

static const int lookup_table_size = 1 << 23;

// DO NOT CHANGE THIS ORDERING
static const int NTHINNING_DIRECTIONS = 6;
static const int UP = 0;
static const int DOWN = 1;
static const int NORTH = 2;
static const int SOUTH = 3;
static const int EAST = 4;
static const int WEST = 5;



// lookup tables

static unsigned char *lut_simple;
static std::unordered_map<long, char> erosions;



// mask variables for bitwise operations

static long long_mask[26];
static unsigned char char_mask[8];
static long n26_offsets[26];
static long n6_offsets[6];
static long reverse_direction[6];

static void set_long_mask(void)
{
    long_mask[ 0] = 0x00000001;
    long_mask[ 1] = 0x00000002;
    long_mask[ 2] = 0x00000004;
    long_mask[ 3] = 0x00000008;
    long_mask[ 4] = 0x00000010;
    long_mask[ 5] = 0x00000020;
    long_mask[ 6] = 0x00000040;
    long_mask[ 7] = 0x00000080;
    long_mask[ 8] = 0x00000100;
    long_mask[ 9] = 0x00000200;
    long_mask[10] = 0x00000400;
    long_mask[11] = 0x00000800;
    long_mask[12] = 0x00001000;
    long_mask[13] = 0x00002000;
    long_mask[14] = 0x00004000;
    long_mask[15] = 0x00008000;
    long_mask[16] = 0x00010000;
    long_mask[17] = 0x00020000;
    long_mask[18] = 0x00040000;
    long_mask[19] = 0x00080000;
    long_mask[20] = 0x00100000;
    long_mask[21] = 0x00200000;
    long_mask[22] = 0x00400000;
    long_mask[23] = 0x00800000;
    long_mask[24] = 0x01000000;
    long_mask[25] = 0x02000000;
}

static void set_char_mask(void)
{
    char_mask[0] = 0x01;
    char_mask[1] = 0x02;
    char_mask[2] = 0x04;
    char_mask[3] = 0x08;
    char_mask[4] = 0x10;
    char_mask[5] = 0x20;
    char_mask[6] = 0x40;
    char_mask[7] = 0x80;
}




static void PopulateOffsets(void)
{
    n26_offsets[0] = -1 * grid_size[OR_Y] * grid_size[OR_X] - grid_size[OR_X] - 1;
    n26_offsets[1] = -1 * grid_size[OR_Y] * grid_size[OR_X] - grid_size[OR_X];
    n26_offsets[2] = -1 * grid_size[OR_Y] * grid_size[OR_X] - grid_size[OR_X] + 1;
    n26_offsets[3] = -1 * grid_size[OR_Y] * grid_size[OR_X] - 1;
    n26_offsets[4] = -1 * grid_size[OR_Y] * grid_size[OR_X];
    n26_offsets[5] = -1 * grid_size[OR_Y] * grid_size[OR_X] + 1;
    n26_offsets[6] = -1 * grid_size[OR_Y] * grid_size[OR_X] + grid_size[OR_X] - 1;
    n26_offsets[7] = -1 * grid_size[OR_Y] * grid_size[OR_X] + grid_size[OR_X];
    n26_offsets[8] = -1 * grid_size[OR_Y] * grid_size[OR_X] + grid_size[OR_X] + 1;

    n26_offsets[9] = -1 * grid_size[OR_X] - 1;
    n26_offsets[10] = -1 * grid_size[OR_X];
    n26_offsets[11] = -1 * grid_size[OR_X] + 1;
    n26_offsets[12] = -1;
    n26_offsets[13] = +1;
    n26_offsets[14] = grid_size[OR_X] - 1;
    n26_offsets[15] = grid_size[OR_X];
    n26_offsets[16] = grid_size[OR_X] + 1;

    n26_offsets[17] = grid_size[OR_Y] * grid_size[OR_X] - grid_size[OR_X] - 1;
    n26_offsets[18] = grid_size[OR_Y] * grid_size[OR_X] - grid_size[OR_X];
    n26_offsets[19] = grid_size[OR_Y] * grid_size[OR_X] - grid_size[OR_X] + 1;
    n26_offsets[20] = grid_size[OR_Y] * grid_size[OR_X] - 1;
    n26_offsets[21] = grid_size[OR_Y] * grid_size[OR_X];
    n26_offsets[22] = grid_size[OR_Y] * grid_size[OR_X] + 1;
    n26_offsets[23] = grid_size[OR_Y] * grid_size[OR_X] + grid_size[OR_X] - 1;
    n26_offsets[24] = grid_size[OR_Y] * grid_size[OR_X] + grid_size[OR_X];
    n26_offsets[25] = grid_size[OR_Y] * grid_size[OR_X] + grid_size[OR_X] + 1;

    // use this order to go UP, DOWN, NORTH, SOUTH, EAST, WEST
    // DO NOT CHANGE THIS ORDERING
    n6_offsets[0] = -1 * grid_size[OR_X];
    n6_offsets[1] = grid_size[OR_X];
    n6_offsets[2] = -1 * grid_size[OR_Y] * grid_size[OR_X];
    n6_offsets[3] = grid_size[OR_Y] * grid_size[OR_X];
    n6_offsets[4] = +1;
    n6_offsets[5] = -1;

    // reverse direction points in the opposite direction
    // so erosions points back to the previous node
    // DO NOT CHANGE THIS ORDERING
    reverse_direction[0] = 1;
    reverse_direction[1] = 0;
    reverse_direction[2] = 3;
    reverse_direction[3] = 2;
    reverse_direction[4] = 5;
    reverse_direction[5] = 4;
}



// very simple double linked list data structure

typedef struct {
    long iv, ix, iy, iz;
    void *next;
    void *prev;
} ListElement;

typedef struct {
    void *first;
    void *last;
} List;

typedef struct {
    long iv, ix, iy, iz;
} Voxel;

typedef struct {
    Voxel v;
    ListElement *ptr;
    void *next;
} Cell;

typedef struct {
    Cell *head;
    Cell *tail;
    int length;
} PointList;

typedef struct {
    ListElement *first;
    ListElement *last;
} DoubleList;

List surface_voxels;



static void NewSurfaceVoxel(long iv, long ix, long iy, long iz)
{
    ListElement *LE = new ListElement();
    LE->iv = iv;
    LE->ix = ix;
    LE->iy = iy;
    LE->iz = iz;

    LE->next = NULL;
    LE->prev = surface_voxels.last;

    if (surface_voxels.last != NULL) ((ListElement *) surface_voxels.last)->next = LE;
    surface_voxels.last = LE;
    if (surface_voxels.first == NULL) surface_voxels.first = LE;
}



static void RemoveSurfaceVoxel(ListElement *LE) 
{
    ListElement *LE2;
    if (surface_voxels.first == LE) surface_voxels.first = LE->next;
    if (surface_voxels.last == LE) surface_voxels.last = LE->prev;

    if (LE->next != NULL) {
        LE2 = (ListElement *)(LE->next);
        LE2->prev = LE->prev;
    }
    if (LE->prev != NULL) {
        LE2 = (ListElement *)(LE->prev);
        LE2->next = LE->next;
    }
    delete LE;
}



static void CreatePointList(PointList *s) 
{
    s->head = NULL;
    s->tail = NULL;
    s->length = 0;
}



static void AddToList(PointList *s, Voxel e, ListElement *ptr) 
{
    Cell *newcell = new Cell();
    newcell->v = e;
    newcell->ptr = ptr;
    newcell->next = NULL;

    if (s->head == NULL) {
        s->head = newcell;
        s->tail = newcell;
        s->length = 1;
    }
    else {
        s->tail->next = newcell;
        s->tail = newcell;
        s->length++;
    }
}



static Voxel GetFromList(PointList *s, ListElement **ptr)
{
    Voxel V;
    Cell *tmp;
    V.iv = -1;
    V.ix = -1;
    V.iy = -1;
    V.iz = -1;
    (*ptr) = NULL;
    if (s->length == 0) return V;
    else {
        V = s->head->v;
        (*ptr) = s->head->ptr;
        tmp = (Cell *) s->head->next;
        delete s->head;
        s->head = tmp;
        s->length--;
        if (s->length == 0) {
            s->head = NULL;
            s->tail = NULL;
        }
        return V;
    }
}



static void DestroyPointList(PointList *s) {
    ListElement *ptr;
    while (s->length) GetFromList(s, &ptr);
}



static void InitializeLookupTables(const char *lookup_table_directory)
{
    char lut_filename[4096];
    FILE *lut_file;

    // read the simple lookup table
    sprintf(lut_filename, "%s/lut_simple.dat", lookup_table_directory);
    lut_simple = new unsigned char[lookup_table_size];
    lut_file = fopen(lut_filename, "rb");
    if (!lut_file) {
        fprintf(stderr, "Failed to read %s\n", lut_filename);
        exit(-1);
    }
    if (fread(lut_simple, 1, lookup_table_size, lut_file) != lookup_table_size) {
        fprintf(stderr, "Failed to read %s\n", lut_filename);
        exit(-1);
    }
    fclose(lut_file);

    // set the mask variables
    set_char_mask();
    set_long_mask();
}



static void CollectSurfaceVoxels(void)
{
    // go through all voxels and check their six neighbors
    for (std::unordered_map<long, char>::iterator it = segment.begin(); it != segment.end(); ++it) {
        // all of these elements are either 1 or 3 and in the segment
        long index = it->first;

        long ix, iy, iz;
        IndexToIndices(index, ix, iy, iz);
        // check the 6 neighbors
        for (long iv = 0; iv < NTHINNING_DIRECTIONS; ++iv) {
            long neighbor_index = index + n6_offsets[iv];

            long ii, ij, ik;
            IndexToIndices(neighbor_index, ii, ij, ik);
            
            // skip the fake boundary elements
            if ((ii == 0) or (ii == grid_size[OR_X] - 1)) continue;
            if ((ij == 0) or (ij == grid_size[OR_Y] - 1)) continue;
            if ((ik == 0) or (ik == grid_size[OR_Z] - 1)) continue;
            
            if (segment.find(neighbor_index) == segment.end()) {
                // this location is a boundary so create a surface voxel and break
                // cannot update it->second if it is synapse so need this test!!
                if (it->second == 1) it->second = 2;

                NewSurfaceVoxel(index, ix, iy, iz);

                // note this location as surface in erosions list
                erosions[index] = -1;
                break;
            }
        }
    }
}



static unsigned int Collect26Neighbors(long ix, long iy, long iz)
{
    unsigned int neighbors = 0;
    long index = IndicesToIndex(ix, iy, iz);
    
    // some of these lookups will create a new entry but the region is 
    // shrinking so memory overhead is minimal
    for (long iv = 0; iv < 26; ++iv) {
        if (segment[index + n26_offsets[iv]]) neighbors |= long_mask[iv];
    }

    return neighbors;
}



static bool Simple26_6(unsigned int neighbors)
{
    return lut_simple[(neighbors >> 3)] & char_mask[neighbors % 8];
}



static void DetectSimpleBorderPoints(PointList *deletable_points, int direction)
{
    ListElement *LE = (ListElement *)surface_voxels.first;
    while (LE != NULL) {
        long iv = LE->iv;
        long ix = LE->ix;
        long iy = LE->iy;
        long iz = LE->iz;

        // not a synapse endpoint (need this here since endpoints are on the list of surfaces)
        // this will only be called on things on the surface already so already in unordered_map
        if (segment[iv] == 2) {
            long value = 0;
            // is the neighbor in the corresponding direction not in the segment
            // some of these keys will not exist but will default to 0 value
            // the search region retracts in from the boundary so limited memory overhead
            // the n6_offsets are in the order UP, DOWN, NORTH, SOUTH, EAST, WEST
            value = segment[iv + n6_offsets[direction]];

            // see if the required point belongs to a different segment
            if (!value) {
                unsigned int neighbors = Collect26Neighbors(ix, iy, iz);

                // deletable point
                if (Simple26_6(neighbors)) {
                    Voxel voxel;
                    voxel.iv = iv;
                    voxel.ix = ix;
                    voxel.iy = iy;
                    voxel.iz = iz;
                    AddToList(deletable_points, voxel, LE);
                }
            }
        }
        LE = (ListElement *) LE->next;
    }
}



static long ThinningIterationStep(void)
{
    long changed = 0;

    // iterate through every direction
    for (int direction = 0; direction < NTHINNING_DIRECTIONS; ++direction) {
        PointList deletable_points;
        ListElement *ptr;

        CreatePointList(&deletable_points);
        DetectSimpleBorderPoints(&deletable_points, direction);

        while (deletable_points.length) {
            Voxel voxel = GetFromList(&deletable_points, &ptr);

            long index = voxel.iv;
            long ix = voxel.ix;
            long iy = voxel.iy;
            long iz = voxel.iz;

            unsigned int neighbors = Collect26Neighbors(ix, iy, iz);
            if (Simple26_6(neighbors)) {
                // delete the simple point
                segment[index] = 0;

                // add the new surface voxels
                for (long ip = 0; ip < NTHINNING_DIRECTIONS; ++ip) {
                    long neighbor_index = index + n6_offsets[ip];

                    // previously not on the surface but is in the object
                    if (segment[neighbor_index] && erosions.find(neighbor_index) == erosions.end()) {

                        long iu, iv, iw;
                        IndexToIndices(neighbor_index, iu, iv, iw);
                        NewSurfaceVoxel(neighbor_index, iu, iv, iw);
                        // do not convert isthmuses to surface points
                        if (segment[neighbor_index] == 1) segment[neighbor_index] = 2;
                        // use ip here so it points to voxel before it 
                        // need to reverse the direction so that we are point back at the previous voxel
                        // DO NOT CHANGE
                        erosions[neighbor_index] = reverse_direction[ip];
                    }
                }

                // remove this from the surface voxels
                RemoveSurfaceVoxel(ptr);
                changed += 1;
            }
        }
        DestroyPointList(&deletable_points);
    }


    // return the number of changes
    return changed;
}



static void SequentialThinning(const char *prefix, long label)
{
    // create a vector of surface voxels
    CollectSurfaceVoxels();
    int iteration = 0;
    long changed = 0;
    do {
        changed = ThinningIterationStep();
        iteration++;
        printf("  Iteration %d deleted %ld points\n", iteration, changed);
    } while (changed);
}




void CppSkeletonGeneration(const char *prefix, long label, const char *lookup_table_directory)
{
    // start timing statistics
    clock_t start_time = clock();

    // create (and clear) the global variables
    segment = std::unordered_map<long, char>(10000000);
    // not needed here but for refinement so need to initialize anyway
    synapses = std::unordered_set<long>(); 
    erosions = std::unordered_map<long, char>(10000000);

    // initialize all of the lookup tables
    InitializeLookupTables(lookup_table_directory);

    // populate the point clouds with segment voxels and anchor points
    CppPopulatePointCloud(prefix, "segmentations", label);
    CppPopulatePointCloud(prefix, "synapses", label);
    CppPopulatePointCloud(prefix,  "somae", label);
    
    // get the number of points
    long initial_points = segment.size();
    printf("Label %ld initial points: %ld\n", label, initial_points);

    // can  use offsets since all paramters are offset by 1
    // needs to happen after PopulatePointCloud()
    PopulateOffsets();

    // call the sequential thinning algorithm
    SequentialThinning(prefix, label);

    // count the number of remaining points
    long num = 0;
    ListElement *LE = (ListElement *) surface_voxels.first;
    while (LE != NULL) {
        num++;
        LE = (ListElement *)LE->next;
    }

    // create an output file for the points
    char output_filename[4096];
    sprintf(output_filename, "skeletons/%s/%06ld.pts", prefix, label);

    FILE *wfp = fopen(output_filename, "wb");
    if (!wfp) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }

    // unpad the grid size
    grid_size[OR_Z] -= 2;
    grid_size[OR_Y] -= 2;
    grid_size[OR_X] -= 2;

    // write the number of elements
    if (fwrite(&(grid_size[OR_Z]), sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
    if (fwrite(&(grid_size[OR_Y]), sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
    if (fwrite(&(grid_size[OR_X]), sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
    if (fwrite(&num, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }

    printf("  Remaining voxels: %ld\n", num);

    while (surface_voxels.first != NULL) {
        // get the surface voxels
        ListElement *LE = (ListElement *) surface_voxels.first;

        // get the coordinates for this skeleton point in the non-cropped segmentation
        long iz = LE->iz - 1;
        long iy = LE->iy - 1;
        long ix = LE->ix - 1;
        long iv = iz * grid_size[OR_X] * grid_size[OR_Y] + iy * grid_size[OR_X] + ix;

        if (fwrite(&iv, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }

        // remove this voxel
        RemoveSurfaceVoxel(LE);
    }

    // close the I/O files
    fclose(wfp);
       
    delete[] lut_simple;

    double total_time = (double) (clock() - start_time) / CLOCKS_PER_SEC;

    char time_filename[4096];
    sprintf(time_filename, "timings/skeletons/%s-%06ld.time", prefix, label);

    FILE *tfp = fopen(time_filename, "wb");
    if (!tfp) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }

    // write the number of points and the total time to file
    if (fwrite(&initial_points, sizeof(long), 1, tfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }
    if (fwrite(&total_time, sizeof(double), 1, tfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }

    // close file
    fclose(tfp);
}



void CppGenerateWidths(const char *prefix, long label, double resolution[3])
{
    // start timing statistics
    clock_t start_time = clock();

    // create (and clear) the global variable segment to read in thinning only
    segment = std::unordered_map<long, char>();

    // populate the point clouds with segment voxels and anchor points
    CppPopulatePointCloud(prefix, "skeletons", label);

    // can  use offsets since all paramters are offset by 1
    // needs to happen after PopulatePointCloud()
    PopulateOffsets();

    // create an output file for the points
    char output_filename[4096];
    sprintf(output_filename, "radii/%s/%06ld.pts", prefix, label);

    FILE *wfp = fopen(output_filename, "wb");
    if (!wfp) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }

    // unpad the grid size
    grid_size[OR_Z] -= 2;
    grid_size[OR_Y] -= 2;
    grid_size[OR_X] -= 2;

    // write the number of elements
    long nelements = segment.size();
    if (fwrite(&(grid_size[OR_Z]), sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
    if (fwrite(&(grid_size[OR_Y]), sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
    if (fwrite(&(grid_size[OR_X]), sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
    if (fwrite(&nelements, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }

    // repad the grid size for iterations
    grid_size[OR_Z] += 2;
    grid_size[OR_Y] += 2;
    grid_size[OR_X] += 2;

    // iterate throught the skeleton endpoints and trace back to the surface
    for (std::unordered_map<long, char>::iterator it = segment.begin(); it != segment.end(); ++it) {
        long current_index = it->first;
        
        // convert to cartesian coordinates
        long ix, iy, iz;
        IndexToIndices(current_index, ix, iy, iz);

        while (erosions[current_index] != -1) {
            // get the direction to the next node
            long direction = erosions[current_index];

            // update the current index
            current_index = current_index + n6_offsets[direction];
        }

        // convert to cartesian coordinates
        long iu, iv, iw;
        IndexToIndices(current_index, iu, iv, iw);

        double xdiff = resolution[OR_X] * (ix - iu);
        double ydiff = resolution[OR_Y] * (iy - iv);
        double zdiff = resolution[OR_Z] * (iz - iw);

        double distance = sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);

        // unpad the cartesian coordinates
        ix--; iy--; iz--;
        // find the corrected linear index
        long index = iz * (grid_size[OR_X] - 2) * (grid_size[OR_Y] - 2) + iy * (grid_size[OR_X] - 2) + ix;

        if (fwrite(&index, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
        if (fwrite(&distance, sizeof(double), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
    }

    // close the file
    fclose(wfp);

    double total_time = (double) (clock() - start_time) / CLOCKS_PER_SEC;

    char time_filename[4096];
    sprintf(time_filename, "timings/radii/%s-%06ld.time", prefix, label);

    FILE *tfp = fopen(time_filename, "wb");
    if (!tfp) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }

    // write the number of points and the total time to file
    if (fwrite(&nelements, sizeof(long), 1, tfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }
    if (fwrite(&total_time, sizeof(double), 1, tfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }

    // close file
    fclose(tfp);
}
