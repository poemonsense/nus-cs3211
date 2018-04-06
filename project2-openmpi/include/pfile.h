#ifndef PFILE_H
#define PFILE_H

#include "spec.h"
#include "pool.h"
#include "mympi.h"

#define MAX_SPECFILE_LINE_SIZE 64

#define TIME_SLOT_OFFSET      11
#define TIME_STEP_OFFSET      10
#define HORIZON_OFFSET        9
#define GRID_SIZE_OFFSET      10
#define SMALL_PTC_NUM_OFFSET  24
#define SMALL_PTC_MASS_OFFSET 19
#define SAMLL_PTC_RAD_OFFSET  21
#define LARGE_PTC_NUM_OFFSET  24

/**
 * Definition of a pixel in PPM file in RGB format
 */
typedef struct {
    uint8_t r;
    uint8_t g;
    uint8_t b;
} RGBPixel;
/**
 * Definitions for MPI DataType of RGBPixel
 */
static const int RGBPIXEL_COUNT = 3;
static const int RGBPIXEL_BLOCK_LENGTH[] = {1, 1, 1};
static const MPI_Aint RGBPIXEL_DISPLACE[] = {
    offsetof(RGBPixel, r),
    offsetof(RGBPixel, g),
    offsetof(RGBPixel, b)
};
static const MPI_Datatype RGBPIXEL_ELEM_TYPES[] = {
    MPI_UINT8_T, MPI_UINT8_T, MPI_UINT8_T
};

typedef struct {
    uint64_t size;
    RGBPixel *pixels;
} PPMFile;
static const int PPMFILE_COUNT = 1;
static const int PPMFILE_BLOCK_LENGTH[] = {1};
static const MPI_Aint PPMFILE_DISPLACE[] = {
    offsetof(PPMFile, size)
};
static const MPI_Datatype PPMFILE_ELEM_TYPES[] = {
    MPI_UINT64_T
};

int read_spec_from_file(const char *filename, Spec *pool);
int print2ppm(const Pool *pool, const char *path);
int ppm2file(const PPMFile *ppm, const char *path);
int pool2ppm(const Pool *pool, PPMFile *ppm);

#endif