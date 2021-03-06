#ifndef PFILE_H
#define PFILE_H

#include "spec.h"
#include "pool.h"

#ifndef POOL_SEQ
#include "mympi.h"
#endif

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
 * Structure describing one pixel in ppm files in RGB format, 
 * including red, green and blue channel. 
 */
typedef struct {
    uint8_t r;
    uint8_t g;
    uint8_t b;
} RGBPixel;

/**
 * Definitions for MPI DataType of RGBPixel
 */
#ifndef POOL_SEQ
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
#endif


/**
 * Structure describing a ppm file, including size (width and height) of
 * the photo and pixels (stored row by row). 
 */
typedef struct {
    uint64_t size;
    RGBPixel *pixels;
} PPMFile;

#ifndef POOL_SEQ
static const int PPMFILE_COUNT = 1;
static const int PPMFILE_BLOCK_LENGTH[] = {1};
static const MPI_Aint PPMFILE_DISPLACE[] = {
    offsetof(PPMFile, size)
};
static const MPI_Datatype PPMFILE_ELEM_TYPES[] = {
    MPI_UINT64_T
};
#endif


/**
 * Read specification from a file in the format as shown in Hugh's PDF. 
 * spec needs to be allocated first but spec->gs.large_ptc will 
 * be allocated dynamically according to the input specification file. 
 * Returns non-zero numbers when failure. 
 */
int read_spec_from_file(const char *filename, Spec *spec);


/**
 * Output pool into a ppm file whose path is specified by path. 
 * path needs to be a directory. 
 * Returns non-zero numbers when failure. 
 */
int print2ppm(const Pool *pool, const char *path);


/**
 * Output ppm to file whose path is specified by path. 
 * path needs to be a directory. 
 * Returns non-zero numbers when failure. 
 */
int ppm2file(const PPMFile *ppm, const char *path);


/**
 * Convert pool into ppm format stored in memory. 
 * ppm needs to be allocated first but ppm->pixels will 
 * be allocated dynamically according to size of the pool. 
 * Returns non-zero numbers when failure. 
 */
int pool2ppm(const Pool *pool, PPMFile *ppm);

#endif