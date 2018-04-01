#ifndef PFILE_H
#define PFILE_H

// #define PFILE_DEBUG

#include "spec.h"
#include "pool.h"

#define MAX_SPECFILE_LINE_SIZE 64

#define TIME_SLOT_OFFSET      11
#define TIME_STEP_OFFSET      10
#define HORIZON_OFFSET        9
#define GRID_SIZE_OFFSET      10
#define SMALL_PTC_NUM_OFFSET  24
#define SMALL_PTC_MASS_OFFSET 19
#define SAMLL_PTC_RAD_OFFSET  21
#define LARGE_PTC_NUM_OFFSET  24

int read_spec_from_file(const char *filename, Spec *pool);
int print2ppm(const Pool *pool, const char *path);

#endif