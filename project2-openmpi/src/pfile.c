#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "spec.h"
#include "pfile.h"

#ifdef POOL_DEBUG
#include "logging.h"
#endif

int read_spec_from_file(const char *filename, Spec *spec) {
    FILE *spec_file = fopen(filename, "r");
    if (!spec_file) {
        fprintf(stderr, "[%s:%d] Cannot open specfile %s\n", __FILE__, __LINE__, filename);
        return -1;
    }
    // buffer for reading every line of the file 
    char buf[MAX_SPECFILE_LINE_SIZE];
    if (!fgets(buf, MAX_SPECFILE_LINE_SIZE, spec_file))
        return -1;
    spec->cs.time_slot = atoi(buf + TIME_SLOT_OFFSET);
    if (!fgets(buf, MAX_SPECFILE_LINE_SIZE, spec_file))
        return -1;
    spec->cs.time_step = atof(buf + TIME_STEP_OFFSET);
    if (!fgets(buf, MAX_SPECFILE_LINE_SIZE, spec_file))
        return -1;
    spec->cs.horizon = atoi(buf + HORIZON_OFFSET);
    if (!fgets(buf, MAX_SPECFILE_LINE_SIZE, spec_file))
        return -1;
    spec->gs.size = atoi(buf + GRID_SIZE_OFFSET);
    if (!fgets(buf, MAX_SPECFILE_LINE_SIZE, spec_file))
        return -1;
    spec->gs.small_num = atoi(buf + SMALL_PTC_NUM_OFFSET);
    if (!fgets(buf, MAX_SPECFILE_LINE_SIZE, spec_file))
        return -1;
    spec->gs.small_mass = atof(buf + SMALL_PTC_MASS_OFFSET);
    if (!fgets(buf, MAX_SPECFILE_LINE_SIZE, spec_file))
        return -1;
    spec->gs.small_rad = atof(buf + SAMLL_PTC_RAD_OFFSET);
    if (!fgets(buf, MAX_SPECFILE_LINE_SIZE, spec_file))
        return -1;
    spec->gs.large_num = atoi(buf + LARGE_PTC_NUM_OFFSET);
    spec->gs.large_ptc = (Particle *)calloc(spec->gs.large_num, sizeof(Particle));
    for (int i = 0; i < spec->gs.large_num; i++) {
        if (!fgets(buf, MAX_SPECFILE_LINE_SIZE, spec_file))
            return -1;
        sscanf(buf, "%lf %lf %lf %lf", &spec->gs.large_ptc[i].rad, &spec->gs.large_ptc[i].mass, 
                &spec->gs.large_ptc[i].loc.x, &spec->gs.large_ptc[i].loc.y);
    }
    #ifdef POOL_DEBUG
    __debug_print_spec(-1, spec);
    #endif
    return 0;
}

int print2ppm(const Pool *pool, const char *path) {
    #ifdef POOL_DEBUG
    INFO("Writing Pool at 0x%lx to ppm file %s", (uint64_t)pool, path)
    #endif
    FILE *ppm_file = fopen(path, "w");
    if (!ppm_file) {
        fprintf(stderr, "[%s:%d] Cannot open ppmfile %s\n", __FILE__, __LINE__, path);
        return -1;
    }
    uint32_t entry_num = pool->size * pool->size;
    // compute R channel
    uint8_t *red = (uint8_t *)calloc(entry_num, sizeof(uint8_t));
    memset(red, 0, entry_num * sizeof(uint8_t));
    for (uint32_t i = 0, num = pool->small_num; i < num; i++) {
        uint32_t col = (uint32_t)pool->small_ptc[i].x;
        uint32_t row = (uint32_t)pool->small_ptc[i].y;
        uint32_t pos = pool->size * row + col;
        if (red[pos] != 255)
            red[pos]++;
    }
    // compute B channel
    uint8_t *blue = (uint8_t *)calloc(entry_num, sizeof(uint8_t));
    memset(blue, 0, entry_num * sizeof(uint8_t));
    for (uint32_t i = 0, num = pool->large_num; i < num; i++) {
        // in case that large ball has gone out of the region
        uint32_t xrange[2] = { 0, pool->size - 1 };
        if (pool->large_ptc[i].loc.x > pool->large_ptc[i].rad)
            xrange[0] = (uint32_t)(pool->large_ptc[i].loc.x - pool->large_ptc[i].rad);
        if (pool->large_ptc[i].loc.x + pool->large_ptc[i].rad < pool->size - 1)
            xrange[1] = (uint32_t)(pool->large_ptc[i].loc.x + pool->large_ptc[i].rad);
        uint32_t yrange[2] = { 0, pool->size - 1 };
        if (pool->large_ptc[i].loc.y > pool->large_ptc[i].rad)
            yrange[0] = (uint32_t)(pool->large_ptc[i].loc.y - pool->large_ptc[i].rad);
        if (pool->large_ptc[i].loc.y + pool->large_ptc[i].rad < pool->size - 1)
            yrange[1] = (uint32_t)(pool->large_ptc[i].loc.y + pool->large_ptc[i].rad);
        for (uint32_t row = yrange[0]; row <= yrange[1]; row++) {
            for (uint32_t col = xrange[0]; col <= xrange[1]; col++) {
                double dis = (row - pool->large_ptc[i].loc.y) * (row - pool->large_ptc[i].loc.y)
                    + (col - pool->large_ptc[i].loc.x) * (col - pool->large_ptc[i].loc.x);
                if (dis <= pool->large_ptc[i].rad)
                    blue[row * pool->size + col] = 255;
            }
        }
    }
    // write to ppm file
    fprintf(ppm_file, "P3\n%d %d\n255\n", pool->size, pool->size);
    for (uint32_t row = 0; row < pool->size; row++) {
        for (uint32_t col = 0; col < pool->size; col++) {
            uint32_t idx = row * pool->size + col;
            if (blue[idx])
                fprintf(ppm_file, "%d %d %d\t", 0, 0, blue[idx]);
            else
                fprintf(ppm_file, "%d %d %d\t", (red[idx] > 0) ? 255: 0, 0, 0);
        }
        fprintf(ppm_file, "\n");
    }
    free(red);
    free(blue);
    fclose(ppm_file);
    return 0;
}