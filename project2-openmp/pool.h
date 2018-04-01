#ifndef POOL_H
#define POOL_H

#define POOL_DEBUG

#include "spec.h"
#include <stdint.h>

typedef struct {
    uint32_t  size;
    uint32_t  small_num;
    Mass      small_mass;
    Radius    small_rad;
    Location *small_ptc;
    uint32_t  large_num;
    Particle *large_ptc;
} Pool;

int run(int rank, int size, int argc, char *argv[]);

void init_pool(int rank, const Spec *spec);
int init_params(int argc, char *argv[], Spec *spec);


#endif