#ifndef POOL_H
#define POOL_H

#define POOL_DEBUG

#include <stdint.h>
#include "spec.h"
#include "mpi.h"

typedef struct {
    uint32_t  size;
    uint32_t  small_num;
    Mass      small_mass;
    Radius    small_rad;
    uint32_t  large_num;
    Location *small_ptc;
    Particle *large_ptc;
} Pool;

static const int POOL_COUNT = 5;
static const int POOL_BLOCK_LENGTH[] = {1, 1, 1, 1, 1, 1};
static const MPI_Aint POOL_DISPLACE[] = {
    offsetof(Pool, size),
    offsetof(Pool, small_num),
    offsetof(Pool, small_mass),
    offsetof(Pool, small_rad),
    offsetof(Pool, large_num)
};
static const MPI_Datatype POOL_ELEM_TYPES[] = {
    MPI_UINT32_T, MPI_UINT32_T, MPI_FLOAT, MPI_FLOAT, MPI_UINT32_T
};

int run(int rank, int size, int argc, char *argv[]);

void init_pool(int rank, const Spec *spec);
int init_params(int rank, int argc, char *argv[], Spec *spec);

/**
 * print Pool to STDOUT (only for debug use)
 */
void __debug_print_pool(int rank, const Pool *pool);

#endif