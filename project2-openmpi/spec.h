#ifndef SPEC_H
#define SPEC_H

#include <stdint.h>
#include <stddef.h>
#include "mpi.h"

typedef float Radius;
typedef float Mass;
typedef struct {
    float x;
    float y;
} Location;

typedef struct {
    Radius   rad;
    Mass     mass;
    Location loc;
} Particle;

typedef struct {
    uint32_t  size;
    uint32_t  small_num;
    Mass      small_mass;
    Radius    small_rad;
    uint32_t  large_num;
    Particle *large_ptc;
} GridSpec;

typedef struct {
    int      time_slot;
    float    time_step;
    int      horizon;
} CompSpec;

static const int COMPSPEC_COUNT = 3;
static const int COMPSPEC_BLOCK_LENGTH[] = {1, 1, 1};
static const MPI_Aint COMPSPEC_DISPLACE[] = {
    offsetof(CompSpec, time_slot),
    offsetof(CompSpec, time_step),
    offsetof(CompSpec, horizon)
};
static const MPI_Datatype COMPSPEC_ELEM_TYPES[] = {
    MPI_INT, MPI_FLOAT, MPI_INT
};

typedef struct {
    CompSpec cs;
    GridSpec gs;
} Spec;

void __debug_print_spec(const Spec *);

#endif