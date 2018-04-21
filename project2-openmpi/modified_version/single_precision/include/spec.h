#ifndef SPEC_H
#define SPEC_H

#include <stdint.h>
#include <stddef.h>
#ifndef POOL_SEQ
#include "mpi.h"
#endif

typedef float Radius;
typedef float Mass;

/**
 * Define the location of a particle on the board
 */
typedef struct {
    float x;
    float y;
} Location;

/**
 * Definitions for MPI DataType of Location
 */
#ifndef POOL_SEQ
static const int LOCATION_COUNT = 2;
static const int LOCATION_BLOCK_LENGTH[] = {1, 1};
static const MPI_Aint LOCATION_DISPLACE[] = {
    offsetof(Location, x),
    offsetof(Location, y)
};
static const MPI_Datatype LOCATION_ELEM_TYPES[] = {
    MPI_FLOAT, MPI_FLOAT
};
#endif

/**
 * Define the information for a particle on the pool
 */
typedef struct {
    Radius   rad;
    Mass     mass;
    Location loc;
} Particle;

/**
 * Definitions for MPI DataType of Particle
 */
#ifndef POOL_SEQ
static const int PARTICLE_COUNT = 4;
static const int PARTICLE_BLOCK_LENGTH[] = {1, 1, 1, 1};
static const MPI_Aint PARTICLE_DISPLACE[] = {
    offsetof(Particle, rad),
    offsetof(Particle, mass),
    offsetof(Particle, loc) + offsetof(Location, x),
    offsetof(Particle, loc) + offsetof(Location, y)
};
static const MPI_Datatype PARTICLE_ELEM_TYPES[] = {
    MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT
};
#endif

/**
 * Define the specification for a Grid
 */
typedef struct {
    uint32_t  size;
    uint32_t  small_num;
    Mass      small_mass;
    Radius    small_rad;
    uint32_t  large_num;
    Particle *large_ptc;
} GridSpec;

/**
 * Define the specification for the computation
 */
typedef struct {
    int   time_slot;
    float time_step;
    int   horizon;
} CompSpec;

/**
 * Definitions for MPI DataType of CompSpec
 */
#ifndef POOL_SEQ
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
#endif

/**
 * Define the specification for one run
 */
typedef struct {
    CompSpec cs;
    GridSpec gs;
} Spec;

#ifdef POOL_DEBUG

/**
 * print Spec to STDOUT (only for debug use)
 */
void __debug_print_spec(int rank, const Spec *);

/**
 * print CompSpcec to STDOUT (only for debug use)
 */
void __debug_print_compspec(int rank, const CompSpec *cs);

/**
 * print an array of Location to STDOUT (only for debug use)
 */
void __debug_print_location(int rank, const Location *loc, uint32_t count);

/**
 * print an array of Particle to STDOUT (only for debug use)
 */
void __debug_print_particle(int rank, const Particle *ptc, uint32_t count);

#endif

#endif
