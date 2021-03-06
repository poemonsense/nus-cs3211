#ifndef POOL_H
#define POOL_H

// #define POOL_DEBUG
// #define POOL_SEQ

#include <stdint.h>
#include "spec.h"

#ifndef POOL_SEQ
#include "mpi.h"
#endif


/**
 * Structure describing the pool table, including size (also width and height), 
 * number of small particles, mass, radius of small particles, number of large particles, 
 * location of every small particle and information on every large particle. 
 */
typedef struct {
    /* size of the pool table       */
    uint32_t  size;
    /* number of small particles    */
    uint32_t  small_num;
    /* mass of small particles      */
    Mass      small_mass;
    /* radius of small particles    */
    Radius    small_rad;
    /* number of large particles    */
    uint32_t  large_num;
    /* locations of small particles */
    Location *small_ptc;
    /* info about large particles, including mass, radius, positions., position.y  */
    Particle *large_ptc;
} Pool;

#ifndef POOL_SEQ
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
    MPI_UINT32_T, MPI_UINT32_T, MPI_DOUBLE, MPI_DOUBLE, MPI_UINT32_T
};
#endif

/**
 * describing the Acceleration of a particle, which doesn't need to be 
 * transferred between processes
 */
typedef double Accel_t;
typedef struct {
    Accel_t ax;
    Accel_t ay;
} Accel;

/**
 * Macro for computing the acceleration (acc1 and acc2)
 * Not final result since G constant has not been applied
 * Need to ensure that no name crashes in the local scope
 * 
 * vec: Location
 * mass1, mass2: double
 * acc1, acc2: Accel
 */
/* #define COMPUTE_ACCEL_2PTC_RAW(vec, mass1, mass2, acc1, acc2) {    \
    double __accel_dist_sqr = vec.x * vec.x + vec.y * vec.y;  \
    double __accel_dist = sqrt(__accel_dist_sqr);             \
    double __accel_den = __accel_dist * __accel_dist_sqr;     \
    double __accel_const_1 = mass2 / __accel_den;             \
    double __accel_const_2 = mass1 / __accel_den;             \
    acc1.ax += __accel_const_1 * vec.x;                       \
    acc1.ay += __accel_const_1 * vec.y;                       \
    acc2.ax -= __accel_const_2 * vec.x;                       \
    acc2.ay -= __accel_const_2 * vec.y;                       \
    } */

#define COMPUTE_ACCEL_2PTC_RAW(vec, mass1, mass2, acc1, acc2) {    \
    double __accel_dist_sqr = vec.x * vec.x + vec.y * vec.y;  \
    double __accel_const_1 = mass2 / __accel_dist_sqr;        \
    double __accel_const_2 = mass1 / __accel_dist_sqr;        \
    acc1.ax += __accel_const_1 * vec.x;                       \
    acc1.ay += __accel_const_1 * vec.y;                       \
    acc2.ax -= __accel_const_2 * vec.x;                       \
    acc2.ay -= __accel_const_2 * vec.y;                       \
    }

/**
 * Macro for computing the acceleration that ptc2 put on ptc1
 * Not final result since G constant has not been applied
 * Need to ensure that no name crashes in the local scope
 * 
 * vec: Location
 * mass1, mass2: double
 * acc1, acc2: Accel
 */
/* #define COMPUTE_ACCEL_RAW(vec, mass2, acc1) {    \
    double __accel_dist_sqr = vec.x * vec.x + vec.y * vec.y;  \
    double __accel_dist = sqrt(__accel_dist_sqr);             \
    double __accel_den = __accel_dist * __accel_dist_sqr;     \
    double __accel_const_1 = mass2 / __accel_den;             \
    acc1.ax += __accel_const_1 * vec.x;                       \
    acc1.ay += __accel_const_1 * vec.y;                       \
    } */

#define COMPUTE_ACCEL_RAW(vec, mass2, acc1) {    \
    double __accel_dist_sqr = vec.x * vec.x + vec.y * vec.y;  \
    double __accel_const_1 = mass2 / __accel_dist_sqr;        \
    acc1.ax += __accel_const_1 * vec.x;                       \
    acc1.ay += __accel_const_1 * vec.y;                       \
    }

/**
 * Gravatational constant for acceleration computation
 * To make the effect more notable, use 6.67 instead of 6.67*(10^-11)
 */
#define __ACCEL_G_CONSTANT 6.67 //0.000000000066742

#define ACCEL_RAW_TO_FINAL(acc) {  \
    acc.ax *= __ACCEL_G_CONSTANT;  \
    acc.ay *= __ACCEL_G_CONSTANT;  \
    }

/**
 * Describing the Velocity of a particle, which doesn't need to be 
 * transferred between processes
 */
typedef double Vel_t;
typedef struct {
    Vel_t vx;
    Vel_t vy;
} Velocity;


/**
 * Running the simulation, which begins with reading specification files. 
 */
#ifdef POOL_SEQ
int run(int argc, char *argv[]);
#else
int run(int rank, int size, int argc, char *argv[]);
#endif

#ifdef POOL_DEBUG

#include <assert.h>

/**
 * print Pool to STDOUT (only for debug use)
 */
void __debug_print_pool(int rank, const Pool *pool);

#define __DEBUG_ASSERT(expr) \
    assert(expr);

#endif

#endif