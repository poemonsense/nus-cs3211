#ifndef SPEC_H
#define SPEC_H

#include <stdint.h>

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

typedef struct {
    CompSpec cs;
    GridSpec gs;
} Spec;

void __debug_print_spec(const Spec *);

#endif