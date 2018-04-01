#include <stdio.h>

#include "spec.h"

void __debug_print_spec(const Spec *spec) {
    printf("TimeSlots: %d\n", spec->cs.time_slot);
    printf("TimeStep: %f\n", spec->cs.time_step);
    printf("Horizon: %d\n\n", spec->cs.horizon);
    printf("GridSize: %d\n", spec->gs.size);
    printf("NumberOfSmallParticles: %d\n", spec->gs.small_num);
    printf("SmallParticleMass: %f\n", spec->gs.small_mass);
    printf("SmallParticleRadius: %f\n", spec->gs.small_rad);
    printf("NumberOfLargeParticles: %d\n", spec->gs.large_num);
    for (int i = 0, num = spec->gs.large_num; i < num; i++)
        printf("%.4f %.4f %.4f %.4f\n", spec->gs.large_ptc[i].rad, spec->gs.large_ptc[i].mass, 
                spec->gs.large_ptc[i].loc.x, spec->gs.large_ptc[i].loc.y);
}