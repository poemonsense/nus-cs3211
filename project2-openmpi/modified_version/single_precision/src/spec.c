#include <stdio.h>
#include <stdint.h>

#include "spec.h"

#ifdef POOL_DEBUG
#include "logging.h"

void __debug_print_spec(int rank, const Spec *spec) {
    DEBUG("Process %d: Spec at address 0x%lx", rank, (uint64_t)spec);
    DEBUG("TimeSlots:              %d", spec->cs.time_slot);
    DEBUG("TimeStep:               %.6f", spec->cs.time_step);
    DEBUG("Horizon:                %d", spec->cs.horizon);
    DEBUG("GridSize:               %d", spec->gs.size);
    DEBUG("NumberOfSmallParticles: %d", spec->gs.small_num);
    DEBUG("SmallParticleMass:      %.6f", spec->gs.small_mass);
    DEBUG("SmallParticleRadius:    %.6f", spec->gs.small_rad);
    DEBUG("NumberOfLargeParticles: %d", spec->gs.large_num);
    for (int i = 0, num = spec->gs.large_num; i < num; i++)
        DEBUG("%.4f\t%.4f\t%.4f\t%.4f", spec->gs.large_ptc[i].rad, spec->gs.large_ptc[i].mass, 
                spec->gs.large_ptc[i].loc.x, spec->gs.large_ptc[i].loc.y);
}

void __debug_print_compspec(int rank, const CompSpec *cs) {
    DEBUG("Process %d: CompSpec at address 0x%lx", rank, (uint64_t)cs);
    DEBUG("time_slot: %d", cs->time_slot);
    DEBUG("time_step: %.6f", cs->time_step);
    DEBUG("horizon:   %d", cs->horizon);
}

void __debug_print_location(int rank, const Location *loc, uint32_t count) {
    DEBUG("Process %d: Location[%u] at address 0x%lx", rank, count, (uint64_t)loc);
    for (uint32_t i = 0; i < count; i++) {
        DEBUG("Location[%u].x: %.4f", i, loc[i].x);
        DEBUG("Location[%u].y: %.4f", i, loc[i].y);        
    }
}

void __debug_print_particle(int rank, const Particle *ptc, uint32_t count) {
    DEBUG("Process %d: Particle[%u] at address 0x%lx", rank, count, (uint64_t)ptc);
    for (uint32_t i = 0; i < count; i++) {
        DEBUG("Particle[%u].rad:   %.4f", i, ptc[i].rad);
        DEBUG("Particle[%u].mass:  %.4f", i, ptc[i].mass);
        DEBUG("Particle[%u].loc.x: %.4f", i, ptc[i].loc.x);        
        DEBUG("Particle[%u].loc.y: %.4f", i, ptc[i].loc.y);        
    }
}

#endif