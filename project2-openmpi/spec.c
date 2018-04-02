#include <stdio.h>
#include <stdint.h>

#include "spec.h"
#include "logging.h"

void __debug_print_spec(int rank, const Spec *spec) {
    INFO("Process %d: Spec at address 0x%lx", rank, (uint64_t)spec);
    INFO("TimeSlots:              %d", spec->cs.time_slot);
    INFO("TimeStep:               %.6f", spec->cs.time_step);
    INFO("Horizon:                %d", spec->cs.horizon);
    INFO("GridSize:               %d", spec->gs.size);
    INFO("NumberOfSmallParticles: %d", spec->gs.small_num);
    INFO("SmallParticleMass:      %.6f", spec->gs.small_mass);
    INFO("SmallParticleRadius:    %.6f", spec->gs.small_rad);
    INFO("NumberOfLargeParticles: %d", spec->gs.large_num);
    for (int i = 0, num = spec->gs.large_num; i < num; i++)
        INFO("%.4f\t%.4f\t%.4f\t%.4f", spec->gs.large_ptc[i].rad, spec->gs.large_ptc[i].mass, 
                spec->gs.large_ptc[i].loc.x, spec->gs.large_ptc[i].loc.y);
}

void __debug_print_compspec(int rank, const CompSpec *cs) {
    INFO("Process %d: CompSpec at address 0x%lx", rank, (uint64_t)cs);
    INFO("time_slot: %d", cs->time_slot);
    INFO("time_step: %.6f", cs->time_step);
    INFO("horizon:   %d", cs->horizon);
}

void __debug_print_location(int rank, const Location *loc, uint32_t count) {
    INFO("Process %d: Location[%u] at address 0x%lx", rank, count, (uint64_t)loc);
    for (uint32_t i = 0; i < count; i++) {
        INFO("Location[%u].x: %.4f", i, loc[i].x);
        INFO("Location[%u].y: %.4f", i, loc[i].y);        
    }
}

void __debug_print_particle(int rank, const Particle *ptc, uint32_t count) {
    INFO("Process %d: Particle[%u] at address 0x%lx", rank, count, (uint64_t)ptc);
    for (uint32_t i = 0; i < count; i++) {
        INFO("Particle[%u].rad:   %.4f", i, ptc[i].rad);
        INFO("Particle[%u].mass:  %.4f", i, ptc[i].mass);
        INFO("Particle[%u].loc.x: %.4f", i, ptc[i].loc.x);        
        INFO("Particle[%u].loc.y: %.4f", i, ptc[i].loc.y);        
    }
}