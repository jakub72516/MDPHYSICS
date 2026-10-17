// Deterministic Nose-Hoover thermostat (single chain, intensive quantities).
//
// The friction update works with the mean kinetic energy per atom so that
// extensive quantities never overflow the fixed-point range for large systems.
#pragma once

#include "mdsim/fixed.hpp"
#include "mdsim/particles.hpp"

namespace mdsim {

struct Thermostat {
    Fix64 kbT;     // kB * target temperature (kcal/mol)
    Fix64 Q;       // thermostat mass (per atom): kbT * tau^2
    Fix64 tau;     // relaxation time (fs)
    Fix64 xi = Fix64(0); // friction coefficient (1/fs)

    void setup(Fix64 target_temp_K, Fix64 relax_fs);

    // Apply the drag v *= exp(-xi * dt) and advance the friction variable.
    void apply(ParticleData& p, Fix64 dt);
};

} // namespace mdsim
