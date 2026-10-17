// Deterministic velocity-Verlet integrator (fixed timestep, no adaptivity).
#pragma once

#include "mdsim/fixed.hpp"
#include "mdsim/particles.hpp"

namespace mdsim {

struct Integrator {
    Fix64 dt;          // timestep (fs)
    Fix64 half_dt;     // 0.5 * dt

    // First half-step: v += 0.5 dt a; then x += dt v (with periodic wrap).
    void kickAndDrift(ParticleData& p, Fix64 box);

    // Second half-step: v += 0.5 dt a (forces must be up to date).
    void kick(ParticleData& p);

    // Position update: x += dt v, with periodic wrap.
    void drift(ParticleData& p, Fix64 box);
};

} // namespace mdsim
