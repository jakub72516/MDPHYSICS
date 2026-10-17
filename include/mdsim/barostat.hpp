// Deterministic isotropic barostat (Berendsen-style pressure control).
//
// The instantaneous pressure is formed from the exact integer state (virial
// + kinetic energy via 128-bit accumulators); the control expression and the
// cube-root scaling are evaluated with the deterministic fixed-point
// transcendental routines, so the whole loop is bit-identical across runs.
#pragma once

#include "mdsim/fixed.hpp"
#include "mdsim/particles.hpp"
#include "mdsim/u128.hpp"

namespace mdsim {

struct Barostat {
    double target_pressure_atm = 1.0;
    double tau_fs = 1000.0;      // barostat relaxation time

    // Compute pressure, rescale box + positions, and return the new box.
    Fix64 apply(ParticleData& p, Fix64 box, Fix64 dt, U128 ke_raw);
};

} // namespace mdsim
