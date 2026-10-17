// Deterministic SHAKE/RATTLE rigid-bond constraints.
//
// Constraints are enforced by a fixed number of Gauss-Seidel sweeps over the
// bonds in a strict (construction) order, so results are bit-identical across
// runs regardless of thread count or platform.
#pragma once

#include <vector>

#include "mdsim/fixed.hpp"
#include "mdsim/particles.hpp"

namespace mdsim {

struct Constraint {
    int a;
    int b;
    Fix64 d;        // target bond length (A)
    Fix64 d2;       // d^2
};

struct Shake {
    std::vector<Constraint> bonds;
    int max_iterations = 64;
    Fix64 tolerance = Fix64::fromDouble(1e-10);

    // Position constraints (SHAKE). Positions must be inside the box.
    void apply(ParticleData& p, Fix64 box);

    // Velocity constraints (RATTLE).
    void applyVelocities(ParticleData& p, Fix64 box);
};

} // namespace mdsim
