// Deterministic Ewald summation for long-range electrostatics.
//
// Real-space (erfc) terms are handled by the short-range force loop; this
// module adds the reciprocal-space force/energy and the self term. The
// reciprocal-space sum runs over a fixed, ordered set of k-vectors and uses
// the fixed-point sin/cos (CORDIC) and exp routines, so it is bit-identical
// everywhere. For scalability beyond a few hundred charges, a PME/FFT back-end
// would replace the direct k-sum (see docs/design.md).
#pragma once

#include "mdsim/fixed.hpp"
#include "mdsim/particles.hpp"
#include "mdsim/u128.hpp"

namespace mdsim {

struct EwaldParams {
    Fix64 alpha;      // splitting parameter (1/A)
    Fix64 coulomb;    // C (kcal/(mol*A) per e^2)
    int kmax = 0;     // max |n| in reciprocal sum
    double tol = 1e-6;

    void setup(Fix64 alpha0, Fix64 box, Fix64 coulomb0);
};

// Add reciprocal-space forces to p.fx/fy/fz; return U_rec (raw, I128).
I128 ewaldReciprocal(ParticleData& p, Fix64 box, const EwaldParams& ep);

// Self term U_self = -C*alpha/sqrt(pi) * sum q_i^2 (raw, I128).
I128 ewaldSelf(const ParticleData& p, const EwaldParams& ep);

} // namespace mdsim
