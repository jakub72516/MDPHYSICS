// Energy accumulation helpers (deterministic, exact via 128-bit sums).
#pragma once

#include "mdsim/u128.hpp"
#include "mdsim/fixed.hpp"
#include "mdsim/particles.hpp"

namespace mdsim {

// Total kinetic energy as an exact Q16.48 raw sum (extensive; may exceed
// 64 bits for very large systems, hence U128).
//
// inv_m = accel_scale / m, so the kinetic energy in kcal/mol is
//   0.5 * v^2 / inv_m = 0.5 * m * v^2 / accel_scale,
// which is exactly the energy whose time derivative is F.v under the
// integrator a = F * inv_m (this keeps energy conservation exact).
inline U128 sumKineticRaw(const ParticleData& p) {
    U128 s{0, 0};
    for (int i = 0; i < p.n; ++i) {
        int64_t v2 = fadd(fadd(fmul(p.vx[i], p.vx[i]),
                               fmul(p.vy[i], p.vy[i])),
                          fmul(p.vz[i], p.vz[i]));
        int64_t half = Fix64::ONE >> 1;
        int64_t ke = fmul(half, fdiv(v2, p.inv_m[i]));
        s = uadd(s, U128{static_cast<std::uint64_t>(ke), 0});
    }
    return s;
}

// Mean kinetic energy per atom (intensive; fits in Fix64).
inline Fix64 meanKinetic(const ParticleData& p) {
    if (p.n == 0) return Fix64(0);
    U128 s = sumKineticRaw(p);
    std::uint64_t avg = udiv128_64(s, static_cast<std::uint64_t>(p.n));
    return Fix64::fromRaw(static_cast<std::int64_t>(avg));
}

} // namespace mdsim
