#include "mdsim/rng.hpp"

#include "mdsim/transcendental.hpp"

namespace mdsim {

Fix64 gaussian(Rng& rng) {
    // Box-Muller transform:
    //   z = sqrt(-2 ln u1) * cos(2 pi u2)
    // with u1 in (0, 1] and u2 in [0, 1). All arithmetic is fixed point.
    Fix64 u = rng.uniform();                 // [0, 1)
    Fix64 u1 = Fix64::fromRaw(Fix64::ONE - u.raw()); // (0, 1]
    Fix64 u2 = rng.uniform();                // [0, 1)

    Fix64 lnu1 = flog(u1);
    Fix64 neg_two = Fix64::fromRaw(-(Fix64::ONE << 1));
    Fix64 s = (neg_two * lnu1).sqrt();

    Fix64 two_pi = Fix64::fromRaw(1768559438007110LL);
    Fix64 ang = two_pi * u2;
    Fix64 c = fcos(ang);
    return s * c;
}

} // namespace mdsim
