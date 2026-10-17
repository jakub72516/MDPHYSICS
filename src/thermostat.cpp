#include "mdsim/thermostat.hpp"

#include "mdsim/energy.hpp"
#include "mdsim/transcendental.hpp"
#include "mdsim/units.hpp"

namespace mdsim {

void Thermostat::setup(Fix64 target_temp_K, Fix64 relax_fs) {
    tau = relax_fs;
    Fix64 kb = Fix64::fromDouble(units::KB);
    kbT = kb * target_temp_K;
    Q = kbT * (tau * tau);
}

void Thermostat::apply(ParticleData& p, Fix64 dt) {
    // Target kinetic energy per atom is (3/2) kB T; d(xi)/dt = (ke - target)/Q.
    Fix64 one_half = Fix64::fromDouble(0.5);
    Fix64 target = (kbT + kbT) + kbT * one_half; // (3/2) kbT
    Fix64 ke = meanKinetic(p);
    Fix64 dxi = (ke - target) / Q;
    xi = xi + dxi * dt;

    // Drag: v *= exp(-xi * dt).
    Fix64 s = fexp(-(xi * dt));
    int64_t sv = s.raw();
    for (int i = 0; i < p.n; ++i) {
        p.vx[i] = fmul(p.vx[i], sv);
        p.vy[i] = fmul(p.vy[i], sv);
        p.vz[i] = fmul(p.vz[i], sv);
    }
}

} // namespace mdsim
