#include "mdsim/particles.hpp"

#include <cmath>

#include "mdsim/energy.hpp"
#include "mdsim/rng.hpp"
#include "mdsim/units.hpp"

namespace mdsim {

// Place particles on a simple-cubic grid filling the box [0,L)^3.
// Returns the number of grid cells per side.
int placeGrid(ParticleData& p, Fix64 L) {
    int n = p.n;
    int nx = 1;
    while (nx * nx * nx < n) ++nx;
    int64_t Lr = L.raw();
    Fix64 spacing = L / Fix64::fromDouble(static_cast<double>(nx));

    int idx = 0;
    for (int iz = 0; iz < nx && idx < n; ++iz) {
        for (int iy = 0; iy < nx && idx < n; ++iy) {
            for (int ix = 0; ix < nx && idx < n; ++ix) {
                p.x[idx] = fmul(spacing.raw(), Fix64::fromDouble(ix + 0.5).raw());
                p.y[idx] = fmul(spacing.raw(), Fix64::fromDouble(iy + 0.5).raw());
                p.z[idx] = fmul(spacing.raw(), Fix64::fromDouble(iz + 0.5).raw());
                ++idx;
            }
        }
    }
    return nx;
}

// Maxwell-Boltzmann velocities (deterministic, seeded RNG).
void initVelocities(ParticleData& p, Rng& rng, Fix64 kbT) {
    for (int i = 0; i < p.n; ++i) {
        // per-component variance = kbT * accel_scale / m = kbT * inv_m.
        Fix64 scale = Fix64::fromRaw(fsqrt(fmul(kbT.raw(), p.inv_m[i])));
        p.vx[i] = fmul(gaussian(rng).raw(), scale.raw());
        p.vy[i] = fmul(gaussian(rng).raw(), scale.raw());
        p.vz[i] = fmul(gaussian(rng).raw(), scale.raw());
    }
}

// Remove net center-of-mass momentum.
void removeDrift(ParticleData& p) {
    I128 px = i128_zero(), py = i128_zero(), pz = i128_zero();
    I128 mtot = i128_zero();
    for (int i = 0; i < p.n; ++i) {
        px = i128_add_i64(px, fmul(p.m[i], p.vx[i]));
        py = i128_add_i64(py, fmul(p.m[i], p.vy[i]));
        pz = i128_add_i64(pz, fmul(p.m[i], p.vz[i]));
        mtot = i128_add_i64(mtot, p.m[i]);
    }
    std::int64_t mt = static_cast<std::int64_t>(mtot.lo);
    if (mt <= 0) return;
    std::int64_t vcx = i128_div_i64(px, mt);
    std::int64_t vcy = i128_div_i64(py, mt);
    std::int64_t vcz = i128_div_i64(pz, mt);
    for (int i = 0; i < p.n; ++i) {
        p.vx[i] -= vcx;
        p.vy[i] -= vcy;
        p.vz[i] -= vcz;
    }
}

// Rescale velocities to the exact target kinetic energy (3/2) N kB T.
void rescaleToTemperature(ParticleData& p, Fix64 kbT) {
    if (p.n == 0) return;
    U128 cur = sumKineticRaw(p);
    // target = (3/2) N kB T
    double cur_d = i128_to_double(I128{cur.lo, static_cast<std::int64_t>(cur.hi)}) / 281474976710656.0;
    double target_d = 1.5 * p.n * kbT.toDouble();
    double scale = std::sqrt(target_d / cur_d);
    int64_t s = Fix64::fromDouble(scale).raw();
    for (int i = 0; i < p.n; ++i) {
        p.vx[i] = fmul(p.vx[i], s);
        p.vy[i] = fmul(p.vy[i], s);
        p.vz[i] = fmul(p.vz[i], s);
    }
}

} // namespace mdsim
