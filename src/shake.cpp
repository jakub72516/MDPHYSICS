#include "mdsim/shake.hpp"

#include "mdsim/forces.hpp" // minImage

namespace mdsim {

void Shake::apply(ParticleData& p, Fix64 box) {
    const int64_t L = box.raw();
    const int64_t two = Fix64::ONE << 1;
    const int nb = static_cast<int>(bonds.size());

    for (int iter = 0; iter < max_iterations; ++iter) {
        bool all_ok = true;
        for (int k = 0; k < nb; ++k) {
            const Constraint& c = bonds[k];
            int a = c.a, b = c.b;

            int64_t dx = p.x[a] - p.x[b]; minImage(dx, L);
            int64_t dy = p.y[a] - p.y[b]; minImage(dy, L);
            int64_t dz = p.z[a] - p.z[b]; minImage(dz, L);

            int64_t r2 = fadd(fadd(fmul(dx, dx), fmul(dy, dy)), fmul(dz, dz));
            int64_t diff = r2 - c.d2.raw();
            int64_t adiff = diff < 0 ? -diff : diff;
            if (adiff < tolerance.raw()) continue;
            all_ok = false;

            int64_t S = p.inv_m[a] + p.inv_m[b];
            // g = (r2 - d2) / (2 * S * r2)
            int64_t denom = fmul(fmul(two, S), r2);
            int64_t g = fdiv(diff, denom);

            int64_t ga = fmul(p.inv_m[a], g);
            int64_t gb = fmul(p.inv_m[b], g);

            p.x[a] -= fmul(ga, dx);
            p.y[a] -= fmul(ga, dy);
            p.z[a] -= fmul(ga, dz);

            p.x[b] += fmul(gb, dx);
            p.y[b] += fmul(gb, dy);
            p.z[b] += fmul(gb, dz);
        }
        if (all_ok) break;
    }
}

void Shake::applyVelocities(ParticleData& p, Fix64 box) {
    const int64_t L = box.raw();
    const int nb = static_cast<int>(bonds.size());

    for (int iter = 0; iter < max_iterations; ++iter) {
        bool all_ok = true;
        for (int k = 0; k < nb; ++k) {
            const Constraint& c = bonds[k];
            int a = c.a, b = c.b;

            int64_t dx = p.x[a] - p.x[b]; minImage(dx, L);
            int64_t dy = p.y[a] - p.y[b]; minImage(dy, L);
            int64_t dz = p.z[a] - p.z[b]; minImage(dz, L);
            int64_t r2 = fadd(fadd(fmul(dx, dx), fmul(dy, dy)), fmul(dz, dz));

            int64_t vrx = p.vx[a] - p.vx[b];
            int64_t vry = p.vy[a] - p.vy[b];
            int64_t vrz = p.vz[a] - p.vz[b];

            int64_t dot = fadd(fadd(fmul(dx, vrx), fmul(dy, vry)), fmul(dz, vrz));
            int64_t adot = dot < 0 ? -dot : dot;
            if (adot < tolerance.raw()) continue;
            all_ok = false;

            int64_t S = p.inv_m[a] + p.inv_m[b];
            int64_t g = fdiv(dot, fmul(S, r2));

            int64_t ga = fmul(p.inv_m[a], g);
            int64_t gb = fmul(p.inv_m[b], g);

            p.vx[a] -= fmul(ga, dx);
            p.vy[a] -= fmul(ga, dy);
            p.vz[a] -= fmul(ga, dz);

            p.vx[b] += fmul(gb, dx);
            p.vy[b] += fmul(gb, dy);
            p.vz[b] += fmul(gb, dz);
        }
        if (all_ok) break;
    }
}

} // namespace mdsim
