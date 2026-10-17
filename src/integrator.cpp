#include "mdsim/integrator.hpp"

namespace mdsim {

void Integrator::kickAndDrift(ParticleData& p, Fix64 box) {
    kick(p);
    drift(p, box);
}

void Integrator::drift(ParticleData& p, Fix64 box) {
    const int64_t dt = this->dt.raw();
    const int64_t L = box.raw();
    for (int i = 0; i < p.n; ++i) {
        int64_t nx = fadd(p.x[i], fmul(dt, p.vx[i]));
        int64_t ny = fadd(p.y[i], fmul(dt, p.vy[i]));
        int64_t nz = fadd(p.z[i], fmul(dt, p.vz[i]));

        while (nx < 0) nx += L;
        while (nx >= L) nx -= L;
        while (ny < 0) ny += L;
        while (ny >= L) ny -= L;
        while (nz < 0) nz += L;
        while (nz >= L) nz -= L;

        p.x[i] = nx;
        p.y[i] = ny;
        p.z[i] = nz;
    }
}

void Integrator::kick(ParticleData& p) {
    const int64_t hdt = half_dt.raw();
    for (int i = 0; i < p.n; ++i) {
        p.vx[i] = fadd(p.vx[i], fmul(hdt, fmul(p.fx[i], p.inv_m[i])));
        p.vy[i] = fadd(p.vy[i], fmul(hdt, fmul(p.fy[i], p.inv_m[i])));
        p.vz[i] = fadd(p.vz[i], fmul(hdt, fmul(p.fz[i], p.inv_m[i])));
    }
}

} // namespace mdsim
