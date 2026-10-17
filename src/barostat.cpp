#include "mdsim/barostat.hpp"

#include "mdsim/transcendental.hpp"

namespace mdsim {

// 1 atm in internal pressure units (kcal/(mol*A^3)).
static const double ATM_SCALE = 1.4584e-5;

Fix64 Barostat::apply(ParticleData& p, Fix64 box, Fix64 dt, U128 ke_raw) {
    // Virial: W = (1/3) sum_i r_i . F_i  (exact 128-bit accumulation).
    I128 virial = i128_zero();
    for (int i = 0; i < p.n; ++i) {
        int64_t rdotf = fadd(fadd(fmul(p.x[i], p.fx[i]), fmul(p.y[i], p.fy[i])),
                             fmul(p.z[i], p.fz[i]));
        virial = i128_add_i64(virial, rdotf);
    }
    int64_t w_raw = i128_div_i64(virial, 3);

    // Pressure in double (deterministic: fixed operation order, no libm).
    double ke = i128_to_double(i128_from_u128(ke_raw));
    double w = static_cast<double>(w_raw) / 281474976710656.0;
    double vol = box.toDouble();
    vol = vol * vol * vol;
    double P = ((2.0 / 3.0) * ke + w) / vol;

    double Ptarget = target_pressure_atm * ATM_SCALE;
    double base = 1.0 - (dt.toDouble() / tau_fs) * (Ptarget - P);
    if (base < 0.1) base = 0.1;
    if (base > 1.9) base = 1.9;

    // mu = base^(1/3), computed with fixed-point exp/log (deterministic).
    Fix64 basef = Fix64::fromDouble(base);
    Fix64 three = Fix64::fromDouble(3.0);
    Fix64 mu = fexp(flog(basef) / three);

    // Scale positions and box.
    int64_t m = mu.raw();
    for (int i = 0; i < p.n; ++i) {
        int64_t nx = fmul(p.x[i], m);
        int64_t ny = fmul(p.y[i], m);
        int64_t nz = fmul(p.z[i], m);
        p.x[i] = nx;
        p.y[i] = ny;
        p.z[i] = nz;
    }
    return Fix64::fromRaw(fmul(box.raw(), m));
}

} // namespace mdsim
