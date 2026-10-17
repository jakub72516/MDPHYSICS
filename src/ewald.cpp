#include "mdsim/ewald.hpp"

#include "mdsim/transcendental.hpp"

namespace mdsim {

void EwaldParams::setup(Fix64 alpha0, Fix64 box, Fix64 coulomb0) {
    alpha = alpha0;
    coulomb = coulomb0;
    // g_n = exp(-(pi/(alpha L))^2 * n^2); choose kmax so g < tol.
    Fix64 pi = Fix64::fromRaw(884279719003555LL);
    Fix64 pi_over_alphaL = pi / (alpha * box);
    Fix64 kfac = pi_over_alphaL * pi_over_alphaL;
    int64_t tol_raw = Fix64::fromDouble(tol).raw();
    int km = 1;
    while (km < 64) {
        Fix64 arg = Fix64::fromDouble(static_cast<double>(km * km)) * kfac;
        Fix64 g = fexp(-arg);
        if (g.raw() < tol_raw) break;
        ++km;
    }
    kmax = km;
}

I128 ewaldSelf(const ParticleData& p, const EwaldParams& ep) {
    // U_self = -(C alpha / sqrt(pi)) * sum q_i^2
    Fix64 pi = Fix64::fromRaw(884279719003555LL);
    Fix64 pref = -(ep.coulomb * ep.alpha) / pi.sqrt();
    I128 acc = i128_zero();
    for (int i = 0; i < p.n; ++i) {
        int64_t q2 = fmul(p.q[i], p.q[i]);
        acc = i128_add_i64(acc, fmul(pref.raw(), q2));
    }
    return acc;
}

I128 ewaldReciprocal(ParticleData& p, Fix64 box, const EwaldParams& ep) {
    // Enumerate half-space k-vectors n = (nx,ny,nz), 0 < n^2 <= kmax^2.
    // Formulas (tinfoil boundaries):
    //   U_rec = 2 * (C/(2 pi L)) * sum g_n/n^2 * (A_n^2 + B_n^2)
    //   F_i   = -2 * (2 C/L^2) q_i * sum g_n/n^2 * [B_n cos - A_n sin] * n
    // with g_n = exp(-(pi^2 n^2)/(alpha^2 L^2)),
    //      A_n = sum q_j cos(theta_j), B_n = sum q_j sin(theta_j),
    //      theta_j = (2 pi/L) (nx x_j + ny y_j + nz z_j).
    const int64_t Lr = box.raw();
    const Fix64 two_pi = Fix64::fromRaw(1768559438007110LL);
    Fix64 pi = Fix64::fromRaw(884279719003555LL);

    Fix64 two_pi_over_L = two_pi / box;
    Fix64 pi_over_alphaL = pi / (ep.alpha * box);
    Fix64 kfac = pi_over_alphaL * pi_over_alphaL;

    Fix64 C_over_2piL = ep.coulomb / (two_pi * box);   // C/(2 pi L)
    Fix64 twoC_over_L2 = (ep.coulomb + ep.coulomb) / (box * box); // 2C/L^2

    const int km = ep.kmax;
    I128 energy = i128_zero();

    std::vector<Fix64> theta(p.n);
    std::vector<Fix64> cosi(p.n), sini(p.n);

    for (int nx = 0; nx <= km; ++nx) {
        int ny_min = (nx == 0) ? 0 : -km;
        for (int ny = ny_min; ny <= km; ++ny) {
            int nz_min = (nx == 0 && ny == 0) ? 1 : -km;
            for (int nz = nz_min; nz <= km; ++nz) {
                int n2 = nx * nx + ny * ny + nz * nz;
                if (n2 > km * km) continue;

                Fix64 g_over_n2 = fexp(-(Fix64::fromDouble((double)n2) * kfac)) /
                                  Fix64::fromDouble((double)n2);

                // Structure factor components.
                I128 A = i128_zero();
                I128 B = i128_zero();
                for (int i = 0; i < p.n; ++i) {
                    Fix64 dot = Fix64::fromDouble((double)nx) * Fix64::fromRaw(p.x[i]) +
                                Fix64::fromDouble((double)ny) * Fix64::fromRaw(p.y[i]) +
                                Fix64::fromDouble((double)nz) * Fix64::fromRaw(p.z[i]);
                    theta[i] = two_pi_over_L * dot;
                    Fix64 s, c;
                    fsincos(theta[i], s, c);
                    cosi[i] = c; sini[i] = s;
                    A = i128_add_i64(A, fmul(p.q[i], c.raw()));
                    B = i128_add_i64(B, fmul(p.q[i], s.raw()));
                }

                // Energy term: 2 * C_over_2piL * g_over_n2 * (A^2 + B^2).
                int64_t a = static_cast<std::int64_t>(A.lo);
                int64_t b = static_cast<std::int64_t>(B.lo);
                int64_t a2 = fmul(a, a);
                int64_t b2 = fmul(b, b);
                int64_t ab = fadd(a2, b2);
                int64_t eterm = fmul(fmul(C_over_2piL.raw(), g_over_n2.raw()), ab);
                energy = i128_add_i64(energy, fadd(eterm, eterm)); // factor 2

                // Forces.
                Fix64 nxx = Fix64::fromDouble((double)nx);
                Fix64 nyy = Fix64::fromDouble((double)ny);
                Fix64 nzz = Fix64::fromDouble((double)nz);
                for (int i = 0; i < p.n; ++i) {
                    // bracket = B*cos - A*sin
                    int64_t br = fsub(fmul(b, cosi[i].raw()), fmul(a, sini[i].raw()));
                    int64_t coeff = fmul(fmul(twoC_over_L2.raw(), g_over_n2.raw()), fmul(p.q[i], br));
                    // factor 2 and negative sign folded into coeff
                    p.fx[i] -= fmul(fadd(coeff, coeff), nxx.raw());
                    p.fy[i] -= fmul(fadd(coeff, coeff), nyy.raw());
                    p.fz[i] -= fmul(fadd(coeff, coeff), nzz.raw());
                }
            }
        }
    }
    return energy;
}

} // namespace mdsim
