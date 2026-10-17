#include "mdsim/forces.hpp"

#include <cstdint>

#include "mdsim/transcendental.hpp"

namespace mdsim {

PairParams makePairParams(const std::vector<AtomType>& types, Fix64 cutoff) {
    PairParams pp;
    pp.ntypes = static_cast<int>(types.size());
    int nt = pp.ntypes;
    pp.sigma2.assign(nt * nt, 0);
    pp.eps.assign(nt * nt, 0);
    pp.min_r2.assign(nt * nt, 0);
    pp.shift.assign(nt * nt, 0);

    // Soft-core exclusion: 0.7 * sigma_ij (a hard repulsive floor that keeps
    // the r^-12 term inside the Fix64 range).
    Fix64 excl_frac = Fix64::fromDouble(0.7);
    Fix64 half = Fix64::fromDouble(0.5);
    Fix64 four(4 * Fix64::ONE);
    Fix64 cutoff2 = cutoff * cutoff;

    for (int i = 0; i < nt; ++i) {
        for (int j = 0; j < nt; ++j) {
            int k = i * nt + j;
            Fix64 sig = (types[i].sigma + types[j].sigma) * half; // arithmetic mean
            Fix64 eps = (types[i].epsilon * types[j].epsilon).sqrt(); // geometric
            pp.sigma2[k] = (sig * sig).raw();
            pp.eps[k] = eps.raw();
            Fix64 rmin = sig * excl_frac;
            pp.min_r2[k] = (rmin * rmin).raw();

            // LJ shift U(rc) so the potential vanishes at the cutoff.
            Fix64 s2 = (sig * sig) / cutoff2;
            Fix64 s6 = s2 * s2 * s2;
            Fix64 s12 = s6 * s6;
            pp.shift[k] = (four * eps * (s12 - s6)).raw();
        }
    }
    return pp;
}

void computeForces(ParticleData& p, const CellList& cl,
                   const PairParams& pp, const ForceParams& fp, ThreadPool* pool) {
    const int nc = cl.nc;
    const int ncells = cl.ncells;
    const int64_t L = cl.box.raw();
    const int64_t rc2 = fp.cutoff2.raw();
    const int64_t C = fp.coulomb.raw();

    const bool has_coul = (fp.es != Electrostatics::None);
    const bool ewald = (fp.es == Electrostatics::Ewald);
    const int64_t alpha = fp.ewald_alpha.raw();
    const int64_t twoa_sqrtpi = fp.two_alpha_over_sqrtpi.raw();

    const Fix64 twentyfour(24 * Fix64::ONE);

    int nworkers = pool ? pool->size() : 1;
    int nblocks = nworkers;
    if (nblocks > ncells) nblocks = ncells;
    if (nblocks < 1) nblocks = 1;

    // Determine min_r2 for a single-type system as a fast path.
    const std::vector<int64_t>& sigma2 = pp.sigma2;
    const std::vector<int64_t>& eps = pp.eps;
    const std::vector<int64_t>& minr2 = pp.min_r2;

    pool->runBlocks(nblocks, [&](int block) {
        int c0 = static_cast<int>((static_cast<long long>(block) * ncells) / nblocks);
        int c1 = static_cast<int>((static_cast<long long>(block + 1) * ncells) / nblocks);

        for (int c = c0; c < c1; ++c) {
            int cz = c / (nc * nc);
            int cy = (c / nc) % nc;
            int cx = c % nc;

            for (int ii = cl.cellBegin(c); ii < cl.cellEnd(c); ++ii) {
                int i = cl.order[ii];
                const int64_t xi = p.x[i], yi = p.y[i], zi = p.z[i];
                int ti = p.type[i];
                int64_t fxi = 0, fyi = 0, fzi = 0;

                for (int dnz = -1; dnz <= 1; ++dnz) {
                    int nz = cz + dnz; if (nz < 0) nz += nc; if (nz >= nc) nz -= nc;
                    for (int dny = -1; dny <= 1; ++dny) {
                        int ny = cy + dny; if (ny < 0) ny += nc; if (ny >= nc) ny -= nc;
                        for (int dnx = -1; dnx <= 1; ++dnx) {
                            int nx = cx + dnx; if (nx < 0) nx += nc; if (nx >= nc) nx -= nc;
                            int cn = (nz * nc + ny) * nc + nx;

                            for (int jj = cl.cellBegin(cn); jj < cl.cellEnd(cn); ++jj) {
                                int j = cl.order[jj];
                                if (j == i) continue;

                                int64_t dx = xi - p.x[j]; minImage(dx, L);
                                int64_t dy = yi - p.y[j]; minImage(dy, L);
                                int64_t dz = zi - p.z[j]; minImage(dz, L);

                                int64_t r2 = fadd(fadd(fmul(dx, dx), fmul(dy, dy)), fmul(dz, dz));
                                if (r2 <= 0) continue;

                                int tj = p.type[j];
                                int k = ti * pp.ntypes + tj;

                                if (r2 < minr2[k]) continue; // soft-core exclusion
                                if (r2 >= rc2) continue;

                                int64_t inv_r2 = finv(r2);

                                // Lennard-Jones.
                                int64_t s2 = fmul(sigma2[k], inv_r2);
                                int64_t s6 = fmul(s2, fmul(s2, s2));
                                int64_t s12 = fmul(s6, s6);
                                int64_t lj_factor = fmul(twentyfour.raw(),
                                    fmul(eps[k], fmul(fsub(fadd(s12, s12), s6), inv_r2)));

                                // Electrostatics (real space).
                                int64_t f = lj_factor;
                                if (has_coul) {
                                    int64_t qq = fmul(p.q[i], p.q[j]);
                                    if (qq != 0) {
                                        if (ewald) {
                                            // C q_i q_j [ erfc(a r) / r^3 + (2a/sqrt(pi)) e^{-a^2 r^2} / r^2 ]
                                            int64_t inv_r = fsqrt(inv_r2);
                                            int64_t r = fsqrt(r2);
                                            int64_t ar = fmul(alpha, r);
                                            int64_t erfc_ar = ferfc(Fix64(ar)).raw();
                                            int64_t a2r2 = fmul(ar, ar);
                                            int64_t eterm = fexp(Fix64(-a2r2)).raw();
                                            int64_t t1 = fmul(fmul(erfc_ar, inv_r2), inv_r);
                                            int64_t t2 = fmul(fmul(twoa_sqrtpi, eterm), inv_r2);
                                            f = fadd(f, fmul(fmul(C, qq), fadd(t1, t2)));
                                        } else {
                                            // C q_i q_j / r^3
                                            int64_t inv_r = fsqrt(inv_r2);
                                            int64_t t = fmul(fmul(inv_r2, inv_r), qq);
                                            f = fadd(f, fmul(C, t));
                                        }
                                    }
                                }

                                fxi += fmul(f, dx);
                                fyi += fmul(f, dy);
                                fzi += fmul(f, dz);
                            }
                        }
                    }
                }

                p.fx[i] = fxi;
                p.fy[i] = fyi;
                p.fz[i] = fzi;
            }
        }
    });
}

} // namespace mdsim
