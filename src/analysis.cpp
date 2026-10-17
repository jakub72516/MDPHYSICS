#include "mdsim/analysis.hpp"

#include "mdsim/transcendental.hpp"
#include "mdsim/forces.hpp" // minImage

namespace mdsim {

I128 potentialEnergy(const ParticleData& p, const NeighborList& nl,
                     const PairParams& pp, const ForceParams& fp,
                     const EwaldParams* ewald, Fix64 box) {
    const int64_t L = box.raw();
    const int64_t C = fp.coulomb.raw();
    const bool has_coul = (fp.es != Electrostatics::None);
    const bool ew = (fp.es == Electrostatics::Ewald);
    const int64_t alpha = fp.ewald_alpha.raw();
    const Fix64 four(4 * Fix64::ONE);

    I128 energy = i128_zero();

    for (const NeighborPair& pr : nl.pairs) {
        int i = pr.i, j = pr.j;
        int64_t dx = p.x[i] - p.x[j]; minImage(dx, L);
        int64_t dy = p.y[i] - p.y[j]; minImage(dy, L);
        int64_t dz = p.z[i] - p.z[j]; minImage(dz, L);
        int64_t r2 = fadd(fadd(fmul(dx, dx), fmul(dy, dy)), fmul(dz, dz));

        int k = p.type[i] * pp.ntypes + p.type[j];

        // LJ energy (same exclusion/cutoff as the force loop).
        if (r2 > 0 && r2 >= pp.min_r2[k] && r2 < fp.cutoff2.raw()) {
            int64_t inv_r2 = finv(r2);
            int64_t s2 = fmul(pp.sigma2[k], inv_r2);
            int64_t s6 = fmul(s2, fmul(s2, s2));
            int64_t s12 = fmul(s6, s6);
            int64_t u = fmul(four.raw(), fmul(pp.eps[k], fsub(s12, s6)));
            energy = i128_add_i64(energy, fsub(u, pp.shift[k]));

            if (has_coul) {
                int64_t qq = fmul(p.q[i], p.q[j]);
                if (qq != 0) {
                    int64_t inv_r = fsqrt(inv_r2);
                    if (ew) {
                        int64_t r = fsqrt(r2);
                        int64_t er = ferfc(Fix64(fmul(alpha, r))).raw();
                        int64_t ue = fmul(fmul(C, qq), fmul(er, inv_r));
                        energy = i128_add_i64(energy, ue);
                    } else {
                        int64_t uc = fmul(fmul(C, qq), inv_r);
                        energy = i128_add_i64(energy, uc);
                    }
                }
            }
        }
    }

    if (ew && ewald) {
        I128 self = ewaldSelf(p, *ewald);
        energy = i128_add(energy, self);
    }
    return energy;
}

void Analysis::recordInitialPositions(const ParticleData& p) {
    init_x.assign(p.x.begin(), p.x.end());
    init_y.assign(p.y.begin(), p.y.end());
    init_z.assign(p.z.begin(), p.z.end());
    msd_ready = true;
}

void Analysis::initRdf(Fix64 max_dist) {
    rdf_max = max_dist;
    rdf_counts.assign(rdf_bins, 0);
    rdf_pair_count = 0;
}

void Analysis::accumulateRdf(const ParticleData& p, const NeighborList& nl, Fix64 box) {
    const int64_t L = box.raw();
    int64_t rmax = rdf_max.raw();
    for (const NeighborPair& pr : nl.pairs) {
        int64_t dx = p.x[pr.i] - p.x[pr.j]; minImage(dx, L);
        int64_t dy = p.y[pr.i] - p.y[pr.j]; minImage(dy, L);
        int64_t dz = p.z[pr.i] - p.z[pr.j]; minImage(dz, L);
        int64_t r2 = fadd(fadd(fmul(dx, dx), fmul(dy, dy)), fmul(dz, dz));
        int64_t r = fsqrt(r2);
        if (r >= rmax) continue;
        int bin = static_cast<int>((r * rdf_bins) / rmax);
        if (bin < 0) bin = 0;
        if (bin >= rdf_bins) bin = rdf_bins - 1;
        rdf_counts[bin]++;
        rdf_pair_count++;
    }
}

Fix64 Analysis::meanSquaredDisplacement(const ParticleData& p, Fix64 box) const {
    if (!msd_ready || p.n == 0) return Fix64(0);
    const int64_t L = box.raw();
    I128 acc = i128_zero();
    for (int i = 0; i < p.n; ++i) {
        int64_t dx = p.x[i] - init_x[i]; minImage(dx, L);
        int64_t dy = p.y[i] - init_y[i]; minImage(dy, L);
        int64_t dz = p.z[i] - init_z[i]; minImage(dz, L);
        int64_t d2 = fadd(fadd(fmul(dx, dx), fmul(dy, dy)), fmul(dz, dz));
        acc = i128_add_i64(acc, d2);
    }
    std::uint64_t avg = udiv128_64(U128{acc.lo, static_cast<std::uint64_t>(acc.hi)},
                                   static_cast<std::uint64_t>(p.n));
    return Fix64::fromRaw(static_cast<std::int64_t>(avg));
}

} // namespace mdsim
