// Short-range force evaluation: Lennard-Jones + real-space electrostatics.
//
// The pair loop uses the deterministic cell list (full-neighbor per owned
// cell), so each particle's force is summed in a fixed order independent of
// thread count. A short-range exclusion radius prevents fixed-point overflow
// of the repulsive r^-12 term (see docs/design.md).
#pragma once

#include <vector>

#include "mdsim/fixed.hpp"
#include "mdsim/particles.hpp"
#include "mdsim/cell_list.hpp"
#include "mdsim/thread_pool.hpp"

namespace mdsim {

enum class Electrostatics { None, Cutoff, Ewald };

struct PairParams {
    int ntypes = 0;
    std::vector<std::int64_t> sigma2;  // sigma_ij^2 (ntypes*ntypes)
    std::vector<std::int64_t> eps;     // epsilon_ij
    std::vector<std::int64_t> min_r2;  // soft-core exclusion radius squared
    std::vector<std::int64_t> shift;   // LJ shift U(rc) (ntypes*ntypes)

    int idx(int ti, int tj) const { return ti * ntypes + tj; }
};

struct ForceParams {
    Fix64 cutoff;
    Fix64 cutoff2;
    Fix64 coulomb = Fix64::fromDouble(332.0637);
    Electrostatics es = Electrostatics::None;

    // Ewald real-space parameters (used only when es == Ewald).
    Fix64 ewald_alpha;
    Fix64 two_alpha_over_sqrtpi; // 2*alpha/sqrt(pi)
};

// Build per-type-pair tables from a type list.
PairParams makePairParams(const std::vector<AtomType>& types, Fix64 cutoff);

// Compute short-range forces into p.fx/fy/fz. Parallelized deterministically
// over cell blocks when pool is non-null.
void computeForces(ParticleData& p, const CellList& cl,
                   const PairParams& pp, const ForceParams& fp, ThreadPool* pool);

// Minimum-image displacement along one axis (d is x_i - x_j, box side L).
inline void minImage(std::int64_t& d, std::int64_t L) {
    std::int64_t half = L >> 1;
    if (d > half) d -= L;
    else if (d < -half) d += L;
}

} // namespace mdsim
