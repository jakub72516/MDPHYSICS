// On-the-fly analysis: energy, mean-square displacement, radial distribution.
//
// All analysis arithmetic uses the same fixed-point core as the simulation,
// and pair-based quantities iterate the deterministic neighbor list, so the
// results are bit-identical across runs.
#pragma once

#include <vector>

#include "mdsim/fixed.hpp"
#include "mdsim/particles.hpp"
#include "mdsim/neighbor_list.hpp"
#include "mdsim/forces.hpp"
#include "mdsim/ewald.hpp"
#include "mdsim/u128.hpp"

namespace mdsim {

// Total potential energy (LJ + electrostatics + Ewald terms) as exact raw sum.
I128 potentialEnergy(const ParticleData& p, const NeighborList& nl,
                     const PairParams& pp, const ForceParams& fp,
                     const EwaldParams* ewald, Fix64 box);

struct Analysis {
    std::vector<std::int64_t> init_x, init_y, init_z;
    bool msd_ready = false;

    // RDF
    int rdf_bins = 200;
    Fix64 rdf_max;
    std::vector<std::int64_t> rdf_counts;
    std::int64_t rdf_pair_count = 0;

    void recordInitialPositions(const ParticleData& p);
    void initRdf(Fix64 max_dist);
    void accumulateRdf(const ParticleData& p, const NeighborList& nl, Fix64 box);
    // Mean-square displacement (per atom, minimum-image).
    Fix64 meanSquaredDisplacement(const ParticleData& p, Fix64 box) const;
};

} // namespace mdsim
