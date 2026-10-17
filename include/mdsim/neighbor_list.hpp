// Deterministic neighbor-pair list (unique pairs i < j within cutoff).
//
// Pairs are enumerated with a fixed half-stencil over the sorted cell list,
// so the ordering (and the set) of pairs is bit-identical across runs. Used
// for potential-energy evaluation and radial distribution analysis.
#pragma once

#include <vector>

#include "mdsim/particles.hpp"
#include "mdsim/cell_list.hpp"
#include "mdsim/forces.hpp"

namespace mdsim {

struct NeighborPair {
    int i;
    int j;
};

struct NeighborList {
    std::vector<NeighborPair> pairs;

    void build(const ParticleData& p, const CellList& cl,
               const PairParams& pp, const ForceParams& fp);
};

} // namespace mdsim
