// Deterministic cell list for periodic neighbor search.
//
// Particles are assigned to cubic cells and sorted by cell key using a stable
// LSD radix sort. Because the sort is stable and the input is ordered by
// particle index, the ordering within every cell (and hence the entire
// neighbor search) is bit-identical on every run and every platform.
#pragma once

#include <vector>

#include "mdsim/fixed.hpp"
#include "mdsim/particles.hpp"

namespace mdsim {

struct CellList {
    Fix64 box;                 // cubic box side L
    int nc = 0;                // cells per dimension
    int ncells = 0;            // nc^3

    std::vector<int> order;        // particle indices sorted by cell key
    std::vector<int> cell_start;   // size ncells+1 (CSR offsets into order)

    void build(const ParticleData& p, Fix64 L, Fix64 rc);

    // Particle indices of cell c: order[cell_start[c]] .. order[cell_start[c+1]].
    int cellBegin(int c) const { return cell_start[c]; }
    int cellEnd(int c) const { return cell_start[c + 1]; }

    // Cell coordinate of particle i (using already-wrapped position).
    static void cellCoord(int64_t x, int64_t y, int64_t z,
                          int64_t L, int nc, int& cx, int& cy, int& cz);
};

// Stable LSD radix sort of `order` by `keys` (unsigned 32-bit).
void radix_sort_u32(const std::vector<std::uint32_t>& keys, std::vector<int>& order);

} // namespace mdsim
