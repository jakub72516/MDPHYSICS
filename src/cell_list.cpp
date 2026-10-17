#include "mdsim/cell_list.hpp"

#include <algorithm>
#include <cstdint>

namespace mdsim {

void radix_sort_u32(const std::vector<std::uint32_t>& keys, std::vector<int>& order) {
    const int n = static_cast<int>(keys.size());
    order.resize(n);
    for (int i = 0; i < n; ++i) order[i] = i;

    std::vector<int> tmp(n);
    for (int shift = 0; shift < 32; shift += 8) {
        int count[256] = {0};
        for (int i = 0; i < n; ++i) count[(keys[order[i]] >> shift) & 0xFF]++;
        int pos[256];
        pos[0] = 0;
        for (int b = 1; b < 256; ++b) pos[b] = pos[b - 1] + count[b - 1];
        for (int i = 0; i < n; ++i) {
            int b = (keys[order[i]] >> shift) & 0xFF;
            tmp[pos[b]++] = order[i];
        }
        order.swap(tmp);
    }
}

void CellList::cellCoord(int64_t x, int64_t y, int64_t z,
                         int64_t L, int nc, int& cx, int& cy, int& cz) {
    // ratio in [0,1): (x / L). floor(ratio * nc).
    auto coord = [&](int64_t v) -> int {
        int64_t ratio = fdiv(v, L); // [0, 2^48)
        int64_t c = (ratio * nc) >> Fix64::FRAC;
        if (c < 0) c = 0;
        if (c >= nc) c = nc - 1;
        return static_cast<int>(c);
    };
    cx = coord(x);
    cy = coord(y);
    cz = coord(z);
}

void CellList::build(const ParticleData& p, Fix64 L, Fix64 rc) {
    box = L;
    int64_t Lr = L.raw();
    // Number of cells per side: floor(L / rc), at least 1.
    nc = static_cast<int>(fdiv(Lr, rc.raw()) >> Fix64::FRAC);
    if (nc < 1) nc = 1;
    if (nc > 1024) nc = 1024; // guard against pathological tiny cutoffs
    ncells = nc * nc * nc;

    const int n = p.n;
    std::vector<std::uint32_t> keys(n);
    for (int i = 0; i < n; ++i) {
        int cx, cy, cz;
        cellCoord(p.x[i], p.y[i], p.z[i], Lr, nc, cx, cy, cz);
        std::uint32_t key = static_cast<std::uint32_t>((cz * nc + cy) * nc + cx);
        keys[i] = key;
    }

    std::vector<int> order(n);
    radix_sort_u32(keys, order);

    cell_start.assign(ncells + 1, 0);
    std::vector<int> cell_count(ncells, 0);
    for (int i = 0; i < n; ++i) cell_count[keys[order[i]]]++;

    int acc = 0;
    for (int c = 0; c < ncells; ++c) {
        cell_start[c] = acc;
        acc += cell_count[c];
    }
    cell_start[ncells] = acc;

    this->order = std::move(order);
}

} // namespace mdsim
