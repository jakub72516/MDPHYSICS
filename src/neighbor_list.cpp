#include "mdsim/neighbor_list.hpp"

namespace mdsim {

void NeighborList::build(const ParticleData& p, const CellList& cl,
                         const PairParams& pp, const ForceParams& fp) {
    const int nc = cl.nc;
    const int ncells = cl.ncells;
    const int64_t L = cl.box.raw();
    const int64_t rc2 = fp.cutoff2.raw();

    pairs.clear();
    pairs.reserve(static_cast<std::size_t>(p.n) * 32);

    // Half-partition by linear cell index: for each cell c, visit the full
    // 27-neighbor stencil but only emit pairs whose neighbor cell index is
    // >= c (equal cell handled with j > i). With nc >= 3 every neighbor cell
    // is reached via a unique wrapping offset, so each unordered pair is
    // emitted exactly once.
    for (int c = 0; c < ncells; ++c) {
        int cz = c / (nc * nc);
        int cy = (c / nc) % nc;
        int cx = c % nc;

        for (int dnz = -1; dnz <= 1; ++dnz) {
            int nz = cz + dnz; if (nz < 0) nz += nc; if (nz >= nc) nz -= nc;
            for (int dny = -1; dny <= 1; ++dny) {
                int ny = cy + dny; if (ny < 0) ny += nc; if (ny >= nc) ny -= nc;
                for (int dnx = -1; dnx <= 1; ++dnx) {
                    int nx = cx + dnx; if (nx < 0) nx += nc; if (nx >= nc) nx -= nc;
                    int cn = (nz * nc + ny) * nc + nx;
                    if (cn < c) continue;

                    if (cn == c) {
                        for (int ii = cl.cellBegin(c); ii < cl.cellEnd(c); ++ii) {
                            int i = cl.order[ii];
                            for (int jj = ii + 1; jj < cl.cellEnd(c); ++jj) {
                                int j = cl.order[jj];
                                int64_t dx = p.x[i] - p.x[j]; minImage(dx, L);
                                int64_t dy = p.y[i] - p.y[j]; minImage(dy, L);
                                int64_t dz = p.z[i] - p.z[j]; minImage(dz, L);
                                int64_t r2 = fadd(fadd(fmul(dx, dx), fmul(dy, dy)), fmul(dz, dz));
                                if (r2 <= 0) continue;
                                int k = p.type[i] * pp.ntypes + p.type[j];
                                if (r2 < pp.min_r2[k]) continue;
                                if (r2 >= rc2) continue;
                                pairs.push_back(NeighborPair{i, j});
                            }
                        }
                    } else {
                        for (int ii = cl.cellBegin(c); ii < cl.cellEnd(c); ++ii) {
                            int i = cl.order[ii];
                            for (int jj = cl.cellBegin(cn); jj < cl.cellEnd(cn); ++jj) {
                                int j = cl.order[jj];
                                int64_t dx = p.x[i] - p.x[j]; minImage(dx, L);
                                int64_t dy = p.y[i] - p.y[j]; minImage(dy, L);
                                int64_t dz = p.z[i] - p.z[j]; minImage(dz, L);
                                int64_t r2 = fadd(fadd(fmul(dx, dx), fmul(dy, dy)), fmul(dz, dz));
                                if (r2 <= 0) continue;
                                int k = p.type[i] * pp.ntypes + p.type[j];
                                if (r2 < pp.min_r2[k]) continue;
                                if (r2 >= rc2) continue;
                                pairs.push_back(NeighborPair{i, j});
                            }
                        }
                    }
                }
            }
        }
    }
}

} // namespace mdsim
