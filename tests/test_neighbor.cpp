// Tests for deterministic neighbor search (cell list + pair list).
#include <cstdio>
#include <vector>

#include "mdsim/fixed.hpp"
#include "mdsim/particles.hpp"
#include "mdsim/rng.hpp"
#include "mdsim/cell_list.hpp"
#include "mdsim/neighbor_list.hpp"
#include "mdsim/forces.hpp"

using namespace mdsim;

static int fails = 0;
#define CHECK(cond) do { if (!(cond)) { std::printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond); ++fails; } } while (0)

int main() {
    // Small deterministic system: 512 atoms on a grid.
    ParticleData p;
    p.resize(512);
    for (int i = 0; i < p.n; ++i) {
        p.type[i] = 0;
        p.q[i] = 0;
        p.m[i] = Fix64::fromDouble(39.948).raw();
        p.inv_m[i] = fdiv(Fix64::fromDouble(4.184e-4).raw(), p.m[i]);
    }
    Fix64 L = Fix64::fromDouble(20.0);
    placeGrid(p, L);

    std::vector<AtomType> types(1);
    types[0].sigma = Fix64::fromDouble(3.405);
    types[0].epsilon = Fix64::fromDouble(0.238);
    types[0].charge = Fix64::fromDouble(0.0);
    types[0].mass = Fix64::fromDouble(39.948);
    PairParams pp = makePairParams(types, Fix64::fromDouble(6.0));
    ForceParams fp;
    fp.cutoff = Fix64::fromDouble(6.0);
    fp.cutoff2 = fp.cutoff * fp.cutoff;
    fp.es = Electrostatics::None;

    // Cell list determinism: build twice, identical.
    CellList c1, c2;
    c1.build(p, L, fp.cutoff);
    c2.build(p, L, fp.cutoff);
    CHECK(c1.order == c2.order);
    CHECK(c1.cell_start == c2.cell_start);

    // Pair list determinism + ordering (i < j).
    NeighborList n1, n2;
    n1.build(p, c1, pp, fp);
    n2.build(p, c2, pp, fp);
    CHECK(n1.pairs.size() == n2.pairs.size());
    bool order_ok = true;
    for (std::size_t k = 0; k < n1.pairs.size(); ++k) {
        if (n1.pairs[k].i != n2.pairs[k].i || n1.pairs[k].j != n2.pairs[k].j) order_ok = false;
        if (n1.pairs[k].i == n1.pairs[k].j) order_ok = false;
    }
    CHECK(order_ok);
    CHECK(n1.pairs.size() > 0);

    // All pairs within cutoff (allow fp rounding slack).
    bool within = true;
    int64_t Lr = L.raw();
    int64_t rc2 = fp.cutoff2.raw() + (Fix64::ONE >> 8);
    for (const auto& pr : n1.pairs) {
        int64_t dx = p.x[pr.i] - p.x[pr.j]; minImage(dx, Lr);
        int64_t dy = p.y[pr.i] - p.y[pr.j]; minImage(dy, Lr);
        int64_t dz = p.z[pr.i] - p.z[pr.j]; minImage(dz, Lr);
        int64_t r2 = fadd(fadd(fmul(dx, dx), fmul(dy, dy)), fmul(dz, dz));
        if (r2 > rc2) within = false;
    }
    CHECK(within);

    if (fails == 0) std::printf("test_neighbor: ALL PASS\n");
    else std::printf("test_neighbor: %d FAILURES\n", fails);
    return fails == 0 ? 0 : 1;
}
