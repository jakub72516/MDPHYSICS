// Structure-of-arrays particle container and atom-type definitions.
//
// All dynamical state (positions, velocities, forces) is stored in a SoA
// layout of contiguous std::vector<std::int64_t> (Q16.48), which is cache
// friendly and SIMD-ready. Static per-particle data (type, charge, mass) is
// stored separately.
#pragma once

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

#include "mdsim/fixed.hpp"

namespace mdsim {

struct AtomType {
    Fix64 sigma;     // A
    Fix64 epsilon;   // kcal/mol
    Fix64 charge;    // e
    Fix64 mass;      // g/mol
    std::string symbol; // element symbol (for XYZ output)
    int id;

    AtomType() : id(0) {}
};

struct ParticleData {
    int n = 0;

    // Dynamic state (Q16.48).
    std::vector<std::int64_t> x, y, z;
    std::vector<std::int64_t> vx, vy, vz;
    std::vector<std::int64_t> fx, fy, fz;

    // Static per-particle attributes.
    std::vector<std::int32_t> type;   // index into type table
    std::vector<std::int64_t> q;      // charge (e)
    std::vector<std::int64_t> m;      // mass (g/mol)
    std::vector<std::int64_t> inv_m;  // accel_scale / mass (acceleration scale)

    void resize(int count) {
        n = count;
        x.resize(n); y.resize(n); z.resize(n);
        vx.resize(n); vy.resize(n); vz.resize(n);
        fx.resize(n); fy.resize(n); fz.resize(n);
        type.resize(n);
        q.resize(n);
        m.resize(n);
        inv_m.resize(n);
    }

    void zeroForces() {
        std::fill(fx.begin(), fx.end(), 0);
        std::fill(fy.begin(), fy.end(), 0);
        std::fill(fz.begin(), fz.end(), 0);
    }
};

class Rng;

// Setup helpers (defined in particles.cpp).
int placeGrid(ParticleData& p, Fix64 L);
void initVelocities(ParticleData& p, Rng& rng, Fix64 kbT);
void removeDrift(ParticleData& p);
void rescaleToTemperature(ParticleData& p, Fix64 kbT);

} // namespace mdsim
