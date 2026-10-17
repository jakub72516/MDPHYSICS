// Top-level simulation orchestration.
#pragma once

#include <cstdint>
#include <string>

namespace mdsim {

struct SimConfig {
    int n_atoms = 108;
    int steps = 1000;
    double dt_fs = 1.0;
    std::uint64_t seed = 42;
    double density = 0.02;         // atoms/A^3
    double temperature_K = 300.0;
    std::string output;            // checkpoint/output path
    std::string traj;              // XYZ trajectory output path (optional)
    int threads = 1;
    int electrostatics = 0;        // 0=none, 1=cutoff, 2=ewald
    bool thermostat = false;
    bool barostat = false;
    std::string system = "argon";  // "argon" | "nacl"
    double cutoff_A = 10.0;
    int print_every = 10;
};

struct SimResult {
    std::string sha256;
    double energy_initial = 0.0;
    double energy_final = 0.0;
    double drift = 0.0;            // relative energy drift
    double temperature = 0.0;
    double msd = 0.0;
    double msteps_per_sec = 0.0;   // million atom-steps / s
    int steps = 0;
};

SimResult runSimulation(const SimConfig& cfg);

} // namespace mdsim
