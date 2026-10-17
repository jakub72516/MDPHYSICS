// Reproducibility test: run a short 10 ps simulation and verify that the
// SHA-256 hash of the final state matches the reference value. The reference
// is computed from the deterministic fixed-point core, so it must match on
// any platform, compiler, or thread count.
#include <cstdio>
#include <string>

#include "mdsim/simulation.hpp"

using namespace mdsim;

// Reference SHA-256 (computed from the deterministic fixed-point core).
static const char* REFERENCE_HASH = "b5cd8805aabf433582077510164664773eaea822ff00fa1b8ebe309cdc426766";

int main() {
    SimConfig cfg;
    cfg.n_atoms = 108;
    cfg.steps = 1000;      // 10 ps at 1 fs
    cfg.dt_fs = 1.0;
    cfg.seed = 42;
    cfg.density = 0.02;
    cfg.cutoff_A = 8.0;
    cfg.system = "argon";
    cfg.electrostatics = 0;
    cfg.threads = 1;
    cfg.print_every = 0;

    SimResult r = runSimulation(cfg);

    std::printf("SHA-256: %s\n", r.sha256.c_str());
    std::printf("energy initial = %.8f  final = %.8f  drift = %.6e\n",
                r.energy_initial, r.energy_final, r.drift);

    if (std::string(REFERENCE_HASH).empty()) {
        std::printf("REFERENCE_HASH is empty; run completed but not verified.\n");
        return 0;
    }
    if (r.sha256 == REFERENCE_HASH) {
        std::printf("test_repro: PASS (hash matches reference)\n");
        return 0;
    }
    std::printf("test_repro: FAIL\n  got      %s\n  expected %s\n",
                r.sha256.c_str(), REFERENCE_HASH);
    return 1;
}
