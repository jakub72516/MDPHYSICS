#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "mdsim/simulation.hpp"

using namespace mdsim;

static void usage(const char* argv0) {
    std::printf(
        "Usage: %s [options]\n"
        "Deterministic molecular dynamics simulator.\n\n"
        "Options:\n"
        "  -N <atoms>       number of atoms (default 108)\n"
        "  -steps <n>       number of timesteps (default 1000)\n"
        "  -dt <fs>         timestep in fs (default 1.0)\n"
        "  -seed <int>      random seed (default 42)\n"
        "  -output <file>   checkpoint/output file\n"
        "  -traj <file>     XYZ trajectory output file\n"
        "  -system <name>   argon | nacl (default argon)\n"
        "  -density <rho>   density in atoms/A^3 (default 0.02)\n"
        "  -temp <K>        temperature in K (default 300)\n"
        "  -cutoff <A>      pair cutoff in A (default 10)\n"
        "  -es <mode>       electrostatics: 0 none, 1 cutoff, 2 ewald (default 0)\n"
        "  -thermostat      enable Nose-Hoover thermostat\n"
        "  -barostat        enable isotropic barostat\n"
        "  -threads <n>     worker threads (default: hardware)\n"
        "  -print <n>       energy report interval (default 10)\n"
        "  -h, --help       show this help\n",
        argv0);
}

int main(int argc, char** argv) {
    SimConfig cfg;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&](const char* name) -> const char* {
            if (i + 1 < argc) return argv[++i];
            std::fprintf(stderr, "error: missing value for %s\n", name);
            std::exit(2);
        };

        if (a == "-N") cfg.n_atoms = std::atoi(next("-N"));
        else if (a == "-steps") cfg.steps = std::atoi(next("-steps"));
        else if (a == "-dt") cfg.dt_fs = std::atof(next("-dt"));
        else if (a == "-seed") cfg.seed = std::strtoull(next("-seed"), nullptr, 10);
        else if (a == "-output") cfg.output = next("-output");
        else if (a == "-traj") cfg.traj = next("-traj");
        else if (a == "-system") cfg.system = next("-system");
        else if (a == "-density") cfg.density = std::atof(next("-density"));
        else if (a == "-temp") cfg.temperature_K = std::atof(next("-temp"));
        else if (a == "-cutoff") cfg.cutoff_A = std::atof(next("-cutoff"));
        else if (a == "-es") cfg.electrostatics = std::atoi(next("-es"));
        else if (a == "-thermostat") cfg.thermostat = true;
        else if (a == "-barostat") cfg.barostat = true;
        else if (a == "-threads") cfg.threads = std::atoi(next("-threads"));
        else if (a == "-print") cfg.print_every = std::atoi(next("-print"));
        else if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
        else {
            std::fprintf(stderr, "error: unknown option '%s'\n", a.c_str());
            usage(argv[0]);
            return 2;
        }
    }

    SimResult r = runSimulation(cfg);

    std::printf("simulation complete\n");
    std::printf("  atoms            : %d\n", cfg.n_atoms);
    std::printf("  steps            : %d\n", cfg.steps);
    std::printf("  dt               : %.4f fs\n", cfg.dt_fs);
    std::printf("  seed             : %llu\n", (unsigned long long)cfg.seed);
    std::printf("  energy initial   : %.8f kcal/mol\n", r.energy_initial);
    std::printf("  energy final     : %.8f kcal/mol\n", r.energy_final);
    std::printf("  energy drift     : %.6e\n", r.drift);
    std::printf("  temperature      : %.3f K\n", r.temperature);
    std::printf("  MSD              : %.6f A^2\n", r.msd);
    std::printf("  performance      : %.3f M atom-steps/s\n", r.msteps_per_sec);
    std::printf("  SHA-256 state    : %s\n", r.sha256.c_str());

    return 0;
}
