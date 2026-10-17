#include "mdsim/simulation.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#include "mdsim/fixed.hpp"
#include "mdsim/particles.hpp"
#include "mdsim/rng.hpp"
#include "mdsim/transcendental.hpp"
#include "mdsim/cell_list.hpp"
#include "mdsim/neighbor_list.hpp"
#include "mdsim/forces.hpp"
#include "mdsim/ewald.hpp"
#include "mdsim/integrator.hpp"
#include "mdsim/thermostat.hpp"
#include "mdsim/barostat.hpp"
#include "mdsim/shake.hpp"
#include "mdsim/thread_pool.hpp"
#include "mdsim/checkpoint.hpp"
#include "mdsim/sha256.hpp"
#include "mdsim/analysis.hpp"
#include "mdsim/energy.hpp"
#include "mdsim/units.hpp"

namespace mdsim {

namespace {

// Deterministic cubic box side: L = exp( (ln N - ln density) / 3 ).
Fix64 boxSide(int n_atoms, double density) {
    Fix64 lnN = flog(Fix64::fromRaw(static_cast<std::int64_t>(n_atoms)));
    Fix64 ln2 = Fix64::fromRaw(195103586505167LL);
    lnN = lnN + Fix64::fromRaw(48 * ln2.raw());
    Fix64 lnV = lnN - flog(Fix64::fromDouble(density));
    return fexp(lnV / Fix64::fromDouble(3.0));
}

double toDouble(I128 v) { return i128_to_double(v) / 281474976710656.0; }

} // namespace

SimResult runSimulation(const SimConfig& cfg) {
    SimResult result;
    result.steps = cfg.steps;

    // --- Atom types -------------------------------------------------------
    std::vector<AtomType> types;
    if (cfg.system == "nacl") {
        AtomType na, cl;
        na.sigma = Fix64::fromDouble(2.583); na.epsilon = Fix64::fromDouble(0.10);
        na.charge = Fix64::fromDouble(1.0);  na.mass = Fix64::fromDouble(22.989769);
        na.id = 0; na.symbol = "Na";
        cl.sigma = Fix64::fromDouble(4.404); cl.epsilon = Fix64::fromDouble(0.10);
        cl.charge = Fix64::fromDouble(-1.0); cl.mass = Fix64::fromDouble(35.453);
        cl.id = 1; cl.symbol = "Cl";
        types = {na, cl};
    } else {
        AtomType ar;
        ar.sigma = Fix64::fromDouble(3.405); ar.epsilon = Fix64::fromDouble(0.238);
        ar.charge = Fix64::fromDouble(0.0);  ar.mass = Fix64::fromDouble(39.948);
        ar.id = 0; ar.symbol = "Ar";
        types = {ar};
    }

    // --- Geometry & timestep ---------------------------------------------
    Fix64 box = boxSide(cfg.n_atoms, cfg.density);
    double rc = cfg.cutoff_A;
    // The cell list requires at least 3 cells per dimension (nc >= 3) so that
    // every neighbor cell is reached by a unique periodic offset; enforce
    // cutoff <= L/3.
    if (rc > box.toDouble() / 3.0) rc = box.toDouble() / 3.0;
    Fix64 cutoff = Fix64::fromDouble(rc);
    Fix64 dt = Fix64::fromDouble(cfg.dt_fs);

    // --- Particles --------------------------------------------------------
    ParticleData p;
    p.resize(cfg.n_atoms);
    bool charged = (cfg.electrostatics != 0);
    for (int i = 0; i < cfg.n_atoms; ++i) {
        int t = (cfg.system == "nacl") ? (i & 1) : 0;
        p.type[i] = t;
        p.q[i] = charged ? types[t].charge.raw() : 0;
        p.m[i] = types[t].mass.raw();
        Fix64 accel = Fix64::fromDouble(units::ACCEL_SCALE);
        p.inv_m[i] = fdiv(accel.raw(), p.m[i]);
    }
    placeGrid(p, box);

    // --- Initial velocities ----------------------------------------------
    Rng rng(cfg.seed);
    Fix64 kbT = Fix64::fromDouble(units::KB) * Fix64::fromDouble(cfg.temperature_K);
    initVelocities(p, rng, kbT);
    removeDrift(p);
    rescaleToTemperature(p, kbT);

    // --- Force tables -----------------------------------------------------
    PairParams pp = makePairParams(types, cutoff);
    ForceParams fp;
    fp.cutoff = cutoff;
    fp.cutoff2 = cutoff * cutoff;
    fp.coulomb = Fix64::fromDouble(units::COULOMB);
    fp.es = (cfg.electrostatics == 2) ? Electrostatics::Ewald
          : (cfg.electrostatics == 1) ? Electrostatics::Cutoff
          : Electrostatics::None;

    EwaldParams ewald;
    if (fp.es == Electrostatics::Ewald) {
        Fix64 alpha = Fix64::fromDouble(0.3);
        ewald.setup(alpha, box, fp.coulomb);
        fp.ewald_alpha = alpha;
        Fix64 pi = Fix64::fromRaw(884279719003555LL);
        fp.two_alpha_over_sqrtpi = (alpha + alpha) / pi.sqrt();
    }

    // --- Runtime ----------------------------------------------------------
    int nthreads = cfg.threads;
    if (nthreads < 1) nthreads = hardwareConcurrency();
    ThreadPool pool(nthreads);

    Integrator integ;
    integ.dt = dt;
    integ.half_dt = dt * Fix64::fromDouble(0.5);

    Thermostat thermo;
    if (cfg.thermostat) thermo.setup(Fix64::fromDouble(cfg.temperature_K), Fix64::fromDouble(100.0));

    Barostat baro;
    baro.target_pressure_atm = 1.0;
    baro.tau_fs = 1000.0;

    Analysis analysis;
    analysis.recordInitialPositions(p);
    analysis.initRdf(cutoff * Fix64::fromDouble(0.5));

    std::vector<std::string> symbols;
    for (const AtomType& t : types) symbols.push_back(t.symbol);
    if (!cfg.traj.empty()) {
        std::FILE* tf = std::fopen(cfg.traj.c_str(), "w");
        if (tf) std::fclose(tf);
        appendXyzFrame(cfg.traj, p, symbols, box, 0);
    }

    // --- Initial forces ---------------------------------------------------
    CellList cl;
    cl.build(p, box, cutoff);
    computeForces(p, cl, pp, fp, &pool);
    if (fp.es == Electrostatics::Ewald) ewaldReciprocal(p, box, ewald);

    // --- Initial energy ---------------------------------------------------
    NeighborList nl;
    nl.build(p, cl, pp, fp);
    I128 pe0 = potentialEnergy(p, nl, pp, fp, (fp.es == Electrostatics::Ewald) ? &ewald : nullptr, box);
    if (fp.es == Electrostatics::Ewald) pe0 = i128_add(pe0, ewaldSelf(p, ewald));
    double ke0 = i128_to_double(i128_from_u128(sumKineticRaw(p))) / 281474976710656.0;
    double e0 = ke0 + toDouble(pe0);
    result.energy_initial = e0;

    // --- Main loop --------------------------------------------------------
    auto t_start = std::chrono::steady_clock::now();
    for (int step = 1; step <= cfg.steps; ++step) {
        integ.kick(p);

        if (cfg.thermostat) thermo.apply(p, dt);

        integ.drift(p, box);

        cl.build(p, box, cutoff);
        computeForces(p, cl, pp, fp, &pool);
        if (fp.es == Electrostatics::Ewald) ewaldReciprocal(p, box, ewald);

        integ.kick(p);

        if (cfg.barostat) {
            box = baro.apply(p, box, dt, sumKineticRaw(p));
        }

        if (cfg.print_every > 0 && (step % cfg.print_every) == 0) {
            nl.build(p, cl, pp, fp);
            if (!cfg.traj.empty()) appendXyzFrame(cfg.traj, p, symbols, box, step);
            I128 pe = potentialEnergy(p, nl, pp, fp, (fp.es == Electrostatics::Ewald) ? &ewald : nullptr, box);
            if (fp.es == Electrostatics::Ewald) pe = i128_add(pe, ewaldSelf(p, ewald));
            double ke = i128_to_double(i128_from_u128(sumKineticRaw(p))) / 281474976710656.0;
            double E = ke + toDouble(pe);
            double T = (2.0 / 3.0) * meanKinetic(p).toDouble() / units::KB;
            if (step == cfg.steps) {
                result.energy_final = E;
                result.temperature = T;
                result.drift = (E - e0) / (e0 == 0.0 ? 1.0 : e0);
            }
        }
    }
    auto t_end = std::chrono::steady_clock::now();
    double seconds = std::chrono::duration<double>(t_end - t_start).count();
    result.msteps_per_sec = (double(cfg.n_atoms) * cfg.steps) / seconds / 1e6;

    // --- Final energy (if not already captured) ---------------------------
    if (cfg.print_every <= 0 || (cfg.steps % cfg.print_every) != 0) {
        nl.build(p, cl, pp, fp);
        I128 pe = potentialEnergy(p, nl, pp, fp, (fp.es == Electrostatics::Ewald) ? &ewald : nullptr, box);
        if (fp.es == Electrostatics::Ewald) pe = i128_add(pe, ewaldSelf(p, ewald));
        double ke = i128_to_double(i128_from_u128(sumKineticRaw(p))) / 281474976710656.0;
        double E = ke + toDouble(pe);
        result.energy_final = E;
        result.temperature = (2.0 / 3.0) * meanKinetic(p).toDouble() / units::KB;
        result.drift = (E - e0) / (e0 == 0.0 ? 1.0 : e0);
    }

    result.msd = analysis.meanSquaredDisplacement(p, box).toDouble();

    // --- Reproducibility hash + checkpoint --------------------------------
    std::vector<std::uint8_t> bytes = serializeState(p, box.raw(), cfg.steps,
                                                     cfg.seed, rng.state(),
                                                     thermo.xi.raw());
    Sha256 hasher;
    hasher.update(bytes.data(), bytes.size());
    std::uint8_t digest[32];
    hasher.finalize(digest);
    static const char* hexdig = "0123456789abcdef";
    result.sha256.clear();
    for (int i = 0; i < 32; ++i) {
        result.sha256.push_back(hexdig[digest[i] >> 4]);
        result.sha256.push_back(hexdig[digest[i] & 0xF]);
    }

    if (!cfg.output.empty()) {
        writeCheckpoint(cfg.output, p, box.raw(), cfg.steps, cfg.seed, rng.state(), thermo.xi.raw());
    }

    return result;
}

} // namespace mdsim
