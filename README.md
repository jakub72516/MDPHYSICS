# mdsim — Deterministic Molecular Dynamics

A high-performance, **bit-identical** molecular-dynamics simulator in C++14.
The physics core uses integer fixed-point arithmetic so that a simulation
reproduces the exact same trajectory — verified by a SHA-256 fingerprint of the
final state — on any hardware, OS, compiler, or thread count.

## Features

- Deterministic **Q16.48 fixed-point** math core (`Fix64`, portable 128-bit `U128`).
- Fixed-point transcendentals (exp, log, sqrt, erfc, sin/cos via CORDIC).
- Lennard-Jones + Coulomb forces; **Ewald summation** for long-range electrostatics.
- Deterministic **cell list + radix sort** neighbor search.
- Velocity-Verlet integrator, Nosé-Hoover thermostat, isotropic barostat.
- **SHAKE/RATTLE** rigid-bond constraints.
- Static multi-threaded **domain decomposition** (bit-identical for 1–8 threads).
- Versioned, little-endian **checkpoint/restart** + SHA-256 reproducibility hash.
- On-the-fly energy drift, mean-square displacement, radial distribution.

## Build

Requires CMake ≥ 3.10 and a C++14 compiler.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

> On Windows with a MinGW toolchain you may need to pass the compiler and make
> program explicitly:
> `cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_CXX_COMPILER=g++ -DCMAKE_MAKE_PROGRAM=mingw32-make`

## Usage

```
mdsim -N <atoms> -steps <n> -dt <fs> -seed <int> -output <file> [options]

  -N <atoms>       number of atoms            (default 108)
  -steps <n>       number of timesteps        (default 1000)
  -dt <fs>         timestep in fs             (default 1.0)
  -seed <int>      random seed                (default 42)
  -output <file>   checkpoint/output file
  -traj <file>     XYZ trajectory output
  -system <name>   argon | nacl               (default argon)
  -density <rho>   density in atoms/A^3       (default 0.02)
  -temp <K>        temperature in K           (default 300)
  -cutoff <A>      pair cutoff in A           (default 10)
  -es <mode>       electrostatics: 0 none, 1 cutoff, 2 ewald
  -thermostat      enable Nose-Hoover thermostat
  -barostat        enable isotropic barostat
  -threads <n>     worker threads             (default: hardware)
```

Example (NVE argon, 10 ps):

```sh
mdsim -N 108 -steps 1000 -dt 1.0 -seed 42 -density 0.02 -cutoff 8.0 -traj run.xyz
```

Example (NaCl with Ewald):

```sh
mdsim -N 216 -steps 500 -dt 0.5 -seed 7 -system nacl -es 2 -cutoff 8.0 -thermostat
```

## Reproducing the reference result

The test suite includes a 10 ps simulation whose final-state SHA-256 must match
the value in `tests/test_repro.cpp` and `reference/ref_hash.txt`:

```
mdsim -N 108 -steps 1000 -dt 1.0 -seed 42 -density 0.02 -cutoff 8.0
```

Expected SHA-256: see `reference/ref_hash.txt`. Run it on any machine, with any
`-threads` value, and the hash will be identical.

## Analysis

`tools/postprocess.py` reads an XYZ trajectory and reports the radial
distribution function and mean-square displacement:

```sh
python tools/postprocess.py run.xyz --rdf --msd --box 17.5441
```

## Documentation

- `docs/design.md` — deterministic design choices, unit system, data layout.

## Project layout

```
include/mdsim/   headers (fixed-point core, physics, IO, analysis)
src/             implementation + CLI
tests/           unit + reproducibility tests (CTest)
tools/           Python post-processing
docs/            design documentation
reference/       reference SHA-256 hash
```
