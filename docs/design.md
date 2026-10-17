# Design: Deterministic Molecular Dynamics

`mdsim` is a deterministic, high-performance molecular-dynamics simulator.
The central guarantee is **bit-identical trajectories** on any IEEE-754
platform, regardless of CPU, SIMD width, OS, compiler, or thread count. This
document describes the deterministic design choices and how they are achieved.

## 1. Deterministic fixed-point arithmetic

The single most important decision is to perform **all** physics arithmetic in
integer fixed-point rather than floating point. Floating-point addition is not
associative, so even reordering terms (or letting the compiler contract a
`mul`/`add` into an FMA) changes the result. Integer arithmetic is exact and
order-independent, which is exactly what determinism requires.

- `Fix64` (`include/mdsim/fixed.hpp`) is a signed 64-bit fixed-point type in
  **Q16.48** format: the value is `raw / 2^48`, giving a range of
  `[-32768, 32768)` and a resolution of `2^-48` (~3.55e-15), comparable to a
  double near unity.

- All multiplications use a portable 128-bit intermediate (`U128`,
  `include/mdsim/u128.hpp`) implemented with 32-bit limb arithmetic, because
  the reference toolchain (MinGW 32-bit GCC 6.3) has no `__int128` and MSVC
  likewise lacks it. Division, reciprocal, and square root are implemented with
  explicit 128-bit long division / Newton's method so results are identical on
  every compiler.

- Rounding is round-to-nearest, ties away from zero, applied identically in
  `fmul`, `fdiv`, and `finv`.

### Transcendentals

`exp`, `log`, `sqrt`, `erfc`, and `sin`/`cos` (see
`src/transcendental.cpp`) are computed with fixed iteration counts from pure
integer arithmetic:

- `exp` / `log`: range reduction + Taylor series.
- `erfc`: Abramowitz–Stegun 7.1.26 polynomial.
- `sin`/`cos`: **CORDIC** (48 fixed rotations), which avoids the overflow that
  a naive power series hits at fixed point and is naturally deterministic.

Every constant is hard-coded as its exact Q16.48 integer value, generated once
with arbitrary-precision arithmetic (mpmath), so no libm call influences the
trajectory. The only libm use is a correctly-rounded `sqrt` during initial
velocity rescaling.

## 2. Deterministic random numbers

`Rng` (`include/mdsim/rng.hpp`) is a 64-bit **xorshift64\*** stream seeded with
`splitmix64`. Initial velocities are drawn via a Box–Muller transform
implemented entirely in fixed point. No stochastic term is used in the dynamics
(the Nosé–Hoover thermostat is deterministic); the RNG only seeds initial
conditions, which are then captured by the checkpoint.

## 3. Deterministic neighbor search

- Particles are assigned to cubic cells and sorted by cell key with a stable
  **LSD radix sort** (`src/cell_list.cpp`). Because the sort is stable and the
  input is ordered by particle index, the intra-cell ordering is identical on
  every run.

- The cell list requires `nc >= 3` cells per dimension so that every neighbor
  cell is reached by a unique periodic offset; the simulation enforces this by
  clamping the cutoff to `L/3`.

- The neighbor-pair list (`src/neighbor_list.cpp`) enumerates unique pairs
  `i < j` with a linear-index half-partition over the sorted cells, giving a
  fixed pair ordering.

## 4. Deterministic parallelism

`ThreadPool` (`src/thread_pool.cpp`) is a static scheduler: workers pull
contiguous *blocks* from an atomic counter and results are combined in fixed
order. (It is backed by Win32 threads because the reference toolchain has no
`std::thread`; on POSIX it uses `std::thread`.)

Forces are computed with a **full-neighbor, owner-cell** decomposition
(`src/forces.cpp`): each block owns a contiguous range of cells and computes
the *complete* force on its particles by iterating the 27-cell stencil in a
fixed order. Because a particle's force is written only by its owning block and
its summation order is fixed, the result is independent of thread count (this
is verified by the reproducibility test across 1–8 threads).

## 5. Unit system

| Quantity | Unit |
|----------|------|
| length   | angstrom (Å) |
| time     | femtosecond (fs) |
| mass     | g/mol (Dalton) |
| energy   | kcal/mol |
| charge   | elementary charge e |

The integrator uses `acceleration = ACCEL_SCALE * force / mass` with
`ACCEL_SCALE = 4.184e-4` (1 kcal = 4184 J, 1 Å = 1e-10 m, 1 fs = 1e-15 s).
The Coulomb prefactor is `332.0637 kcal·Å/(mol·e²)`.

## 6. Data layout

Dynamical state is stored **Structure-of-Arrays** (`ParticleData`,
`include/mdsim/particles.hpp`): contiguous `int64` arrays for `x/y/z`,
`vx/vy/vz`, and `fx/fy/fz`. Static attributes (`type`, `q`, `m`, `inv_m`) are
separate arrays. This layout is cache-friendly and SIMD-ready.

## 7. Checkpoint / reproducibility

`src/checkpoint.cpp` serializes the reproducible state (box, step, seed, RNG
state, thermostat friction, and all particle data) into a versioned,
little-endian byte stream. The same stream is hashed with a self-contained
SHA-256 (`src/sha256.cpp`) to produce the reproducibility fingerprint. A
checkpoint written on one machine restores the identical trajectory elsewhere.

## 8. Known limitations / scope

- **Fixed-point range**: `Q16.48` limits coordinates/forces to `±32768` and
  can overflow the repulsive `r^-12` term for separations below ~0.7σ; a
  soft-core exclusion radius prevents this (never reached in stable runs).
- **Extensive quantities** (total energy, virial) overflow 64 bits for very
  large systems; they are accumulated in 128-bit (`U128`/`I128`) and reduced to
  intensive (per-atom) values where needed.
- **Ewald** is a direct reciprocal-space k-sum (O(N·K)); it is deterministic but
  intended for systems up to a few hundred charges. A PME/FFT back-end would
  scale it to millions of atoms.
- The barostat is an isotropic Berendsen-style controller (simpler than a full
  Parrinello–Rahman cell); the thermostat is a single Nosé–Hoover chain link.
- SIMD (AVX2/AVX-512) intrinsics are not yet used; the fixed-point inner loop
  is scalar. This is the primary target for the performance milestone.
