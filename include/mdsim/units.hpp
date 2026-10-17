// Unit system and physical constants.
//
// Units (LAMMPS "real"-style):
//   length : angstrom (A)
//   time   : femtosecond (fs)
//   mass   : g/mol (a.k.a. atomic mass unit / Dalton)
//   energy : kcal/mol
//   charge : elementary charge e
//
// Derived conversion used by the integrator (see docs/design.md):
//   acceleration [A/fs^2] = ACCEL_SCALE * force[kcal/(mol*A)] / mass[g/mol]
// where ACCEL_SCALE = 4.184e-4 (1 kcal = 4184 J, 1 A = 1e-10 m, 1 fs = 1e-15 s).
//
// These are tunable parameters stored as decimal double literals. They are
// converted to Fix64 exactly once at start-up via Fix64::fromDouble, which is
// deterministic because decimal-to-binary parsing is correctly rounded on all
// conforming compilers.
#pragma once

namespace mdsim {
namespace units {

// Force -> acceleration conversion (A/fs^2 per (kcal/(mol*A))/(g/mol)).
constexpr double ACCEL_SCALE = 4.184e-4;

// Coulomb prefactor U = COULOMB * q_i q_j / r [kcal/mol] (q in e, r in A).
constexpr double COULOMB = 332.0637;

// Boltzmann constant [kcal/(mol*K)].
constexpr double KB = 0.0019872041;

} // namespace units
} // namespace mdsim
