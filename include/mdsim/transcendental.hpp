// Deterministic fixed-point transcendental functions.
//
// exp, log, sqrt, erfc, sin/cos are implemented with pure integer arithmetic
// (Taylor / Horner / CORDIC) and a fixed iteration count, so results are
// bit-identical on every platform. No libm functions are used at run time:
// every constant is hard-coded as its exact Q16.48 integer representation,
// generated once with arbitrary-precision arithmetic (see docs/design.md).
#pragma once

#include "mdsim/fixed.hpp"

namespace mdsim {

// e^x  (x must keep the result within Fix64 range; fine for x <= ~10).
Fix64 fexp(Fix64 x);

// ln(x) for x > 0.
Fix64 flog(Fix64 x);

// erfc(x) for all real x.
Fix64 ferfc(Fix64 x);

// sin(x) and cos(x) for all real x (CORDIC, fixed 48 iterations).
void fsincos(Fix64 x, Fix64& sin_out, Fix64& cos_out);

inline Fix64 fsin(Fix64 x) { Fix64 s, c; fsincos(x, s, c); return s; }
inline Fix64 fcos(Fix64 x) { Fix64 s, c; fsincos(x, s, c); return c; }

// Round a Q16.48 value to the nearest integer (ties away from zero).
inline std::int64_t fround_int(std::int64_t v) {
    if (v >= 0) return (v + (Fix64::ONE >> 1)) >> Fix64::FRAC;
    return -(((-v) + (Fix64::ONE >> 1)) >> Fix64::FRAC);
}

} // namespace mdsim
