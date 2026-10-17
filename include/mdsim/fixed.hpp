// Fix64: deterministic signed 64-bit fixed-point arithmetic (Q16.48).
//
// A Fix64 value `v` represents the real number v / 2^48. This gives a
// representable range of [-32768, 32768) and a resolution of 2^-48 (~3.55e-15),
// comparable to IEEE-754 double precision near unity.
//
// Every arithmetic operation is implemented with explicit integer math
// (via U128) so that the results are bit-identical on every platform,
// regardless of compiler, FPU precision, SIMD width, or OS. This is the
// foundation of the project's determinism guarantee.
#pragma once

#include <cstdint>
#include <string>

#include "mdsim/u128.hpp"

namespace mdsim {

// Forward declarations of deterministic free functions used by Fix64.
inline std::int64_t fadd(std::int64_t a, std::int64_t b);
inline std::int64_t fsub(std::int64_t a, std::int64_t b);
inline std::int64_t fmul(std::int64_t a, std::int64_t b);
inline std::int64_t fdiv(std::int64_t a, std::int64_t b);
inline std::int64_t finv(std::int64_t a);
inline std::int64_t fsqrt(std::int64_t a);

class Fix64 {
public:
    static const int FRAC = 48;
    static const std::int64_t ONE = std::int64_t(1) << FRAC;
    static const std::int64_t MAX_V = std::int64_t(0x7FFFFFFFFFFFFFFFLL);
    static const std::int64_t MIN_V = -MAX_V - 1;

    std::int64_t v;

    Fix64() : v(0) {}
    explicit Fix64(std::int64_t raw) : v(raw) {}

    static Fix64 fromRaw(std::int64_t raw) { return Fix64(raw); }
    std::int64_t raw() const { return v; }

    static Fix64 fromDouble(double d);
    double toDouble() const;
    std::string toString(int decimals = 6) const;

    // ---- comparison ------------------------------------------------------
    bool operator==(const Fix64& o) const { return v == o.v; }
    bool operator!=(const Fix64& o) const { return v != o.v; }
    bool operator<(const Fix64& o) const { return v < o.v; }
    bool operator>(const Fix64& o) const { return v > o.v; }
    bool operator<=(const Fix64& o) const { return v <= o.v; }
    bool operator>=(const Fix64& o) const { return v >= o.v; }

    // ---- unary -----------------------------------------------------------
    Fix64 operator-() const { return Fix64(-v); }
    Fix64 operator+() const { return *this; }

    // ---- arithmetic (delegates to deterministic free functions) ---------
    Fix64 operator+(const Fix64& o) const { return Fix64(fadd(v, o.v)); }
    Fix64 operator-(const Fix64& o) const { return Fix64(fsub(v, o.v)); }
    Fix64 operator*(const Fix64& o) const { return Fix64(fmul(v, o.v)); }
    Fix64 operator/(const Fix64& o) const { return Fix64(fdiv(v, o.v)); }

    Fix64& operator+=(const Fix64& o) { v = fadd(v, o.v); return *this; }
    Fix64& operator-=(const Fix64& o) { v = fsub(v, o.v); return *this; }
    Fix64& operator*=(const Fix64& o) { v = fmul(v, o.v); return *this; }
    Fix64& operator/=(const Fix64& o) { v = fdiv(v, o.v); return *this; }

    Fix64 abs() const { return Fix64(v < 0 ? -v : v); }

    Fix64 inv() const { return Fix64(finv(v)); }
    Fix64 sqrt() const { return Fix64(fsqrt(v)); }
};

// ---- deterministic free functions ----------------------------------------

inline std::int64_t fadd(std::int64_t a, std::int64_t b) { return a + b; }
inline std::int64_t fsub(std::int64_t a, std::int64_t b) { return a - b; }

// a * b, rounded to nearest (ties away from zero).
inline std::int64_t fmul(std::int64_t a, std::int64_t b) {
    bool neg = (a < 0) != (b < 0);
    std::uint64_t ua = a < 0 ? std::uint64_t(0) - std::uint64_t(a) : std::uint64_t(a);
    std::uint64_t ub = b < 0 ? std::uint64_t(0) - std::uint64_t(b) : std::uint64_t(b);
    U128 p = umul64(ua, ub);
    p = uadd(p, U128{std::uint64_t(1) << (Fix64::FRAC - 1), 0});
    p = ushr(p, Fix64::FRAC);
    std::int64_t mag = static_cast<std::int64_t>(p.lo);
    return neg ? (std::int64_t(0) - mag) : mag;
}

// a / b, rounded to nearest (ties away from zero).
inline std::int64_t fdiv(std::int64_t a, std::int64_t b) {
    bool neg = (a < 0) != (b < 0);
    std::uint64_t ua = a < 0 ? std::uint64_t(0) - std::uint64_t(a) : std::uint64_t(a);
    std::uint64_t ub = b < 0 ? std::uint64_t(0) - std::uint64_t(b) : std::uint64_t(b);
    U128 n = ushl(U128{ua, 0}, Fix64::FRAC);
    n = uadd(n, U128{ub >> 1, 0});
    std::uint64_t q = udiv128_64(n, ub);
    std::int64_t mag = static_cast<std::int64_t>(q);
    return neg ? (std::int64_t(0) - mag) : mag;
}

// 1 / a, rounded to nearest (ties away from zero).
inline std::int64_t finv(std::int64_t a) {
    bool neg = a < 0;
    std::uint64_t ua = neg ? std::uint64_t(0) - std::uint64_t(a) : std::uint64_t(a);
    U128 n = ushl(U128{1, 0}, 2 * Fix64::FRAC);
    n = uadd(n, U128{ua >> 1, 0});
    std::uint64_t q = udiv128_64(n, ua);
    std::int64_t mag = static_cast<std::int64_t>(q);
    return neg ? (std::int64_t(0) - mag) : mag;
}

// sqrt(a) for a >= 0.
inline std::int64_t fsqrt(std::int64_t a) {
    if (a <= 0) return 0;
    std::uint64_t ua = static_cast<std::uint64_t>(a);
    U128 n = ushl(U128{ua, 0}, Fix64::FRAC);
    return static_cast<std::int64_t>(isqrt128(n));
}

} // namespace mdsim
