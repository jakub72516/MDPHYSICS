// Portable 128-bit unsigned arithmetic.
//
// The project targets 32-bit and 64-bit platforms (including MSVC and
// MinGW 32-bit, neither of which provides __int128). All operations are
// implemented with explicit 64-bit and 32-bit integer arithmetic so that
// results are bit-identical on every IEEE-754-capable C++ toolchain.
#pragma once

#include <cstdint>

namespace mdsim {

struct U128 {
    std::uint64_t lo;
    std::uint64_t hi;
};

// 64x64 -> 128 unsigned multiply (Hacker's Delight, muldu).
inline U128 umul64(std::uint64_t a, std::uint64_t b) {
    std::uint64_t a0 = a & 0xFFFFFFFFu;
    std::uint64_t a1 = a >> 32;
    std::uint64_t b0 = b & 0xFFFFFFFFu;
    std::uint64_t b1 = b >> 32;

    std::uint64_t w0 = a0 * b0;
    std::uint64_t t  = a1 * b0 + (w0 >> 32);
    std::uint64_t w1 = t & 0xFFFFFFFFu;
    std::uint64_t w2 = t >> 32;
    w1 = a0 * b1 + w1;

    U128 r;
    r.hi = a1 * b1 + w2 + (w1 >> 32);
    r.lo = (w1 << 32) + (w0 & 0xFFFFFFFFu);
    return r;
}

inline U128 uadd(U128 a, U128 b) {
    U128 r;
    r.lo = a.lo + b.lo;
    r.hi = a.hi + b.hi + (r.lo < a.lo ? 1u : 0u);
    return r;
}

inline U128 usub(U128 a, U128 b) {
    U128 r;
    r.lo = a.lo - b.lo;
    r.hi = a.hi - b.hi - (a.lo < b.lo ? 1u : 0u);
    return r;
}

inline int ucmp(U128 a, U128 b) {
    if (a.hi != b.hi) return a.hi < b.hi ? -1 : 1;
    if (a.lo != b.lo) return a.lo < b.lo ? -1 : 1;
    return 0;
}

inline U128 ushl(U128 a, int s) {
    if (s <= 0) return a;
    if (s >= 128) return U128{0, 0};
    if (s >= 64) return U128{0, a.lo << (s - 64)};
    U128 r;
    r.hi = (a.hi << s) | (a.lo >> (64 - s));
    r.lo = a.lo << s;
    return r;
}

inline U128 ushr(U128 a, int s) {
    if (s <= 0) return a;
    if (s >= 128) return U128{0, 0};
    if (s >= 64) return U128{a.hi >> (s - 64), 0};
    U128 r;
    r.lo = (a.lo >> s) | (a.hi << (64 - s));
    r.hi = a.hi >> s;
    return r;
}

// 128 / 64 division returning the 64-bit quotient (floor).
// Requires d != 0. Correct for quotients that fit in 64 bits (which is the
// only case exercised by the fixed-point core).
inline std::uint64_t udiv128_64(U128 n, std::uint64_t d) {
    U128 r{0, 0};
    std::uint64_t q = 0;
    for (int i = 127; i >= 0; --i) {
        std::uint64_t bit = (i >= 64) ? ((n.hi >> (i - 64)) & 1u) : ((n.lo >> i) & 1u);
        std::uint64_t hi = (r.hi << 1) | (r.lo >> 63);
        std::uint64_t lo = (r.lo << 1) | bit;
        r = U128{lo, hi};
        q <<= 1;
        if (r.hi != 0 || r.lo >= d) {
            std::uint64_t lo2 = r.lo - d;
            std::uint64_t hi2 = r.hi - (lo2 > r.lo ? 1u : 0u);
            r = U128{lo2, hi2};
            q |= 1;
        }
    }
    return q;
}

// 128-bit integer square root (floor), returning a 64-bit result.
inline std::uint64_t isqrt128(U128 n) {
    std::uint64_t x = 0;
    for (int b = 63; b >= 0; --b) {
        std::uint64_t trial = x | (std::uint64_t(1) << b);
        if (ucmp(umul64(trial, trial), n) <= 0) x = trial;
    }
    return x;
}

// Signed 128-bit accumulator used for extensive sums (e.g. total energy).
struct I128 {
    std::uint64_t lo;
    std::int64_t hi;
};

inline I128 i128_zero() { return I128{0, 0}; }

inline I128 i128_from_i64(std::int64_t x) {
    return I128{static_cast<std::uint64_t>(x), x < 0 ? -1 : 0};
}

// Add a signed 64-bit value to a signed 128-bit accumulator.
inline I128 i128_add_i64(I128 a, std::int64_t b) {
    std::uint64_t lo = a.lo + static_cast<std::uint64_t>(b);
    bool carry = lo < a.lo;
    std::int64_t hi = a.hi + (b < 0 ? -1 : 0) + (carry ? 1 : 0);
    return I128{lo, hi};
}

// Convert to double (display/analysis only; the accumulator is exact).
inline double i128_to_double(I128 a) {
    return static_cast<double>(a.hi) * 18446744073709551616.0 +
           static_cast<double>(a.lo);
}

inline I128 i128_from_u128(U128 u) {
    return I128{u.lo, static_cast<std::int64_t>(u.hi)};
}

inline I128 i128_neg(I128 a) {
    std::uint64_t lo = ~a.lo;
    std::int64_t hi = ~a.hi;
    lo += 1;
    if (lo == 0) hi += 1;
    return I128{lo, hi};
}

inline I128 i128_add(I128 a, I128 b) {
    std::uint64_t lo = a.lo + b.lo;
    bool carry = lo < a.lo;
    std::int64_t hi = a.hi + b.hi + (carry ? 1 : 0);
    return I128{lo, hi};
}

inline I128 i128_shl1(I128 a) {
    return I128{a.lo << 1, (a.hi << 1) | static_cast<std::int64_t>(a.lo >> 63)};
}

// Signed 128/64 division (truncation). Quotient must fit in 64 bits.
inline std::int64_t i128_div_i64(I128 n, std::int64_t d) {
    bool neg = (n.hi < 0) != (d < 0);
    I128 mag = n.hi < 0 ? i128_neg(n) : n;
    std::uint64_t ud = d < 0 ? (static_cast<std::uint64_t>(-(d + 1)) + 1)
                             : static_cast<std::uint64_t>(d);
    std::uint64_t q = udiv128_64(U128{mag.lo, static_cast<std::uint64_t>(mag.hi)}, ud);
    std::int64_t r = static_cast<std::int64_t>(q);
    return neg ? -r : r;
}

} // namespace mdsim
