// Deterministic pseudo-random number generation.
//
// xorshift64* is used for the stream and splitmix64 for seed expansion. Both
// are pure integer operations, giving identical sequences on every platform.
#pragma once

#include <cstdint>

#include "mdsim/fixed.hpp"

namespace mdsim {

class Rng {
public:
    explicit Rng(std::uint64_t seed) : s_(seed) {
        if (s_ == 0) s_ = 0x9E3779B97F4A7C15ull;
    }

    std::uint64_t next() {
        std::uint64_t x = s_;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        s_ = x;
        return x * 0x2545F4914F6CDD1Dull;
    }

    std::uint64_t state() const { return s_; }
    void setState(std::uint64_t s) { s_ = s; }

    // Uniform in [0, 1) as a Q16.48 Fix64 (top 48 bits of next()).
    Fix64 uniform() {
        std::uint64_t x = next();
        return Fix64::fromRaw(static_cast<std::int64_t>(x >> 16));
    }

private:
    std::uint64_t s_;
};

// splitmix64: deterministic seed mixer / expander.
inline std::uint64_t splitmix64(std::uint64_t& st) {
    std::uint64_t z = (st += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

// Standard-normal sample (Box-Muller), fully in fixed point.
Fix64 gaussian(Rng& rng);

} // namespace mdsim
