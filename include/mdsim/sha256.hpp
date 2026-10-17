// Self-contained SHA-256 (FIPS 180-4), used for reproducibility verification.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace mdsim {

struct Sha256 {
    std::uint32_t state[8];
    std::uint64_t bitlen = 0;
    std::uint8_t buffer[64];
    std::size_t buflen = 0;

    Sha256();
    void update(const std::uint8_t* data, std::size_t len);
    void update(const char* data, std::size_t len) {
        update(reinterpret_cast<const std::uint8_t*>(data), len);
    }
    void finalize(std::uint8_t digest[32]);

    std::string hex() const; // after finalize
private:
    void process(const std::uint8_t* block);
};

} // namespace mdsim
