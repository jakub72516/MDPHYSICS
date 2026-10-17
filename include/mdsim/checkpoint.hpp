// Checkpoint/restart with versioned, little-endian binary format.
//
// State is serialized to a canonical little-endian byte stream so that a
// checkpoint written on any machine can be read back elsewhere and reproduce
// the exact same trajectory. The same serialization feeds the SHA-256
// reproducibility hash.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "mdsim/fixed.hpp"
#include "mdsim/particles.hpp"

namespace mdsim {

// --- little-endian byte helpers ---
void put_i64(std::vector<std::uint8_t>& buf, std::int64_t v);
void put_u64(std::vector<std::uint8_t>& buf, std::uint64_t v);
void put_i32(std::vector<std::uint8_t>& buf, std::int32_t v);
void put_u32(std::vector<std::uint8_t>& buf, std::uint32_t v);
std::int64_t get_i64(const std::uint8_t* p);
std::uint64_t get_u64(const std::uint8_t* p);
std::int32_t get_i32(const std::uint8_t* p);
std::uint32_t get_u32(const std::uint8_t* p);

// Canonical deterministic serialization of the reproducible simulation state.
std::vector<std::uint8_t> serializeState(const ParticleData& p,
                                         std::int64_t box_raw,
                                         std::int64_t step,
                                         std::uint64_t seed,
                                         std::uint64_t rng_state,
                                         std::int64_t xi_raw);

bool writeCheckpoint(const std::string& path, const ParticleData& p,
                     std::int64_t box_raw, std::int64_t step,
                     std::uint64_t seed, std::uint64_t rng_state,
                     std::int64_t xi_raw);

bool readCheckpoint(const std::string& path, ParticleData& p,
                    std::int64_t& box_raw, std::int64_t& step,
                    std::uint64_t& seed, std::uint64_t& rng_state,
                    std::int64_t& xi_raw);

// Append one XYZ trajectory frame (used for post-processing).
bool appendXyzFrame(const std::string& path, const ParticleData& p,
                    const std::vector<std::string>& symbols,
                    Fix64 box, std::int64_t step);

} // namespace mdsim
