#include "mdsim/checkpoint.hpp"

#include <cstdio>
#include <cstring>

namespace mdsim {

void put_i64(std::vector<std::uint8_t>& buf, std::int64_t v) {
    for (int i = 0; i < 8; ++i) buf.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
}
void put_u64(std::vector<std::uint8_t>& buf, std::uint64_t v) {
    for (int i = 0; i < 8; ++i) buf.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
}
void put_i32(std::vector<std::uint8_t>& buf, std::int32_t v) {
    for (int i = 0; i < 4; ++i) buf.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
}
void put_u32(std::vector<std::uint8_t>& buf, std::uint32_t v) {
    for (int i = 0; i < 4; ++i) buf.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
}

std::int64_t get_i64(const std::uint8_t* p) {
    std::uint64_t v = 0;
    for (int i = 0; i < 8; ++i) v |= static_cast<std::uint64_t>(p[i]) << (8 * i);
    return static_cast<std::int64_t>(v);
}
std::uint64_t get_u64(const std::uint8_t* p) {
    std::uint64_t v = 0;
    for (int i = 0; i < 8; ++i) v |= static_cast<std::uint64_t>(p[i]) << (8 * i);
    return v;
}
std::int32_t get_i32(const std::uint8_t* p) {
    std::uint32_t v = 0;
    for (int i = 0; i < 4; ++i) v |= static_cast<std::uint32_t>(p[i]) << (8 * i);
    return static_cast<std::int32_t>(v);
}
std::uint32_t get_u32(const std::uint8_t* p) {
    std::uint32_t v = 0;
    for (int i = 0; i < 4; ++i) v |= static_cast<std::uint32_t>(p[i]) << (8 * i);
    return v;
}

std::vector<std::uint8_t> serializeState(const ParticleData& p,
                                         std::int64_t box_raw,
                                         std::int64_t step,
                                         std::uint64_t seed,
                                         std::uint64_t rng_state,
                                         std::int64_t xi_raw) {
    std::vector<std::uint8_t> buf;
    buf.reserve(static_cast<std::size_t>(p.n) * 64 + 128);

    put_i64(buf, box_raw);
    put_i64(buf, step);
    put_u64(buf, seed);
    put_u64(buf, rng_state);
    put_i64(buf, xi_raw);
    put_i32(buf, p.n);

    for (int i = 0; i < p.n; ++i) {
        put_i64(buf, p.x[i]); put_i64(buf, p.y[i]); put_i64(buf, p.z[i]);
        put_i64(buf, p.vx[i]); put_i64(buf, p.vy[i]); put_i64(buf, p.vz[i]);
        put_i32(buf, p.type[i]);
        put_i64(buf, p.q[i]); put_i64(buf, p.m[i]);
    }
    return buf;
}

namespace {
const std::uint32_t CHK_MAGIC = 0x4D445349u; // "MDSI"
const std::uint32_t CHK_VERSION = 1u;
}

bool writeCheckpoint(const std::string& path, const ParticleData& p,
                     std::int64_t box_raw, std::int64_t step,
                     std::uint64_t seed, std::uint64_t rng_state,
                     std::int64_t xi_raw) {
    std::vector<std::uint8_t> body = serializeState(p, box_raw, step, seed, rng_state, xi_raw);

    std::vector<std::uint8_t> out;
    put_u32(out, CHK_MAGIC);
    put_u32(out, CHK_VERSION);
    put_u32(out, static_cast<std::uint32_t>(body.size()));
    out.insert(out.end(), body.begin(), body.end());

    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    std::size_t w = std::fwrite(out.data(), 1, out.size(), f);
    std::fclose(f);
    return w == out.size();
}

bool readCheckpoint(const std::string& path, ParticleData& p,
                    std::int64_t& box_raw, std::int64_t& step,
                    std::uint64_t& seed, std::uint64_t& rng_state,
                    std::int64_t& xi_raw) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return false;
    std::fseek(f, 0, SEEK_END);
    long sz = std::ftell(f);
    std::fseek(f, 0, SEEK_SET);
    std::vector<std::uint8_t> out(static_cast<std::size_t>(sz));
    if (sz > 0) std::fread(out.data(), 1, static_cast<std::size_t>(sz), f);
    std::fclose(f);

    if (out.size() < 12) return false;
    if (get_u32(out.data()) != CHK_MAGIC) return false;
    if (get_u32(out.data() + 4) != CHK_VERSION) return false;
    std::uint32_t body_size = get_u32(out.data() + 8);
    const std::uint8_t* b = out.data() + 12;
    if (out.size() < 12 + body_size) return false;

    box_raw = get_i64(b + 0);
    step = get_i64(b + 8);
    seed = get_u64(b + 16);
    rng_state = get_u64(b + 24);
    xi_raw = get_i64(b + 32);
    int n = get_i32(b + 40);
    p.resize(n);

    const std::uint8_t* q = b + 44;
    for (int i = 0; i < n; ++i) {
        p.x[i] = get_i64(q); q += 8;
        p.y[i] = get_i64(q); q += 8;
        p.z[i] = get_i64(q); q += 8;
        p.vx[i] = get_i64(q); q += 8;
        p.vy[i] = get_i64(q); q += 8;
        p.vz[i] = get_i64(q); q += 8;
        p.type[i] = get_i32(q); q += 4;
        p.q[i] = get_i64(q); q += 8;
        p.m[i] = get_i64(q); q += 8;
    }
    return true;
}

bool appendXyzFrame(const std::string& path, const ParticleData& p,
                    const std::vector<std::string>& symbols,
                    Fix64 box, std::int64_t step) {
    FILE* f = std::fopen(path.c_str(), "a");
    if (!f) return false;
    std::fprintf(f, "%d\n", p.n);
    std::fprintf(f, "step %lld box %.6f\n", static_cast<long long>(step), box.toDouble());
    for (int i = 0; i < p.n; ++i) {
        const std::string& s = symbols[p.type[i]];
        std::fprintf(f, "%s %.6f %.6f %.6f\n", s.c_str(),
                     Fix64::fromRaw(p.x[i]).toDouble(),
                     Fix64::fromRaw(p.y[i]).toDouble(),
                     Fix64::fromRaw(p.z[i]).toDouble());
    }
    std::fclose(f);
    return true;
}

} // namespace mdsim
