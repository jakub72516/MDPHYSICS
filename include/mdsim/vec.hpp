// Minimal 3-vector over Fix64, used for setup and analysis (the hot path
// uses structure-of-arrays storage, see particles.hpp).
#pragma once

#include "mdsim/fixed.hpp"

namespace mdsim {

struct Vec3 {
    Fix64 x, y, z;

    Vec3() : x(0), y(0), z(0) {}
    Vec3(Fix64 a, Fix64 b, Fix64 c) : x(a), y(b), z(c) {}

    Vec3 operator+(const Vec3& o) const { return Vec3(x + o.x, y + o.y, z + o.z); }
    Vec3 operator-(const Vec3& o) const { return Vec3(x - o.x, y - o.y, z - o.z); }
    Vec3 operator*(Fix64 s) const { return Vec3(x * s, y * s, z * s); }

    Fix64 dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }

    Fix64 normSq() const { return dot(*this); }
};

} // namespace mdsim
