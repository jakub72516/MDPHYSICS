#include "mdsim/fixed.hpp"

#include <cstdio>
#include <cmath>

namespace mdsim {

Fix64 Fix64::fromDouble(double d) {
    // Clamp to the representable range.
    const double scale = 281474976710656.0; // 2^48
    if (d >= 32768.0) d = 32767.999999999;
    if (d <= -32768.0) d = -32768.0;
    double scaled = d * scale;
    std::int64_t r;
    if (scaled >= 0.0)
        r = static_cast<std::int64_t>(scaled + 0.5);
    else
        r = static_cast<std::int64_t>(scaled - 0.5);
    return Fix64(r);
}

double Fix64::toDouble() const {
    return static_cast<double>(v) / 281474976710656.0;
}

std::string Fix64::toString(int decimals) const {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*f", decimals, toDouble());
    return std::string(buf);
}

} // namespace mdsim
