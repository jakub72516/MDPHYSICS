// Unit tests for the deterministic fixed-point math core.
#include <cstdio>
#include <cmath>

#include "mdsim/fixed.hpp"
#include "mdsim/u128.hpp"

using namespace mdsim;

static int fails = 0;
#define CHECK(cond) do { if (!(cond)) { std::printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond); ++fails; } } while (0)

static void neard(const char* name, double got, double want, double tol) {
    double err = std::fabs(got - want);
    if (err > tol) {
        std::printf("FAIL %-16s got=%.12g want=%.12g err=%.3g\n", name, got, want, err);
        ++fails;
    }
}

int main() {
    // 128-bit multiply / divide against known vectors.
    U128 p = umul64(0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull);
    CHECK(p.lo == 1 && p.hi == 0xFFFFFFFFFFFFFFFEull);

    U128 q = umul64(123456789012345ull, 987654321098765ull);
    CHECK(q.lo == 14417890538969770277ull && q.hi == 6609981178ull);

    U128 num{12345678901234567890ull, 6ull};
    CHECK(udiv128_64(num, 987654321ull) == 124563970133ull);

    CHECK(isqrt128(umul64(1000000ull, 1000000ull)) == 1000000ull);

    // Basic fixed-point arithmetic.
    Fix64 a = Fix64::fromDouble(3.0);
    Fix64 b = Fix64::fromDouble(4.0);
    neard("add", (a + b).toDouble(), 7.0, 1e-12);
    neard("sub", (a - b).toDouble(), -1.0, 1e-12);
    neard("mul", (a * b).toDouble(), 12.0, 1e-12);
    neard("div", (a / b).toDouble(), 0.75, 1e-12);
    neard("inv", b.inv().toDouble(), 0.25, 1e-12);
    neard("sqrt", Fix64::fromDouble(2.0).sqrt().toDouble(), std::sqrt(2.0), 1e-12);

    neard("div big", (Fix64::fromDouble(1000.0) / Fix64::fromDouble(0.1)).toDouble(), 10000.0, 1e-6);
    neard("mul big", (Fix64::fromDouble(12345.67) * Fix64::fromDouble(0.0001)).toDouble(), 1.234567, 1e-9);
    neard("negmul", (Fix64::fromDouble(-3.0) * Fix64::fromDouble(4.0)).toDouble(), -12.0, 1e-12);
    neard("negdiv", (Fix64::fromDouble(-3.0) / Fix64::fromDouble(4.0)).toDouble(), -0.75, 1e-12);
    neard("neginv", Fix64::fromDouble(-4.0).inv().toDouble(), -0.25, 1e-12);

    // Determinism of fromDouble.
    double d = 0.123456789012345;
    CHECK(Fix64::fromDouble(d).raw() == Fix64::fromDouble(d).raw());

    if (fails == 0) std::printf("test_fixed: ALL PASS\n");
    else std::printf("test_fixed: %d FAILURES\n", fails);
    return fails == 0 ? 0 : 1;
}
