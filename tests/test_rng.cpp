// Tests for deterministic RNG and transcendental functions.
#include <cstdio>
#include <cmath>

#include "mdsim/rng.hpp"
#include "mdsim/transcendental.hpp"

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
    // RNG determinism: same seed -> same sequence.
    Rng r1(12345);
    Rng r2(12345);
    for (int i = 0; i < 1000; ++i) CHECK(r1.next() == r2.next());

    // Transcendentals vs libm.
    neard("exp(1)", fexp(Fix64::fromDouble(1.0)).toDouble(), std::exp(1.0), 1e-12);
    neard("exp(-2)", fexp(Fix64::fromDouble(-2.0)).toDouble(), std::exp(-2.0), 1e-12);
    neard("exp(3.5)", fexp(Fix64::fromDouble(3.5)).toDouble(), std::exp(3.5), 1e-11);
    neard("log(2)", flog(Fix64::fromDouble(2.0)).toDouble(), std::log(2.0), 1e-12);
    neard("log(0.5)", flog(Fix64::fromDouble(0.5)).toDouble(), std::log(0.5), 1e-12);
    neard("log(7.3)", flog(Fix64::fromDouble(7.3)).toDouble(), std::log(7.3), 1e-12);

    neard("erfc(0)", ferfc(Fix64::fromDouble(0.0)).toDouble(), 1.0, 1e-6);
    neard("erfc(1)", ferfc(Fix64::fromDouble(1.0)).toDouble(), std::erfc(1.0), 1e-5);
    neard("erfc(2)", ferfc(Fix64::fromDouble(2.0)).toDouble(), std::erfc(2.0), 1e-5);
    neard("erfc(-1)", ferfc(Fix64::fromDouble(-1.0)).toDouble(), std::erfc(-1.0), 1e-5);

    // sin/cos over a range of angles.
    for (double ang = -6.0; ang <= 6.0; ang += 0.13) {
        Fix64 s, c;
        fsincos(Fix64::fromDouble(ang), s, c);
        neard("sin", s.toDouble(), std::sin(ang), 1e-9);
        neard("cos", c.toDouble(), std::cos(ang), 1e-9);
    }
    neard("sin(31)", fsin(Fix64::fromDouble(31.0)).toDouble(), std::sin(31.0), 1e-8);
    neard("cos(-25)", fcos(Fix64::fromDouble(-25.0)).toDouble(), std::cos(-25.0), 1e-8);

    if (fails == 0) std::printf("test_rng: ALL PASS\n");
    else std::printf("test_rng: %d FAILURES\n", fails);
    return fails == 0 ? 0 : 1;
}
