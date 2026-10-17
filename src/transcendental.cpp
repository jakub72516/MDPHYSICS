#include "mdsim/transcendental.hpp"

namespace mdsim {

namespace {

// Exact Q16.48 constants (integer representations), generated with mpmath.
const std::int64_t C_PI     = 884279719003555LL;
const std::int64_t C_TWO_PI = 1768559438007110LL;
const std::int64_t C_HALF_PI = 442139859501778LL;
const std::int64_t C_LN2    = 195103586505167LL;
const std::int64_t C_LOG2E  = 406082553034800LL;

const std::int64_t C_CORDIC_K = 170926505739102LL; // 0.607252935008881256

const int CORDIC_N = 48;
const std::int64_t C_ATAN[48] = {
    221069929750889LL, 130505199945453LL, 68955363498242LL, 35002819193903LL,
    17569333089919LL, 8793231387230LL, 4397688649582LL, 2198978517948LL,
    1099506035422LL, 549755114839LL, 274877819563LL, 137438942549LL,
    68719475371LL, 34359738197LL, 17179869163LL, 8589934589LL,
    4294967296LL, 2147483648LL, 1073741824LL, 536870912LL,
    268435456LL, 134217728LL, 67108864LL, 33554432LL,
    16777216LL, 8388608LL, 4194304LL, 2097152LL,
    1048576LL, 524288LL, 262144LL, 131072LL,
    65536LL, 32768LL, 16384LL, 8192LL,
    4096LL, 2048LL, 1024LL, 512LL,
    256LL, 128LL, 64LL, 32LL,
    16LL, 8LL, 4LL, 2LL
};

const std::int64_t C_ERF_P  = 92208697243118LL;
const std::int64_t C_ERF_A1 = 71728153473386LL;
const std::int64_t C_ERF_A2 = -80078712139858LL;
const std::int64_t C_ERF_A3 = 400092399644181LL;
const std::int64_t C_ERF_A4 = -409025932956868LL;
const std::int64_t C_ERF_A5 = 298759068408339LL;

inline Fix64 fromInt(std::int64_t n) { return Fix64(n * Fix64::ONE); }

} // namespace

Fix64 fexp(Fix64 x) {
    // e^x = 2^(x/ln2). Split x = k*ln2 + r with k = round(x/ln2), |r| <= ln2/2.
    std::int64_t k = fround_int(fmul(x.raw(), C_LOG2E));
    Fix64 r = Fix64(fsub(x.raw(), fmul(fromInt(k).raw(), C_LN2)));

    // e^r via Taylor: sum r^n/n!, 14 terms suffices for |r| <= ln2/2.
    Fix64 sum = Fix64(Fix64::ONE);
    Fix64 term = Fix64(Fix64::ONE);
    for (std::int64_t n = 1; n <= 14; ++n) {
        term = Fix64(fdiv(fmul(term.raw(), r.raw()), fromInt(n).raw()));
        sum = Fix64(fadd(sum.raw(), term.raw()));
    }

    // Multiply by 2^k. sum is always positive here.
    std::int64_t raw = sum.raw();
    if (k >= 0) {
        raw = raw << k;
    } else {
        raw = raw >> (-k);
    }
    return Fix64(raw);
}

Fix64 flog(Fix64 x) {
    if (x.raw() <= 0) return Fix64(0);
    // Normalize x = m * 2^e, m in [1, 2).
    std::int64_t raw = x.raw();
    std::int64_t e = 0;
    while (raw >= (Fix64::ONE << 1)) { raw >>= 1; ++e; }
    while (raw < Fix64::ONE) { raw <<= 1; --e; }
    Fix64 m = Fix64::fromRaw(raw);

    // ln(m) = 2*atanh((m-1)/(m+1)) = 2*(t + t^3/3 + t^5/5 + ...).
    Fix64 t = Fix64(fdiv(fsub(m.raw(), Fix64::ONE), fadd(m.raw(), Fix64::ONE)));
    Fix64 t2 = Fix64(fmul(t.raw(), t.raw()));
    Fix64 sum = t;
    Fix64 tp = t;
    for (std::int64_t n = 1; n <= 15; ++n) {
        tp = Fix64(fmul(tp.raw(), t2.raw()));
        sum = Fix64(fadd(sum.raw(), fdiv(tp.raw(), fromInt(2 * n + 1).raw())));
    }
    Fix64 lnm = Fix64(sum.raw() + sum.raw());
    return Fix64(fadd(lnm.raw(), fmul(fromInt(e).raw(), C_LN2)));
}

Fix64 ferfc(Fix64 x) {
    if (x.raw() < 0) {
        // erfc(x) = 2 - erfc(-x)
        return Fix64(fsub(fromInt(2).raw(), ferfc(Fix64(-x.raw())).raw()));
    }
    // Abramowitz & Stegun 7.1.26.
    Fix64 t = Fix64(finv(fadd(Fix64::ONE, fmul(C_ERF_P, x.raw()))));
    Fix64 inner = Fix64(C_ERF_A5);
    inner = Fix64(fadd(fmul(inner.raw(), t.raw()), C_ERF_A4));
    inner = Fix64(fadd(fmul(inner.raw(), t.raw()), C_ERF_A3));
    inner = Fix64(fadd(fmul(inner.raw(), t.raw()), C_ERF_A2));
    inner = Fix64(fadd(fmul(inner.raw(), t.raw()), C_ERF_A1));
    Fix64 x2 = Fix64(fmul(x.raw(), x.raw()));
    Fix64 ex = fexp(Fix64(-x2.raw()));
    return Fix64(fmul(fmul(t.raw(), inner.raw()), ex.raw()));
}

void fsincos(Fix64 x, Fix64& sin_out, Fix64& cos_out) {
    // Reduce to [-pi, pi].
    std::int64_t n = fround_int(fdiv(x.raw(), C_TWO_PI));
    Fix64 r = Fix64(fsub(x.raw(), fmul(fromInt(n).raw(), C_TWO_PI)));

    // Reduce to [-pi/2, pi/2], tracking the sign flip (sin/cos(x+pi) = -sin/-cos).
    bool neg = false;
    if (r.raw() > C_HALF_PI) { r = Fix64(fsub(r.raw(), C_PI)); neg = true; }
    else if (r.raw() < -C_HALF_PI) { r = Fix64(fadd(r.raw(), C_PI)); neg = true; }

    // CORDIC rotation mode.
    std::int64_t xv = C_CORDIC_K;
    std::int64_t yv = 0;
    std::int64_t z = r.raw();
    for (int i = 0; i < CORDIC_N; ++i) {
        std::int64_t xnew, ynew;
        if (z >= 0) {
            xnew = xv - (yv >> i);
            ynew = yv + (xv >> i);
            z -= C_ATAN[i];
        } else {
            xnew = xv + (yv >> i);
            ynew = yv - (xv >> i);
            z += C_ATAN[i];
        }
        xv = xnew;
        yv = ynew;
    }

    if (neg) { xv = -xv; yv = -yv; }
    cos_out = Fix64(xv);
    sin_out = Fix64(yv);
}

} // namespace mdsim
