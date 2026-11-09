/*
 * Self-contained tests of the L1 arithmetic and L2 polynomial kernels of the selected backend.
 *   - exhaustive int16 tests of the reductions, random tests of Montgomery/Barrett/rounding helpers
 *   - NTT round trips, leaf products vs. schoolbook, inversion, rotation, FCL selection
 * Build: make test_kernels [END_N=512|1024] [BACKEND=ref|aarch64|armv7em]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "end_params.h"
#include "end_arith.h"
#include "poly.h"

static uint64_t rng_state = 0x1234567887654321ULL;
static uint32_t rnd(void)
{
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17;
    return (uint32_t)(rng_state >> 11);
}
static int fails = 0;
#define CHECK(cond, ...) do { if (!(cond)) { fails++; printf("FAIL %s:%d: ", __FILE__, __LINE__); printf(__VA_ARGS__); printf("\n"); } } while (0)

static int64_t modp(int64_t x, int64_t p) { x %= p; if (x < 0) x += p; return x; }
static int64_t modinv(int64_t a, int64_t p) { int64_t r = 1, e = p - 2; a = modp(a, p); while (e) { if (e & 1) r = r * a % p; a = a * a % p; e >>= 1; } return r; }

/* schoolbook negacyclic product mod p */
static void schoolbook(int64_t *r, const int16_t *a, const int16_t *b, int64_t p)
{
    static int64_t t[2 * END_N];
    unsigned i, j;
    memset(t, 0, sizeof t);
    for (i = 0; i < END_N; i++)
        for (j = 0; j < END_N; j++)
            t[i + j] += (int64_t)a[i] * b[j];
    for (i = 0; i < END_N; i++)
        r[i] = modp(t[i] - t[i + END_N], p);
}

static void test_reductions(void)
{
    int x;
    for (x = -32768; x <= 32767; x++) {
        CHECK(red257((int16_t)x) == modp(x, 257), "red257(%d)", x);
        int16_t b = barrett3329((int16_t)x);
        CHECK(b > -1665 && b < 1665 && modp(b, 3329) == modp(x, 3329), "barrett3329(%d)=%d", x, b);
        CHECK(canon3329((int16_t)x) == modp(x, 3329), "canon3329(%d)", x);
        b = barrett12289((int16_t)x);
        CHECK(b >= -6145 && b <= 6145 && modp(b, 12289) == modp(x, 12289), "barrett12289(%d)=%d", x, b);
        CHECK(canon12289((int16_t)x) == modp(x, 12289), "canon12289(%d)", x);
    }
    /* Montgomery: a*R^-1 mod p */
    int64_t rinv3329 = modinv(65536, 3329), rinv12289 = modinv(65536, 12289);
    unsigned i;
    for (i = 0; i < 2000000; i++) {
        int32_t a = (int32_t)rnd();
        a = a % (1 << 28);
        CHECK(modp(mont257(a), 257) == modp(a, 257), "mont257(%d)", (int)a);
        int32_t a3 = a % (32768 * 3329);
        int16_t m = mont3329(a3);
        CHECK(m > -3329 && m < 3329 && modp(m, 3329) == modp((int64_t)a3 * rinv3329, 3329), "mont3329(%d)=%d", (int)a3, (int)m);
        int32_t a12 = (int32_t)(rnd() % (2u * 32768u * 12289u)) - 32768 * 12289;
        m = mont12289(a12);
        CHECK(m > -12289 && m < 12289 && modp(m, 12289) == modp((int64_t)a12 * rinv12289, 12289), "mont12289(%d)=%d", (int)a12, (int)m);
        int32_t a27 = (int32_t)(rnd() & ((1u << 28) - 1)) - (1 << 27);
        CHECK(red3329_32(a27) == modp(a27, 3329), "red3329_32(%d)", (int)a27);
        /* rounding helper */
        int32_t v = (int32_t)(rnd() & ((1u << 28) - 1)) - (1 << 27);
        int64_t want = ((int64_t)v + END_QQP);
        want = (want >= 0) ? want / END_2QQP : -((-want + END_2QQP - 1) / END_2QQP);
        CHECK(round_div_2qqp(v) == want, "round_div_2qqp(%d)=%d want %lld", (int)v, (int)round_div_2qqp(v), (long long)want);
        int32_t c = (int32_t)(rnd() % (2u * END_2QQP)) - END_2QQP;
        int32_t cc = center_2qqp(c);
        CHECK(cc >= -END_QQP && cc < END_QQP && modp(cc, END_2QQP) == modp(c, END_2QQP), "center_2qqp(%d)", (int)c);
    }
    /* edge values of the rounding helper */
    {
        int32_t edges[] = { -END_QQP, -END_QQP - 1, -END_QQP + 1, END_QQP, END_QQP - 1, END_QQP + 1, 0, 1, -1,
                            3 * END_QQP, -3 * END_QQP, 3 * END_QQP - 1, -3 * END_QQP + 1, (1 << 27) - 1, -(1 << 27) };
        for (i = 0; i < sizeof edges / sizeof edges[0]; i++) {
            int32_t v = edges[i];
            int64_t want = ((int64_t)v + END_QQP);
            want = (want >= 0) ? want / END_2QQP : -((-want + END_2QQP - 1) / END_2QQP);
            CHECK(round_div_2qqp(v) == want, "edge round_div_2qqp(%d)", (int)v);
        }
    }
    printf("reductions: done (%d failures so far)\n", fails);
}

static void rand_poly(int16_t *a, int lo, int hi)
{
    unsigned i;
    for (i = 0; i < END_N; i++) a[i] = (int16_t)(lo + (int)(rnd() % (unsigned)(hi - lo + 1)));
}

static void test_ntt(void)
{
    static int16_t a[END_N] END_ALIGN, b[END_N] END_ALIGN, c[END_N] END_ALIGN, d[END_N] END_ALIGN;
    static int64_t ref[END_N];
    unsigned it, i;
    for (it = 0; it < 20; it++) {
        /* --- 257 --- */
        rand_poly(a, -32768, 32767);
        memcpy(b, a, sizeof a);
        poly_ntt257(b);
        for (i = 0; i < END_N; i++) CHECK(b[i] >= 0 && b[i] < 257, "ntt257 range");
        poly_intt257(b);
        for (i = 0; i < END_N; i++) CHECK(b[i] == modp(a[i], 257), "ntt257 roundtrip at %u", i);
        rand_poly(a, 0, 256); rand_poly(b, 0, 256);
        schoolbook(ref, a, b, 257);
        memcpy(c, a, sizeof a); memcpy(d, b, sizeof b);
        poly_ntt257(c); poly_ntt257(d); poly_basemul257(c, c, d); poly_intt257(c);
        for (i = 0; i < END_N; i++) CHECK(c[i] == ref[i], "basemul257 at %u: %d vs %lld", i, c[i], (long long)ref[i]);
        /* inversion: ternary-like polynomial */
        rand_poly(a, -1, 1);
        memcpy(c, a, sizeof a); poly_ntt257(c);
        if (poly_baseinv257(d, c) == 0) {
            poly_basemul257(c, c, d);
            for (i = 0; i < END_N; i++) CHECK(c[i] == ((i % END_LEAF) == 0), "baseinv257 at %u", i);
        } else {
            printf("  (ntt257 inversion: sample not invertible, skipped)\n");
        }
        /* --- 3329 --- */
        rand_poly(a, -32768, 32767);
        memcpy(b, a, sizeof a);
        poly_ntt3329(b);
        for (i = 0; i < END_N; i++) CHECK(b[i] >= 0 && b[i] < 3329, "ntt3329 range");
        poly_intt3329(b);
        for (i = 0; i < END_N; i++) CHECK(b[i] > -1665 && b[i] < 1665 && modp(b[i], 3329) == modp(a[i], 3329), "ntt3329 roundtrip at %u", i);
        rand_poly(a, -40, 40); rand_poly(b, -1, 1);
        schoolbook(ref, a, b, 3329);
        memcpy(c, a, sizeof a); memcpy(d, b, sizeof b);
        poly_ntt3329(c); poly_ntt3329(d); poly_basemul3329(c, c, d); poly_intt3329(c);
        for (i = 0; i < END_N; i++) CHECK(modp(c[i], 3329) == ref[i], "basemul3329 at %u", i);
        /* --- 12289 --- */
        rand_poly(a, -32768, 32767);
        memcpy(b, a, sizeof a);
        poly_ntt12289(b);
        for (i = 0; i < END_N; i++) CHECK(b[i] >= 0 && b[i] < 12289, "ntt12289 range");
        poly_intt12289(b);
        for (i = 0; i < END_N; i++) CHECK(b[i] > -6145 && b[i] < 6145 && modp(b[i], 12289) == modp(a[i], 12289), "ntt12289 roundtrip at %u", i);
        rand_poly(a, -257, 256); rand_poly(b, -32768, 32767);
        schoolbook(ref, a, b, 12289);
        memcpy(c, a, sizeof a); memcpy(d, b, sizeof b);
        poly_ntt12289(c); poly_ntt12289(d); poly_pointwise12289(c, c, d); poly_intt12289(c);
        for (i = 0; i < END_N; i++) CHECK(modp(c[i], 12289) == ref[i], "pointwise12289 at %u", i);
    }
    /* layout fingerprint of NTT-257 (compare with golden/end_golden.py ntt257 of a[i] = (i*7+3) mod 257) */
    for (i = 0; i < END_N; i++) a[i] = (int16_t)((i * 7 + 3) % 257);
    poly_ntt257(a);
    {
        uint32_t fp = 0;
        for (i = 0; i < END_N; i++) fp = fp * 2654435761u + (uint32_t)a[i];
        printf("ntt257 layout fingerprint: %08x\n", fp);
    }
    printf("ntt: done (%d failures so far)\n", fails);
}

static void test_vec(void)
{
    static int16_t a[END_N] END_ALIGN, r[END_N] END_ALIGN;
    static uint32_t sc[END_N];
    unsigned it, i, k;
    for (it = 0; it < 50; it++) {
        rand_poly(a, -32768, 32767);
        k = rnd() % END_N;
        poly_rot16(r, a, k);
        for (i = 0; i < END_N; i++) {
            int16_t want = (i >= k) ? a[i - k] : (int16_t)(-a[i - k + END_N]);
            CHECK(r[i] == want, "rot16 k=%u at %u", k, i);
        }
        uint64_t s = 0;
        for (i = 0; i < END_N; i++) s += (uint64_t)((int32_t)a[i] * a[i]);
        CHECK(poly_sqnorm16(a) == s, "sqnorm16");
        for (i = 0; i < END_N; i++) sc[i] = rnd() % 64;
        uint32_t idx[4];
        fcl_select(idx, 2, sc);
        /* naive: sort by (score, index) */
        unsigned best0 = 0, best1 = 0;
        for (i = 1; i < END_N; i++) if (sc[i] < sc[best0]) best0 = i;
        best1 = (best0 == 0) ? 1 : 0;
        for (i = 0; i < END_N; i++) if (i != best0 && sc[i] < sc[best1]) best1 = i;
        CHECK(idx[0] == best0 && idx[1] == best1, "fcl_select: got (%u,%u) want (%u,%u)", (unsigned)idx[0], (unsigned)idx[1], best0, best1);
    }
    printf("vec: done (%d failures so far)\n", fails);
}

int main(void)
{
    printf("test_kernels %s\n", END_NAME);
    test_reductions();
    test_ntt();
    test_vec();
    printf(fails ? "TEST_KERNELS: %d FAILURES\n" : "TEST_KERNELS: ALL OK\n", fails);
    return fails != 0;
}
