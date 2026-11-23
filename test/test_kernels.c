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
        rand_poly(a, -16384, 16384);
        memcpy(b, a, sizeof a);
        poly_ntt257(b);
        for (i = 0; i < END_N; i++) CHECK(b[i] >= 0 && b[i] < 257, "ntt257 range");
        poly_intt257(b);
        for (i = 0; i < END_N; i++) CHECK(b[i] == modp(a[i], 257), "ntt257 roundtrip at %u", i);
        rand_poly(a, 0, 256); rand_poly(b, 0, 256);
        schoolbook(ref, a, b, 257);
        memcpy(c, a, sizeof a); memcpy(d, b, sizeof b);
        poly_ntt257(c); poly_ntt257(d);
        int16_t *product = (it & 1) ? d : c;     /* callers overwrite either operand */
        poly_basemul257(product, c, d); poly_intt257(product);
        for (i = 0; i < END_N; i++) CHECK(product[i] == ref[i], "basemul257 at %u: %d vs %lld", i, product[i], (long long)ref[i]);
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
        rand_poly(a, -16384, 16384);
        memcpy(b, a, sizeof a);
        poly_ntt3329(b);
#ifndef END_M4_ASM
        for (i = 0; i < END_N; i++) CHECK(b[i] >= 0 && b[i] < 3329, "ntt3329 range");
#endif
        for (i = 0; i < END_N; i++) b[i] = canon3329(b[i]);     /* the NTT-3329 representation is backend-defined */
        poly_intt3329(b);
        for (i = 0; i < END_N; i++) CHECK(b[i] > -1665 && b[i] < 1665 && modp(b[i], 3329) == modp(a[i], 3329), "ntt3329 roundtrip at %u", i);
        rand_poly(a, -40, 40); rand_poly(b, -1, 1);
        schoolbook(ref, a, b, 3329);
        memcpy(c, a, sizeof a); memcpy(d, b, sizeof b);
        poly_ntt3329(c); poly_ntt3329(d);
        poly_basemul3329(product, c, d); poly_intt3329(product);
        for (i = 0; i < END_N; i++) CHECK(modp(product[i], 3329) == ref[i], "basemul3329 at %u", i);
        /* --- 12289 --- */
        rand_poly(a, -8191, 8191);
        memcpy(b, a, sizeof a);
        poly_ntt12289(b);
#ifndef END_M4_ASM
        for (i = 0; i < END_N; i++) CHECK(b[i] >= 0 && b[i] < 12289, "ntt12289 range");
#endif
        for (i = 0; i < END_N; i++) b[i] = barrett12289(b[i]);   /* the NTT-12289 representation is backend-defined */
        poly_intt12289(b);
        for (i = 0; i < END_N; i++) CHECK(b[i] > -6145 && b[i] < 6145 && modp(b[i], 12289) == modp(a[i], 12289), "ntt12289 roundtrip at %u", i);
        rand_poly(a, -257, 256); rand_poly(b, -8191, 8191);
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
        printf("ntt257 layout fingerprint: %08x\n", (unsigned)fp);
    }
    printf("ntt: done (%d failures so far)\n", fails);
}

/* ---- exact reference of the 7-layer transforms and of the leaf products (layout of src/ref) ---- */
static void ref_ntt7(int64_t *r, const int16_t *a, const int16_t *zetas, int64_t p)
{
    unsigned len, start, j, k = 1;
    for (j = 0; j < END_N; j++) r[j] = modp(a[j], p);
    for (len = END_N >> 1; len >= END_LEAF; len >>= 1)
        for (start = 0; start < END_N; start += 2 * len) {
            int64_t z = zetas[k++];
            for (j = start; j < start + len; j++) {
                int64_t t = z * r[j + len] % p;
                r[j + len] = modp(r[j] - t, p);
                r[j] = (r[j] + t) % p;
            }
        }
}

static void ref_basemul7(int64_t *r, const int16_t *a, const int16_t *b, const int16_t *zetas, int64_t p)
{
    unsigned blk, i, k;
    for (blk = 0; blk < END_N / END_LEAF; blk++) {
        int64_t z = zetas[64 + (blk >> 1)];
        if (blk & 1) z = p - z;
        for (k = 0; k < END_LEAF; k++) {
            int64_t lo = 0, hi = 0;
            for (i = 0; i < END_LEAF; i++) {
                if (i <= k) lo += (int64_t)a[blk * END_LEAF + i] * b[blk * END_LEAF + k - i];
                else        hi += (int64_t)a[blk * END_LEAF + i] * b[blk * END_LEAF + k + END_LEAF - i];
            }
            r[blk * END_LEAF + k] = modp(lo + z * modp(hi, p), p);
        }
    }
}

/* exact complete transform mod 12289: zeta_k = psi^brv(k), layout of src/ref */
static void ref_ntt12289(int64_t *r, const int16_t *a)
{
    static int32_t zet[END_N], pw[END_N];
    static int ready = 0;
    const int64_t p = 12289, psi = (END_N == 512) ? END_PSI12289_512 : END_PSI12289_1024;
    unsigned len, start, j, k, b;
    if (!ready) {
        pw[0] = 1;
        for (k = 1; k < END_N; k++) pw[k] = (int32_t)(pw[k - 1] * psi % p);
        for (k = 0; k < END_N; k++) {
            unsigned e = 0;
            for (b = 0; b < END_LOGN; b++) e |= ((k >> b) & 1u) << (END_LOGN - 1 - b);
            zet[k] = pw[e];
        }
        ready = 1;
    }
    for (j = 0; j < END_N; j++) r[j] = modp(a[j], p);
    k = 1;
    for (len = END_N >> 1; len >= 1; len >>= 1)
        for (start = 0; start < END_N; start += 2 * len) {
            int64_t z = zet[k++];
            for (j = start; j < start + len; j++) {
                int64_t t = z * r[j + len] % p;
                r[j + len] = modp(r[j] - t, p);
                r[j] = (r[j] + t) % p;
            }
        }
}

/* input patterns that maximise the growth inside lazily reduced transforms */
#define NPATTERN 11
static void pattern_poly(int16_t *a, unsigned pat, int bound)
{
    unsigned i;
    for (i = 0; i < END_N; i++) {
        int v;
        switch (pat) {
        case 0: v = bound; break;
        case 1: v = -bound; break;
        case 2: v = (i & 1) ? -bound : bound; break;
        case 3: v = ((i / END_LEAF) & 1) ? -bound : bound; break;
        case 4: v = (i < END_N / 2) ? bound : -bound; break;
        case 5: v = (rnd() & 1) ? bound : -bound; break;
        case 6: v = (i == 0) ? bound : 0; break;
        case 7: v = -bound + (int)(rnd() % (unsigned)(2 * bound + 1)); break;
        default: v = ((i / (END_LEAF << (pat - 7))) & 1) ? -bound : bound; break;   /* square waves over 2, 4, 8 leaves */
        }
        a[i] = (int16_t)v;
    }
}

static void test_ntt_exact(void)
{
    static int16_t a[END_N] END_ALIGN, b[END_N] END_ALIGN, c[END_N] END_ALIGN;
    static int64_t ref[END_N];
    unsigned pat, it, i;
    for (it = 0; it < 3; it++) {
        for (pat = 0; pat < NPATTERN; pat++) {
            /* forward transforms at the input bound of the contract (|a| <= 16384): values and, mod 257, layout */
            pattern_poly(a, pat, it == 0 ? 16384 : (it == 1 ? 3329 : 257));
            memcpy(b, a, sizeof a);
            poly_ntt257(b);
            ref_ntt7(ref, a, end_zetas257, 257);
            for (i = 0; i < END_N; i++) CHECK(b[i] == ref[i], "ntt257 exact pat %u at %u: %d vs %lld", pat, i, b[i], (long long)ref[i]);
            memcpy(b, a, sizeof a);
            poly_ntt3329(b);
            ref_ntt7(ref, a, end_zetas3329, 3329);
            for (i = 0; i < END_N; i++) CHECK(modp(b[i], 3329) == ref[i], "ntt3329 exact pat %u at %u: %d vs %lld", pat, i, b[i], (long long)ref[i]);
            /* inverse transforms at their input bound (|a| <= 4095): forward(inverse(x)) == x */
            pattern_poly(a, pat, it == 0 ? 4095 : (it == 1 ? 3328 : 256));
            memcpy(b, a, sizeof a);
            poly_intt257(b);
            for (i = 0; i < END_N; i++) CHECK(b[i] >= 0 && b[i] < 257, "intt257 range pat %u", pat);
            ref_ntt7(ref, b, end_zetas257, 257);
            for (i = 0; i < END_N; i++) CHECK(ref[i] == modp(a[i], 257), "intt257 exact pat %u at %u", pat, i);
            memcpy(b, a, sizeof a);
            poly_intt3329(b);
            for (i = 0; i < END_N; i++) CHECK(b[i] > -1665 && b[i] < 1665, "intt3329 range pat %u", pat);
            ref_ntt7(ref, b, end_zetas3329, 3329);
            for (i = 0; i < END_N; i++) CHECK(ref[i] == modp(a[i], 3329), "intt3329 exact pat %u at %u", pat, i);
            /* leaf products on canonical extremes */
            pattern_poly(a, pat, 128); pattern_poly(b, (pat + 3) % NPATTERN, 128);
            for (i = 0; i < END_N; i++) { a[i] = red257(a[i]); b[i] = red257(b[i]); }
            ref_basemul7(ref, a, b, end_zetas257, 257);
            poly_basemul257(c, a, b);
            for (i = 0; i < END_N; i++) CHECK(c[i] == ref[i], "basemul257 exact pat %u at %u: %d vs %lld", pat, i, c[i], (long long)ref[i]);
            pattern_poly(a, pat, 1664); pattern_poly(b, (pat + 3) % NPATTERN, 1664);
            for (i = 0; i < END_N; i++) { a[i] = canon3329(a[i]); b[i] = canon3329(b[i]); }
            ref_basemul7(ref, a, b, end_zetas3329, 3329);
            poly_basemul3329(c, a, b);
            for (i = 0; i < END_N; i++) CHECK(modp(c[i], 3329) == ref[i], "basemul3329 exact pat %u at %u", pat, i);
            /* complete transform mod 12289 at its input bound (|a| <= 8191), forward and inverse */
            pattern_poly(a, pat, it == 0 ? 8191 : (it == 1 ? 6144 : 257));
            memcpy(b, a, sizeof a);
            poly_ntt12289(b);
            ref_ntt12289(ref, a);
            for (i = 0; i < END_N; i++) CHECK(modp(b[i], 12289) == ref[i], "ntt12289 exact pat %u at %u: %d vs %lld", pat, i, b[i], (long long)ref[i]);
            memcpy(b, a, sizeof a);
            poly_intt12289(b);
            for (i = 0; i < END_N; i++) CHECK(b[i] > -6145 && b[i] < 6145, "intt12289 range pat %u", pat);
            ref_ntt12289(ref, b);
            for (i = 0; i < END_N; i++) CHECK(ref[i] == modp(a[i], 12289), "intt12289 exact pat %u at %u", pat, i);
            pattern_poly(a, pat, 6144); pattern_poly(b, (pat + 3) % NPATTERN, 6144);
            for (i = 0; i < END_N; i++) { a[i] = canon12289(a[i]); b[i] = canon12289(b[i]); }
            poly_pointwise12289(c, a, b);
            for (i = 0; i < END_N; i++) CHECK(modp(c[i], 12289) == modp((int64_t)a[i] * b[i], 12289), "pointwise12289 exact pat %u at %u", pat, i);
#ifdef END_M4_ASM
            pattern_poly(a, pat, 32767); pattern_poly(b, (pat + 5) % NPATTERN, 32767);
            if (pat == 1) { for (i = 0; i < END_N; i++) { a[i] = -32768; b[i] = -32768; } }
            memcpy(c, a, sizeof a);
            poly_pointwise12289(c, c, b);
            for (i = 0; i < END_N; i++) CHECK(c[i] >= -6144 && c[i] <= 6144 && modp(c[i], 12289) == modp((int64_t)a[i] * b[i], 12289), "pointwise12289 lazy pat %u at %u", pat, i);
#endif
#ifdef END_M4_ASM
            /* this backend accepts any int16 in the leaf products (lazy NTT outputs); results alias the inputs */
            pattern_poly(a, pat, 32767); pattern_poly(b, (pat + 5) % NPATTERN, 32767);
            if (pat == 1) { for (i = 0; i < END_N; i++) { a[i] = -32768; b[i] = -32768; } }
            ref_basemul7(ref, a, b, end_zetas3329, 3329);
            memcpy(c, a, sizeof a);
            poly_basemul3329(c, c, b);
            for (i = 0; i < END_N; i++) CHECK(c[i] >= -1664 && c[i] <= 1664 && modp(c[i], 3329) == ref[i], "basemul3329 lazy pat %u at %u", pat, i);
            ref_basemul7(ref, a, b, end_zetas257, 257);
            memcpy(c, b, sizeof b);
            poly_basemul257(c, a, c);
            for (i = 0; i < END_N; i++) CHECK(c[i] == ref[i], "basemul257 lazy pat %u at %u", pat, i);
#endif
        }
    }
    printf("ntt exact: done (%d failures so far)\n", fails);
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
        /* FCL: scores with many ties (range 0..1023 over END_N entries); iteration 0 is the crafted case
         * "two tied second-best entries arrive before a unique minimum", which a score-only insertion
         * gets wrong (it keeps index 1 instead of index 0) */
        for (i = 0; i < END_N; i++) sc[i] = rnd() % 1024;
        if (it == 0) { for (i = 0; i < END_N; i++) sc[i] = 100; sc[0] = 5; sc[1] = 5; sc[2] = 3; }
        {
            unsigned m0, j;
            uint32_t idx[END_FCL_MAX], want[END_FCL_MAX];
            for (m0 = 1; m0 <= END_FCL_MAX; m0++) {
                fcl_select(idx, m0, sc);
                /* naive: repeatedly take the smallest (score, index) not taken yet */
                for (j = 0; j < m0; j++) {
                    unsigned b = END_N, t;
                    for (i = 0; i < END_N; i++) {
                        unsigned taken = 0;
                        for (t = 0; t < j; t++) taken |= (want[t] == i);
                        if (!taken && (b == END_N || sc[i] < sc[b])) b = i;
                    }
                    want[j] = b;
                    CHECK(idx[j] == b, "fcl_select m0=%u pos %u: got %u want %u", m0, j, (unsigned)idx[j], b);
                }
            }
        }
    }
    printf("vec: done (%d failures so far)\n", fails);
}

/* ---- inner loops of the decoder against direct formulas ---- */
static void test_deckern(void)
{
    static int16_t a[END_N] END_ALIGN, b[END_N] END_ALIGN, c[END_N] END_ALIGN, d[END_N] END_ALIGN, r[END_N] END_ALIGN;
    static int32_t C2[END_N] END_ALIGN, ct[END_N] END_ALIGN, cu[END_N] END_ALIGN;
    static const int sg[9][2] = { {0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1}, {-1, -1}, {-1, 1}, {1, -1}, {1, 1} };
    const int64_t inv257 = modinv(257, 12289);
    unsigned it, i;
    for (it = 0; it < 36; it++) {
        int s0 = sg[it % 9][0], s1 = sg[it % 9][1];
        uint16_t par = (uint16_t)(it & 1);
        /* lift */
        rand_poly(a, 0, 256); rand_poly(b, 0, 256);
        poly_lift257(r, a, b, par);
        for (i = 0; i < END_N; i++) {
            int64_t x = modp(2 * a[i] - b[i], 257);
            if ((x ^ par) & 1) x += 257;
            if (x >= 257) x -= 514;
            CHECK(r[i] == x, "lift257 at %u: %d vs %d", i, r[i], (int)x);
        }
        /* CRT */
        rand_poly(c, -6144, 6144); rand_poly(d, -257, 256);
        if (it == 0) { for (i = 0; i < END_N; i++) { a[i] = 256; b[i] = 0; c[i] = (i & 1) ? 6144 : -6144; d[i] = (i & 2) ? 256 : -257; } }
        if (it == 1 || it == 2) {
            /* Force Q*k - t12289 to either side of the centring boundary,
             * including both exact ties and both parities of the lift. */
            static const int16_t residues[] = { -6144, -1, 0, 1, 6144 };
            for (i = 0; i < END_N; i++) {
                int k = (int)(i % 3) + 256;
                a[i] = 0;
                c[i] = residues[(i / 3) % 5];
                d[i] = (int16_t)((unsigned)(k + c[i] + par) & 1u);
                b[i] = (int16_t)modp((int64_t)END_QP * (k - d[i]) - c[i], END_Q);
            }
        }
        poly_crt_c2(C2, a, c, d, b, par);
        for (i = 0; i < END_N; i++) {
            int64_t u257 = modp(2 * a[i] - b[i], 257), u12 = modp(c[i], 12289);
            int64_t u = u257 + 257 * modp((u12 - u257) * inv257, 12289), v;
            if ((u ^ par) & 1) u += END_QQP;
            v = modp((int64_t)END_QP * d[i] - u + END_QQP, END_2QQP) - END_QQP;
            CHECK(C2[i] == v, "crt_c2 at %u: %d vs %d", i, (int)C2[i], (int)v);
        }
        /* candidate rounding: C2 in [-qQ, qQ), w any int16 */
        rand_poly(a, -32768, 32767); rand_poly(b, -32768, 32767);
        for (i = 0; i < END_N; i++) C2[i] = (int32_t)(rnd() % (uint32_t)END_2QQP) - END_QQP;
        if (it < 9) { for (i = 0; i < END_N; i++) { C2[i] = (i & 1) ? END_QQP - 1 : -END_QQP; a[i] = (i & 2) ? 32767 : -32768; b[i] = (i & 4) ? 32767 : -32768; } }
        poly_cand_round(r, ct, C2, a, b, s0, s1);
        for (i = 0; i < END_N; i++) {
            int64_t v = (int64_t)C2[i] - 2 * END_Q * ((int64_t)s0 * a[i] + (int64_t)s1 * b[i]);
            int64_t q = (v + END_QQP >= 0) ? (v + END_QQP) / END_2QQP : -((-(v + END_QQP) + END_2QQP - 1) / END_2QQP);
            CHECK(r[i] == q && ct[i] == v - q * END_2QQP, "cand_round at %u: %d %d vs %d %d", i, r[i], (int)ct[i], (int)q, (int)(v - q * END_2QQP));
        }
        /* candidate combination and norm */
        rand_poly(a, -2047, 2047); rand_poly(b, -127, 127); rand_poly(c, -127, 127); rand_poly(d, -2047, 2047);
        {
            uint64_t want = 0, got = poly_cand_combine(r, a, b, c, d, s0, s1);
            for (i = 0; i < END_N; i++) {
                int64_t x = a[i] + 2 * (s0 * b[i] + s1 * c[i]) + 2 * d[i];
                CHECK(r[i] == x, "cand_combine at %u", i);
                want += (uint64_t)(x * x);
            }
            CHECK(got == want, "cand_combine norm");
        }
        /* +-1 check: one polynomial that passes, one with a single wrong coefficient, one random */
        rand_poly(b, -127, 127); rand_poly(c, -127, 127);
        for (i = 0; i < END_N; i++) a[i] = (int16_t)(((rnd() & 1) ? 1 : -1) + 2 * (s0 * b[i] + s1 * c[i]));
        CHECK(poly_cand_check(r, a, b, c, s0, s1) == 0, "cand_check pass");
        for (i = 0; i < END_N; i++) CHECK(r[i] == a[i] - 2 * (s0 * b[i] + s1 * c[i]), "cand_check value at %u", i);
        i = rnd() % END_N;
        a[i] = (int16_t)(a[i] + ((it & 1) ? 4 : -4));          /* +-1 becomes +-3 or +-5 */
        CHECK(poly_cand_check(r, a, b, c, s0, s1) != 0, "cand_check fail");
        rand_poly(a, -16383, 16383);
        CHECK(poly_cand_check(r, a, b, c, s0, s1) != 0, "cand_check random");
        for (i = 0; i < END_N; i++) CHECK(r[i] == a[i] - 2 * (s0 * b[i] + s1 * c[i]), "cand_check value2 at %u", i);
        /* masked copies */
        for (i = 0; i < END_N; i++) { ct[i] = (int32_t)rnd() - (int32_t)rnd(); cu[i] = (int32_t)rnd(); C2[i] = cu[i]; }
        poly_cmov32(cu, ct, (it & 1) ? 0xFFFFFFFFu : 0u);
        for (i = 0; i < END_N; i++) CHECK(cu[i] == ((it & 1) ? ct[i] : C2[i]), "cmov32 at %u", i);
        rand_poly(a, -32768, 32767); rand_poly(b, -32768, 32767);
        memcpy(c, b, sizeof b);
        poly_cmov16(b, a, (uint16_t)((it & 1) ? 0xFFFFu : 0u));
        for (i = 0; i < END_N; i++) CHECK(b[i] == ((it & 1) ? a[i] : c[i]), "cmov16 at %u", i);
    }
    printf("decoder kernels: done (%d failures so far)\n", fails);
}

int main(void)
{
    printf("test_kernels %s\n", END_NAME);
    test_reductions();
    test_ntt();
    test_ntt_exact();
    test_vec();
    test_deckern();
    printf(fails ? "TEST_KERNELS: %d FAILURES\n" : "TEST_KERNELS: ALL OK\n", fails);
    return fails != 0;
}
