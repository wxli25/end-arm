/*
 * END KEM -- AArch64 NEON kernel: NTT mod 3329 (same 7-layer tree as mod 257) and leaf products.
 * Range discipline: forward inputs centred (|a| <= 1665), 7 lazy layers (< 1665 + 7*3329 < 2^15),
 * canonical output; inverse sums Barrett-reduced every layer; products via widening MACs + Montgomery.
 */
#include "poly.h"
#include "neon_common.h"

#define L END_LEAF
#define P 3329

/* forward and inverse transform: three merged layers per pass (neon_ntt7.h).
 * forward: inputs reduced to |a| <= 1665 when they are read, seven lazy layers (< 1665 + 7 * 3329), canonical output.
 * inverse: sums of layer 7 reduced to |a| <= 1665; after three more layers the four sums (< 8 * 3329) are reduced
 *          again, so the last three layers stay below 8 * 3329 as well; exact centred output. */
static inline int16x8_t barrett3329v(int16x8_t x)
{
    return neon_barrett(x, vdupq_n_s16(BARRETT_V_3329), vdupq_n_s16(P));
}
static inline int16x8_t canon3329v(int16x8_t x)
{
    return neon_canon(barrett3329v(x), vdupq_n_s16(P));
}
static inline int16x8_t center3329v(int16x8_t x)
{
    return neon_center_exact(barrett3329v(x), vdupq_n_s16(P), vdupq_n_s16(1664));
}
#define NTT7_FWD        poly_ntt3329
#define NTT7_INV        poly_intt3329
#define NTT7_P          3329
#define NTT7_PINV       PINV16_3329
#define NTT7_ZETAS      end_zetas3329_mont
#define NTT7_ZETAS_INV  end_zetas3329_inv_mont
#define NTT7_NINV       END_INV128_MONT_MOD3329
#define NTT7_FQMUL(a, b) fqmul3329(a, b)
#define NTT7_RED_IN(x)  barrett3329v(x)
#define NTT7_FWD_OUT(x) canon3329v(x)
#define NTT7_RED_SUM(x) barrett3329v(x)
#define NTT7_RED_MID    1
#define NTT7_INV_OUT(x) center3329v(x)
#include "neon_ntt7.h"

/*
 * leaf product coefficient:  lo + z * hi  with 32-bit accumulators lo, hi (|.| < 2^27)
 *   h16 = montred32(hi) = hi R^-1 ;  z-term = h16 * zR (widening) = hi z ;  total = lo + hi z ;
 *   result = montmul(montred32(total) (= total R^-1), R^2) = total  mod p, canonical.
 */
static inline int16x8_t leafcoef3329(int32x4_t lo0, int32x4_t hi0, int32x4_t lo1, int32x4_t hi1,
                                     int16x4_t zR0, int16x4_t zR1)
{
    const int32x4_t pinv32 = vdupq_n_s32(PINV16_3329);
    const int16x4_t p16 = vdup_n_s16(P);
    const int16x8_t pv = vdupq_n_s16(P), bv = vdupq_n_s16(BARRETT_V_3329);
    const int16x8_t r2 = vdupq_n_s16(END_MONT_R2_MOD3329), r2q = vmulq_s16(r2, vdupq_n_s16(PINV16_3329));
    int16x4_t h0 = neon_montred32(hi0, pinv32, p16), h1 = neon_montred32(hi1, pinv32, p16);
    int32x4_t t0 = vmlal_s16(lo0, h0, zR0), t1 = vmlal_s16(lo1, h1, zR1);
    int16x8_t x = vcombine_s16(neon_montred32(t0, pinv32, p16), neon_montred32(t1, pinv32, p16));
    x = neon_montmul(x, r2, r2q, pv);
    return neon_canon(neon_barrett(x, bv, pv), pv);
}

#if L == 4
#define TRANSPOSE4(A0, A1, A2, A3, a0, a1, a2, a3) do {                 \
        int16x8_t t0 = vuzp1q_s16(a0, a1), t1 = vuzp2q_s16(a0, a1);      \
        int16x8_t t2 = vuzp1q_s16(a2, a3), t3 = vuzp2q_s16(a2, a3);      \
        A0 = vuzp1q_s16(t0, t2); A2 = vuzp2q_s16(t0, t2);                \
        A1 = vuzp1q_s16(t1, t3); A3 = vuzp2q_s16(t1, t3);                \
    } while (0)
#define UNTRANSPOSE4(a0, a1, a2, a3, A0, A1, A2, A3) do {               \
        int16x8_t t0 = vzip1q_s16(A0, A2), t2 = vzip2q_s16(A0, A2);      \
        int16x8_t t1 = vzip1q_s16(A1, A3), t3 = vzip2q_s16(A1, A3);      \
        a0 = vzip1q_s16(t0, t1); a1 = vzip2q_s16(t0, t1);                \
        a2 = vzip1q_s16(t2, t3); a3 = vzip2q_s16(t2, t3);                \
    } while (0)
#define MUL(x, y, part) vmull_s16(vget_##part##_s16(x), vget_##part##_s16(y))
#define MLA(acc, x, y, part) vmlal_s16(acc, vget_##part##_s16(x), vget_##part##_s16(y))

void poly_basemul3329(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N])
{
    unsigned i;
    for (i = 0; i < END_N; i += 32) {
        int16x8_t a0 = vld1q_s16(a + i), a1 = vld1q_s16(a + i + 8), a2 = vld1q_s16(a + i + 16), a3 = vld1q_s16(a + i + 24);
        int16x8_t b0 = vld1q_s16(b + i), b1 = vld1q_s16(b + i + 8), b2 = vld1q_s16(b + i + 16), b3 = vld1q_s16(b + i + 24);
        int16x8_t A0, A1, A2, A3, B0, B1, B2, B3, R0, R1, R2, R3;
        int16x4_t z4 = vld1_s16(end_zetas3329_mont + 64 + i / 8);
        int16x8_t zz = vzip1q_s16(vcombine_s16(z4, z4), vnegq_s16(vcombine_s16(z4, z4)));   /* zeta R, [z0 -z0 z1 -z1 ..] */
        int16x4_t zl = vget_low_s16(zz), zh = vget_high_s16(zz);
        int32x4_t lo, hi, hl, hh;
        TRANSPOSE4(A0, A1, A2, A3, a0, a1, a2, a3);
        TRANSPOSE4(B0, B1, B2, B3, b0, b1, b2, b3);
        /* r0 */
        lo = MUL(A0, B0, low); hi = MUL(A0, B0, high);
        hl = MUL(A1, B3, low); hl = MLA(hl, A2, B2, low); hl = MLA(hl, A3, B1, low);
        hh = MUL(A1, B3, high); hh = MLA(hh, A2, B2, high); hh = MLA(hh, A3, B1, high);
        R0 = leafcoef3329(lo, hl, hi, hh, zl, zh);
        /* r1 */
        lo = MUL(A0, B1, low); lo = MLA(lo, A1, B0, low); hi = MUL(A0, B1, high); hi = MLA(hi, A1, B0, high);
        hl = MUL(A2, B3, low); hl = MLA(hl, A3, B2, low); hh = MUL(A2, B3, high); hh = MLA(hh, A3, B2, high);
        R1 = leafcoef3329(lo, hl, hi, hh, zl, zh);
        /* r2 */
        lo = MUL(A0, B2, low); lo = MLA(lo, A1, B1, low); lo = MLA(lo, A2, B0, low);
        hi = MUL(A0, B2, high); hi = MLA(hi, A1, B1, high); hi = MLA(hi, A2, B0, high);
        hl = MUL(A3, B3, low); hh = MUL(A3, B3, high);
        R2 = leafcoef3329(lo, hl, hi, hh, zl, zh);
        /* r3 */
        lo = MUL(A0, B3, low); lo = MLA(lo, A1, B2, low); lo = MLA(lo, A2, B1, low); lo = MLA(lo, A3, B0, low);
        hi = MUL(A0, B3, high); hi = MLA(hi, A1, B2, high); hi = MLA(hi, A2, B1, high); hi = MLA(hi, A3, B0, high);
        hl = vdupq_n_s32(0); hh = vdupq_n_s32(0);
        R3 = leafcoef3329(lo, hl, hi, hh, zl, zh);
        UNTRANSPOSE4(a0, a1, a2, a3, R0, R1, R2, R3);
        vst1q_s16(r + i, a0); vst1q_s16(r + i + 8, a1); vst1q_s16(r + i + 16, a2); vst1q_s16(r + i + 24, a3);
    }
}
#else
static inline void transpose8(int16x8_t v[8])
{
    int16x8_t t[8];
    unsigned i;
    for (i = 0; i < 8; i += 2) { t[i] = vtrn1q_s16(v[i], v[i + 1]); t[i + 1] = vtrn2q_s16(v[i], v[i + 1]); }
    for (i = 0; i < 8; i += 4) {
        v[i]     = vreinterpretq_s16_s32(vtrn1q_s32(vreinterpretq_s32_s16(t[i]),     vreinterpretq_s32_s16(t[i + 2])));
        v[i + 2] = vreinterpretq_s16_s32(vtrn2q_s32(vreinterpretq_s32_s16(t[i]),     vreinterpretq_s32_s16(t[i + 2])));
        v[i + 1] = vreinterpretq_s16_s32(vtrn1q_s32(vreinterpretq_s32_s16(t[i + 1]), vreinterpretq_s32_s16(t[i + 3])));
        v[i + 3] = vreinterpretq_s16_s32(vtrn2q_s32(vreinterpretq_s32_s16(t[i + 1]), vreinterpretq_s32_s16(t[i + 3])));
    }
    for (i = 0; i < 4; i++) {
        t[i]     = vreinterpretq_s16_s64(vtrn1q_s64(vreinterpretq_s64_s16(v[i]), vreinterpretq_s64_s16(v[i + 4])));
        t[i + 4] = vreinterpretq_s16_s64(vtrn2q_s64(vreinterpretq_s64_s16(v[i]), vreinterpretq_s64_s16(v[i + 4])));
    }
    for (i = 0; i < 8; i++) v[i] = t[i];
}

void poly_basemul3329(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N])
{
    unsigned i, k, m;
    for (i = 0; i < END_N; i += 64) {
        int16x8_t A[8], B[8], R[8];
        int16x4_t z4 = vld1_s16(end_zetas3329_mont + 64 + i / 16);
        int16x8_t zz = vzip1q_s16(vcombine_s16(z4, z4), vnegq_s16(vcombine_s16(z4, z4)));
        int16x4_t zl = vget_low_s16(zz), zh = vget_high_s16(zz);
        for (k = 0; k < 8; k++) { A[k] = vld1q_s16(a + i + 8 * k); B[k] = vld1q_s16(b + i + 8 * k); }
        transpose8(A); transpose8(B);
        for (k = 0; k < 8; k++) {
            int32x4_t lo = vdupq_n_s32(0), hi = vdupq_n_s32(0), hl = vdupq_n_s32(0), hh = vdupq_n_s32(0);
            for (m = 0; m < 8; m++) {
                unsigned j = (k + 8 - m) & 7;
                if (m <= k) {
                    lo = vmlal_s16(lo, vget_low_s16(A[m]), vget_low_s16(B[j]));
                    hi = vmlal_s16(hi, vget_high_s16(A[m]), vget_high_s16(B[j]));
                } else {
                    hl = vmlal_s16(hl, vget_low_s16(A[m]), vget_low_s16(B[j]));
                    hh = vmlal_s16(hh, vget_high_s16(A[m]), vget_high_s16(B[j]));
                }
            }
            R[k] = leafcoef3329(lo, hl, hi, hh, zl, zh);
        }
        transpose8(R);
        for (k = 0; k < 8; k++) vst1q_s16(r + i + 8 * k, R[k]);
    }
}
#endif
