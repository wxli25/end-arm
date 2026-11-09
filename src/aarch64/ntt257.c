/*
 * END KEM -- AArch64 NEON kernel: NTT mod 257 (reference layout), basemul.
 *
 * Vector arithmetic: 8 x int16 lanes. Twiddles are in Montgomery form (for 257, R = 1 so z_mont = z).
 * Range discipline: forward inputs reduced to [0,257), 7 lazy CT layers (< 8*257), final red257;
 * inverse sums reduced with red257 every layer, differences through Montgomery.
 */
#include "poly.h"
#include "neon_common.h"

#define L END_LEAF

static inline void ct_layer(int16_t *a, unsigned len, const int16_t *zt, unsigned *k)
{
    const int16x8_t pv = vdupq_n_s16(257), pinv = vdupq_n_s16(PINV16_257);
    unsigned start, j;
    for (start = 0; start < END_N; start += 2 * len) {
        int16x8_t z = vdupq_n_s16(zt[*k]), zq = vmulq_s16(z, pinv);
        (*k)++;
        for (j = start; j < start + len; j += 8) {
            int16x8_t x = vld1q_s16(a + j), y = vld1q_s16(a + j + len);
            int16x8_t t = neon_montmul(y, z, zq, pv);
            vst1q_s16(a + j + len, vsubq_s16(x, t));
            vst1q_s16(a + j, vaddq_s16(x, t));
        }
    }
}

void poly_ntt257(int16_t a[END_N])
{
    unsigned i, len, k = 1;
    for (i = 0; i < END_N; i += 8)
        vst1q_s16(a + i, neon_red257(vld1q_s16(a + i)));
    for (len = END_N >> 1; len >= 8; len >>= 1)
        ct_layer(a, len, end_zetas257_mont, &k);
#if L == 4
    {
        const int16x8_t pv = vdupq_n_s16(257), pinv = vdupq_n_s16(PINV16_257);
        for (i = 0; i < END_N; i += 16, k += 2) {
            int16x8_t v0 = vld1q_s16(a + i), v1 = vld1q_s16(a + i + 8);
            int16x8_t zz = vcombine_s16(vdup_n_s16(end_zetas257_mont[k]), vdup_n_s16(end_zetas257_mont[k + 1]));
            int16x8_t zq = vmulq_s16(zz, pinv);
            int16x8_t lo = neon_zip64_lo(v0, v1), hi = neon_zip64_hi(v0, v1);
            int16x8_t t = neon_montmul(hi, zz, zq, pv);
            int16x8_t ln = vaddq_s16(lo, t), hn = vsubq_s16(lo, t);
            vst1q_s16(a + i, neon_zip64_lo(ln, hn));
            vst1q_s16(a + i + 8, neon_zip64_hi(ln, hn));
        }
    }
#endif
    for (i = 0; i < END_N; i += 8)
        vst1q_s16(a + i, neon_red257(vld1q_s16(a + i)));
}

void poly_intt257(int16_t a[END_N])
{
    const int16x8_t pv = vdupq_n_s16(257), pinv = vdupq_n_s16(PINV16_257);
    unsigned i, len, b, j;
    for (i = 0; i < END_N; i += 8)
        vst1q_s16(a + i, neon_red257(vld1q_s16(a + i)));
#if L == 4
    for (i = 0, b = 0; i < END_N; i += 16, b += 2) {
        int16x8_t v0 = vld1q_s16(a + i), v1 = vld1q_s16(a + i + 8);
        int16x8_t zz = vcombine_s16(vdup_n_s16(end_zetas257_inv_mont[64 + b]), vdup_n_s16(end_zetas257_inv_mont[64 + b + 1]));
        int16x8_t zq = vmulq_s16(zz, pinv);
        int16x8_t lo = neon_zip64_lo(v0, v1), hi = neon_zip64_hi(v0, v1);
        int16x8_t ln = neon_red257(vaddq_s16(lo, hi));
        int16x8_t hn = neon_montmul(vsubq_s16(lo, hi), zz, zq, pv);
        vst1q_s16(a + i, neon_zip64_lo(ln, hn));
        vst1q_s16(a + i + 8, neon_zip64_hi(ln, hn));
    }
#endif
    for (len = 8; len <= END_N >> 1; len <<= 1) {
        unsigned nb = END_N / (2 * len);
        for (b = 0; b < nb; b++) {
            int16x8_t z = vdupq_n_s16(end_zetas257_inv_mont[nb + b]), zq = vmulq_s16(z, pinv);
            unsigned start = 2 * len * b;
            for (j = start; j < start + len; j += 8) {
                int16x8_t x = vld1q_s16(a + j), y = vld1q_s16(a + j + len);
                vst1q_s16(a + j, neon_red257(vaddq_s16(x, y)));
                vst1q_s16(a + j + len, neon_montmul(vsubq_s16(x, y), z, zq, pv));
            }
        }
    }
    {
        int16x8_t z = vdupq_n_s16(END_INV128_MONT_MOD257), zq = vmulq_s16(z, pinv);
        for (i = 0; i < END_N; i += 8)
            vst1q_s16(a + i, neon_red257(neon_montmul(vld1q_s16(a + i), z, zq, pv)));
    }
}

/* int32 lanes -> int16 lanes in [0,257): x = hi16*2^16 + lo16 = lo16 + hi16 (mod 257), then Fermat 8-bit step */
static inline int16x8_t red257_from32(int32x4_t lo, int32x4_t hi)
{
    const int32x4_t m16 = vdupq_n_s32(0xFFFF), m8 = vdupq_n_s32(0xFF);
    int32x4_t s0 = vaddq_s32(vandq_s32(lo, m16), vshrq_n_s32(lo, 16));
    int32x4_t s1 = vaddq_s32(vandq_s32(hi, m16), vshrq_n_s32(hi, 16));
    s0 = vsubq_s32(vandq_s32(s0, m8), vshrq_n_s32(s0, 8));           /* [-257, 255] */
    s1 = vsubq_s32(vandq_s32(s1, m8), vshrq_n_s32(s1, 8));
    return neon_red257(vcombine_s16(vmovn_s32(s0), vmovn_s32(s1)));
}

#if L == 4
/* leaves of 4: a0..a3 hold 8 leaves (2 per vector); transpose to coefficient-major A0..A3 */
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

void poly_basemul257(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N])
{
    unsigned i;
    for (i = 0; i < END_N; i += 32) {
        int16x8_t a0 = vld1q_s16(a + i), a1 = vld1q_s16(a + i + 8), a2 = vld1q_s16(a + i + 16), a3 = vld1q_s16(a + i + 24);
        int16x8_t b0 = vld1q_s16(b + i), b1 = vld1q_s16(b + i + 8), b2 = vld1q_s16(b + i + 16), b3 = vld1q_s16(b + i + 24);
        int16x8_t A0, A1, A2, A3, B0, B1, B2, B3, R0, R1, R2, R3;
        int16x4_t z4 = vld1_s16(end_zetas257 + 64 + i / 8);                 /* zetas of leaf pairs */
        int16x8_t zz = vzip1q_s16(vcombine_s16(z4, z4), vnegq_s16(vcombine_s16(z4, z4)));  /* [z0 -z0 z1 -z1 ...] */
        int32x4_t zl = vmovl_s16(vget_low_s16(zz)), zh = vmovl_s16(vget_high_s16(zz));
        TRANSPOSE4(A0, A1, A2, A3, a0, a1, a2, a3);
        TRANSPOSE4(B0, B1, B2, B3, b0, b1, b2, b3);
#define MUL(x, y, part) vmull_s16(vget_##part##_s16(x), vget_##part##_s16(y))
#define MLA(acc, x, y, part) vmlal_s16(acc, vget_##part##_s16(x), vget_##part##_s16(y))
#define LEAFPROD(R, lo_expr, hi_expr_or_null, part, zpart) do { } while (0)
        {
            int32x4_t lo, hi, h;
            /* r0 = a0 b0 + z (a1 b3 + a2 b2 + a3 b1) */
            h = MUL(A1, B3, low); h = MLA(h, A2, B2, low); h = MLA(h, A3, B1, low);
            lo = vmlaq_s32(MUL(A0, B0, low), h, zl);
            h = MUL(A1, B3, high); h = MLA(h, A2, B2, high); h = MLA(h, A3, B1, high);
            hi = vmlaq_s32(MUL(A0, B0, high), h, zh);
            R0 = red257_from32(lo, hi);
            /* r1 = a0 b1 + a1 b0 + z (a2 b3 + a3 b2) */
            h = MUL(A2, B3, low); h = MLA(h, A3, B2, low);
            lo = MUL(A0, B1, low); lo = MLA(lo, A1, B0, low); lo = vmlaq_s32(lo, h, zl);
            h = MUL(A2, B3, high); h = MLA(h, A3, B2, high);
            hi = MUL(A0, B1, high); hi = MLA(hi, A1, B0, high); hi = vmlaq_s32(hi, h, zh);
            R1 = red257_from32(lo, hi);
            /* r2 = a0 b2 + a1 b1 + a2 b0 + z a3 b3 */
            h = MUL(A3, B3, low);
            lo = MUL(A0, B2, low); lo = MLA(lo, A1, B1, low); lo = MLA(lo, A2, B0, low); lo = vmlaq_s32(lo, h, zl);
            h = MUL(A3, B3, high);
            hi = MUL(A0, B2, high); hi = MLA(hi, A1, B1, high); hi = MLA(hi, A2, B0, high); hi = vmlaq_s32(hi, h, zh);
            R2 = red257_from32(lo, hi);
            /* r3 = a0 b3 + a1 b2 + a2 b1 + a3 b0 */
            lo = MUL(A0, B3, low); lo = MLA(lo, A1, B2, low); lo = MLA(lo, A2, B1, low); lo = MLA(lo, A3, B0, low);
            hi = MUL(A0, B3, high); hi = MLA(hi, A1, B2, high); hi = MLA(hi, A2, B1, high); hi = MLA(hi, A3, B0, high);
            R3 = red257_from32(lo, hi);
        }
#undef MUL
#undef MLA
        UNTRANSPOSE4(a0, a1, a2, a3, R0, R1, R2, R3);
        vst1q_s16(r + i, a0); vst1q_s16(r + i + 8, a1); vst1q_s16(r + i + 16, a2); vst1q_s16(r + i + 24, a3);
    }
}
#else
/* leaves of 8 (END-1024): each vector is one leaf; transpose 8 leaves (8x8 int16) to coefficient-major */
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

void poly_basemul257(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N])
{
    unsigned i, k, m;
    for (i = 0; i < END_N; i += 64) {
        int16x8_t A[8], B[8], R[8];
        int16x4_t z4 = vld1_s16(end_zetas257 + 64 + i / 16);
        int16x8_t zz = vzip1q_s16(vcombine_s16(z4, z4), vnegq_s16(vcombine_s16(z4, z4)));
        int32x4_t zl = vmovl_s16(vget_low_s16(zz)), zh = vmovl_s16(vget_high_s16(zz));
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
            lo = vmlaq_s32(lo, hl, zl);
            hi = vmlaq_s32(hi, hh, zh);
            R[k] = red257_from32(lo, hi);
        }
        transpose8(R);
        for (k = 0; k < 8; k++) vst1q_s16(r + i + 8 * k, R[k]);
    }
}
#endif
