/*
 * END KEM -- AArch64 NEON kernel: complete NTT mod 12289 (Montgomery, R = 2^16).
 * Layers with len >= 8 are whole-vector butterflies; len = 4/2/1 use 64/32/16-bit lane interleaves.
 * Range discipline: centred values, Barrett (|.| <= 6145) every two forward layers and every inverse layer.
 */
#include "poly.h"
#include "neon_common.h"

#define P 12289
#if END_N == 512
#define ZETAS     end_zetas12289_512
#define ZETAS_INV end_zetas12289_inv_512
#define NINV_MONT END_NINV_MONT12289_512
#else
#define ZETAS     end_zetas12289_1024
#define ZETAS_INV end_zetas12289_inv_1024
#define NINV_MONT END_NINV_MONT12289_1024
#endif

static inline void reduce_all(int16_t *a)
{
    const int16x8_t pv = vdupq_n_s16(P), bv = vdupq_n_s16(BARRETT_V_12289);
    unsigned i;
    for (i = 0; i < END_N; i += 8)
        vst1q_s16(a + i, neon_barrett(vld1q_s16(a + i), bv, pv));
}

static inline void ct_layer(int16_t *a, unsigned len, unsigned *k)
{
    const int16x8_t pv = vdupq_n_s16(P), pinv = vdupq_n_s16(PINV16_12289);
    unsigned start, j;
    for (start = 0; start < END_N; start += 2 * len) {
        int16x8_t z = vdupq_n_s16(ZETAS[*k]), zq = vmulq_s16(z, pinv);
        (*k)++;
        for (j = start; j < start + len; j += 8) {
            int16x8_t x = vld1q_s16(a + j), y = vld1q_s16(a + j + len);
            int16x8_t t = neon_montmul(y, z, zq, pv);
            vst1q_s16(a + j + len, vsubq_s16(x, t));
            vst1q_s16(a + j, vaddq_s16(x, t));
        }
    }
}

/* the three intra-vector CT layers (len = 4, 2, 1) on 16 coefficients, with inputs reduced */
static inline void ct_tail(int16_t *a, unsigned k4, unsigned k2, unsigned k1)
{
    const int16x8_t pv = vdupq_n_s16(P), pinv = vdupq_n_s16(PINV16_12289), bv = vdupq_n_s16(BARRETT_V_12289);
    int16x8_t v0 = vld1q_s16(a), v1 = vld1q_s16(a + 8), lo, hi, t, z, zq;
    /* len = 4 */
    z = vcombine_s16(vdup_n_s16(ZETAS[k4]), vdup_n_s16(ZETAS[k4 + 1])); zq = vmulq_s16(z, pinv);
    lo = neon_zip64_lo(v0, v1); hi = neon_zip64_hi(v0, v1);
    t = neon_montmul(hi, z, zq, pv);
    v0 = neon_zip64_lo(vaddq_s16(lo, t), vsubq_s16(lo, t));
    v1 = neon_zip64_hi(vaddq_s16(lo, t), vsubq_s16(lo, t));
    /* len = 2 : L = [a0 a1 a4 a5 b0 b1 b4 b5], H = [a2 a3 a6 a7 b2 b3 b6 b7], twiddles [z z z' z' ...] */
    z = vld1q_s16(ZETAS + k2); z = vzip1q_s16(z, z); zq = vmulq_s16(z, pinv);
    lo = neon_uzp32_1(v0, v1); hi = neon_uzp32_2(v0, v1);
    t = neon_montmul(hi, z, zq, pv);
    v0 = neon_zip32_1(vaddq_s16(lo, t), vsubq_s16(lo, t));
    v1 = neon_zip32_2(vaddq_s16(lo, t), vsubq_s16(lo, t));
    /* reduce before the last layer (values may reach 6145 + 3*12289 otherwise) */
    v0 = neon_barrett(v0, bv, pv); v1 = neon_barrett(v1, bv, pv);
    /* len = 1 : L = even lanes, H = odd lanes, twiddles = 8 consecutive */
    z = vld1q_s16(ZETAS + k1); zq = vmulq_s16(z, pinv);
    lo = vuzp1q_s16(v0, v1); hi = vuzp2q_s16(v0, v1);
    t = neon_montmul(hi, z, zq, pv);
    vst1q_s16(a, vzip1q_s16(vaddq_s16(lo, t), vsubq_s16(lo, t)));
    vst1q_s16(a + 8, vzip2q_s16(vaddq_s16(lo, t), vsubq_s16(lo, t)));
}

void poly_ntt12289(int16_t a[END_N])
{
    const int16x8_t pv = vdupq_n_s16(P), bv = vdupq_n_s16(BARRETT_V_12289);
    unsigned i, len, k = 1, layer = 0;
    reduce_all(a);
    for (len = END_N >> 1; len >= 8; len >>= 1) {
        ct_layer(a, len, &k);
        if ((++layer & 1) == 0)
            reduce_all(a);
    }
    /* here k == N/8 (start of the len = 4 twiddles); len = 2 twiddles start at N/4, len = 1 at N/2 */
    if (layer & 1)
        reduce_all(a);
    for (i = 0; i < END_N; i += 16)
        ct_tail(a + i, END_N / 8 + i / 8, END_N / 4 + i / 4, END_N / 2 + i / 2);
    for (i = 0; i < END_N; i += 8)
        vst1q_s16(a + i, neon_canon(neon_barrett(vld1q_s16(a + i), bv, pv), pv));
}

/* inverse tail: GS layers len = 1, 2, 4 on 16 coefficients; inverse twiddles in consumption order */
static inline void gs_tail(int16_t *a, unsigned i1, unsigned i2, unsigned i4)
{
    const int16x8_t pv = vdupq_n_s16(P), pinv = vdupq_n_s16(PINV16_12289), bv = vdupq_n_s16(BARRETT_V_12289);
    int16x8_t v0 = vld1q_s16(a), v1 = vld1q_s16(a + 8), lo, hi, z, zq, s, d;
    /* len = 1 */
    z = vld1q_s16(ZETAS_INV + i1); zq = vmulq_s16(z, pinv);
    lo = vuzp1q_s16(v0, v1); hi = vuzp2q_s16(v0, v1);
    s = neon_barrett(vaddq_s16(lo, hi), bv, pv); d = neon_montmul(vsubq_s16(lo, hi), z, zq, pv);
    v0 = vzip1q_s16(s, d); v1 = vzip2q_s16(s, d);
    /* len = 2 */
    z = vld1q_s16(ZETAS_INV + i2); z = vzip1q_s16(z, z); zq = vmulq_s16(z, pinv);
    lo = neon_uzp32_1(v0, v1); hi = neon_uzp32_2(v0, v1);
    s = neon_barrett(vaddq_s16(lo, hi), bv, pv); d = neon_montmul(vsubq_s16(lo, hi), z, zq, pv);
    v0 = neon_zip32_1(s, d); v1 = neon_zip32_2(s, d);
    /* len = 4 */
    z = vcombine_s16(vdup_n_s16(ZETAS_INV[i4]), vdup_n_s16(ZETAS_INV[i4 + 1])); zq = vmulq_s16(z, pinv);
    lo = neon_zip64_lo(v0, v1); hi = neon_zip64_hi(v0, v1);
    s = neon_barrett(vaddq_s16(lo, hi), bv, pv); d = neon_montmul(vsubq_s16(lo, hi), z, zq, pv);
    vst1q_s16(a, neon_zip64_lo(s, d));
    vst1q_s16(a + 8, neon_zip64_hi(s, d));
}

void poly_intt12289(int16_t a[END_N])
{
    const int16x8_t pv = vdupq_n_s16(P), pinv = vdupq_n_s16(PINV16_12289), bv = vdupq_n_s16(BARRETT_V_12289);
    unsigned i, len, b, j, idx;
    reduce_all(a);
    /* consumption order of ZETAS_INV: len = 1 -> N/2 entries, len = 2 -> N/4, len = 4 -> N/8, ... */
    for (i = 0; i < END_N; i += 16)
        gs_tail(a + i, i / 2, END_N / 2 + i / 4, END_N / 2 + END_N / 4 + i / 8);
    idx = END_N / 2 + END_N / 4 + END_N / 8;
    for (len = 8; len < END_N; len <<= 1) {
        unsigned nb = END_N / (2 * len);
        for (b = 0; b < nb; b++) {
            int16x8_t z = vdupq_n_s16(ZETAS_INV[idx++]), zq = vmulq_s16(z, pinv);
            unsigned start = 2 * len * b;
            for (j = start; j < start + len; j += 8) {
                int16x8_t x = vld1q_s16(a + j), y = vld1q_s16(a + j + len);
                vst1q_s16(a + j, neon_barrett(vaddq_s16(x, y), bv, pv));
                vst1q_s16(a + j + len, neon_montmul(vsubq_s16(x, y), z, zq, pv));
            }
        }
    }
    {
        int16x8_t z = vdupq_n_s16(NINV_MONT), zq = vmulq_s16(z, pinv), half = vdupq_n_s16(6144);
        for (i = 0; i < END_N; i += 8)
            vst1q_s16(a + i, neon_center_exact(neon_barrett(neon_montmul(vld1q_s16(a + i), z, zq, pv), bv, pv), pv, half));
    }
}

void poly_pointwise12289(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N])
{
    const int16x8_t pv = vdupq_n_s16(P), pinv = vdupq_n_s16(PINV16_12289), bv = vdupq_n_s16(BARRETT_V_12289);
    const int16x8_t r2 = vdupq_n_s16(END_MONT_R2_MOD12289), r2q = vmulq_s16(r2, pinv);
    unsigned i;
    for (i = 0; i < END_N; i += 8) {
        int16x8_t x = vld1q_s16(a + i), y = vld1q_s16(b + i);
        int16x8_t yq = vmulq_s16(y, pinv);
        int16x8_t t = neon_montmul(x, y, yq, pv);          /* a b R^-1 */
        t = neon_montmul(t, r2, r2q, pv);                   /* a b */
        vst1q_s16(r + i, neon_canon(neon_barrett(t, bv, pv), pv));
    }
}
