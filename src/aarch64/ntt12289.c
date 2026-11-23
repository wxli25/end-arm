/*
 * END KEM -- AArch64 NEON kernel: complete NTT mod 12289 (Montgomery, R = 2^16).
 * Layers with len >= 8 are whole-vector butterflies; len = 4/2/1 use 64/32/16-bit lane interleaves.
 * Two wide CT layers share four vectors. Their outputs are Barrett-reduced before storage.
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

/* The two layer lengths are 2*stride and stride. Their twiddles form one subtree
 * (k, 2k, 2k+1). Each of the two additions grows the magnitude by less than P,
 * so a centred starting value stays within 6145 + 2*(P-1) = 30721. */
static inline void ct_layer_pair(int16_t *a, unsigned stride)
{
    const int16x8_t pv = vdupq_n_s16(P), pinv = vdupq_n_s16(PINV16_12289);
    const int16x8_t bv = vdupq_n_s16(BARRETT_V_12289);
    unsigned k = END_N / (4 * stride);
    for (unsigned start = 0; start < END_N; start += 4 * stride, k++) {
        int16x8_t z = vdupq_n_s16(ZETAS[k]), zq = vmulq_s16(z, pinv);
        int16x8_t zl = vdupq_n_s16(ZETAS[2 * k]), zlq = vmulq_s16(zl, pinv);
        int16x8_t zr = vdupq_n_s16(ZETAS[2 * k + 1]), zrq = vmulq_s16(zr, pinv);
        for (unsigned j = 0; j < stride; j += 8) {
            int16_t *p = a + start + j;
            int16x8_t x0 = vld1q_s16(p), x1 = vld1q_s16(p + stride);
            int16x8_t t0 = neon_montmul(vld1q_s16(p + 2 * stride), z, zq, pv);
            int16x8_t t1 = neon_montmul(vld1q_s16(p + 3 * stride), z, zq, pv);
            if (stride == END_N / 4) {
                /* The other two inputs have already passed through a Montgomery product. */
                x0 = neon_barrett(x0, bv, pv);
                x1 = neon_barrett(x1, bv, pv);
            }
            int16x8_t l = vaddq_s16(x0, t0), r = vsubq_s16(x0, t0);
            int16x8_t tl = neon_montmul(vaddq_s16(x1, t1), zl, zlq, pv);
            int16x8_t tr = neon_montmul(vsubq_s16(x1, t1), zr, zrq, pv);
            vst1q_s16(p,              neon_barrett(vaddq_s16(l, tl), bv, pv));
            vst1q_s16(p + stride,     neon_barrett(vsubq_s16(l, tl), bv, pv));
            vst1q_s16(p + 2 * stride, neon_barrett(vaddq_s16(r, tr), bv, pv));
            vst1q_s16(p + 3 * stride, neon_barrett(vsubq_s16(r, tr), bv, pv));
        }
    }
}

/* Short CT layers on a group of 16 coefficients starting at offset i. */
static inline void ct_tail(int16_t *a, unsigned i)
{
    const unsigned k4 = END_N / 8 + i / 8, k2 = END_N / 4 + i / 4, k1 = END_N / 2 + i / 2;
    const int16x8_t pv = vdupq_n_s16(P), pinv = vdupq_n_s16(PINV16_12289), bv = vdupq_n_s16(BARRETT_V_12289);
    int16x8_t v0 = vld1q_s16(a), v1 = vld1q_s16(a + 8), lo, hi, t, z, zq;
#if END_N == 1024
    /* Seven wide layers: the unpaired len=8 butterfly shares the tail's load. */
    z = vdupq_n_s16(ZETAS[END_N / 16 + i / 16]);
    t = neon_montmul(v1, z, vmulq_s16(z, pinv), pv);
    v1 = neon_barrett(vsubq_s16(v0, t), bv, pv);
    v0 = neon_barrett(vaddq_s16(v0, t), bv, pv);
#endif
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
    v0 = neon_canon(neon_barrett(vaddq_s16(lo, t), bv, pv), pv);
    v1 = neon_canon(neon_barrett(vsubq_s16(lo, t), bv, pv), pv);
    vst1q_s16(a, vzip1q_s16(v0, v1));
    vst1q_s16(a + 8, vzip2q_s16(v0, v1));
}

void poly_ntt12289(int16_t a[END_N])
{
    for (unsigned stride = END_N / 4; stride >= 8; stride /= 4)
        ct_layer_pair(a, stride);
    for (unsigned i = 0; i < END_N; i += 16)
        ct_tail(a + i, i);
}

/* Short GS layers on the same group, in inverse twiddle consumption order. */
static inline void gs_tail(int16_t *a, unsigned i)
{
    const unsigned i1 = i / 2, i2 = END_N / 2 + i / 4, i4 = 3 * END_N / 4 + i / 8;
    const int16x8_t pv = vdupq_n_s16(P), pinv = vdupq_n_s16(PINV16_12289), bv = vdupq_n_s16(BARRETT_V_12289);
    int16x8_t v0 = vld1q_s16(a), v1 = vld1q_s16(a + 8), lo, hi, z, zq, s, d;
    v0 = neon_barrett(v0, bv, pv);
    v1 = neon_barrett(v1, bv, pv);
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
    v0 = neon_zip64_lo(s, d);
    v1 = neon_zip64_hi(s, d);
#if END_N == 1024
    z = vdupq_n_s16(ZETAS_INV[7 * END_N / 8 + i / 16]);
    d = neon_montmul(vsubq_s16(v0, v1), z, vmulq_s16(z, pinv), pv);
    v0 = neon_barrett(vaddq_s16(v0, v1), bv, pv);
    v1 = d;
#endif
    vst1q_s16(a, v0);
    vst1q_s16(a + 8, v1);
}

/* Inverse of a wide CT pair. Each GS sum is reduced before the next butterfly,
 * so all additions and subtractions have magnitude < 2P. The last pair folds
 * 1/N into the twiddle on the difference and into the writeback of the sum. */
static inline void gs_layer_pair(int16_t *a, unsigned stride)
{
    const int16x8_t pv = vdupq_n_s16(P), pinv = vdupq_n_s16(PINV16_12289), bv = vdupq_n_s16(BARRETT_V_12289);
    const int16x8_t scale = vdupq_n_s16(NINV_MONT), scaleq = vmulq_s16(scale, pinv);
    const int16x8_t half = vdupq_n_s16(P / 2);
    const unsigned blocks = END_N / (4 * stride);
    for (unsigned b = 0; b < blocks; b++) {
        unsigned k = END_N - 4 * blocks + 2 * b;
        int16x8_t zl = vdupq_n_s16(ZETAS_INV[k]), zlq = vmulq_s16(zl, pinv);
        int16x8_t zr = vdupq_n_s16(ZETAS_INV[k + 1]), zrq = vmulq_s16(zr, pinv);
        int16_t twiddle = ZETAS_INV[END_N - 2 * blocks + b];
        if (blocks == 1)
            twiddle = fqmul12289(twiddle, NINV_MONT);
        int16x8_t z = vdupq_n_s16(twiddle), zq = vmulq_s16(z, pinv);
        for (unsigned j = 0; j < stride; j += 8) {
            int16_t *p = a + 4 * stride * b + j;
            int16x8_t x0 = vld1q_s16(p), x1 = vld1q_s16(p + stride);
            int16x8_t x2 = vld1q_s16(p + 2 * stride), x3 = vld1q_s16(p + 3 * stride);
            int16x8_t l = neon_barrett(vaddq_s16(x0, x1), bv, pv);
            int16x8_t r = neon_barrett(vaddq_s16(x2, x3), bv, pv);
            int16x8_t dl = neon_montmul(vsubq_s16(x0, x1), zl, zlq, pv);
            int16x8_t dr = neon_montmul(vsubq_s16(x2, x3), zr, zrq, pv);
            x0 = vaddq_s16(l, r);
            x1 = vaddq_s16(dl, dr);
            x2 = neon_montmul(vsubq_s16(l, r), z, zq, pv);
            x3 = neon_montmul(vsubq_s16(dl, dr), z, zq, pv);
            if (blocks == 1) {
                x0 = neon_center_exact(neon_montmul(x0, scale, scaleq, pv), pv, half);
                x1 = neon_center_exact(neon_montmul(x1, scale, scaleq, pv), pv, half);
                x2 = neon_center_exact(x2, pv, half);
                x3 = neon_center_exact(x3, pv, half);
            } else {
                x0 = neon_barrett(x0, bv, pv);
                x1 = neon_barrett(x1, bv, pv);
            }
            vst1q_s16(p, x0);
            vst1q_s16(p + stride, x1);
            vst1q_s16(p + 2 * stride, x2);
            vst1q_s16(p + 3 * stride, x3);
        }
    }
}

void poly_intt12289(int16_t a[END_N])
{
    for (unsigned i = 0; i < END_N; i += 16)
        gs_tail(a + i, i);
    /* The short layers stop before len=8 (END-512) or len=16 (END-1024). */
    for (unsigned stride = END_N / 64; stride <= END_N / 4; stride *= 4)
        gs_layer_pair(a, stride);
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
