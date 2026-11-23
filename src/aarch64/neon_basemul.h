/* END KEM -- AArch64 products in the NTT leaves modulo X^L - z.
 * Inputs are canonical. Montgomery multiplication by the twiddle keeps each
 * wrapped coefficient of b in (-p,p), so every output is one widening dot product. */
#ifndef END_NEON_BASEMUL_H
#define END_NEON_BASEMUL_H

#include "neon_common.h"

/* |sum| <= L(p-1)^2 <= 88604672 for both supported primes. With reciprocal = round(2^31/p),
 * SQRDMULH estimates sum/p to within 1/2 + 8(p-1)^2/2^32 (< 1).
 * Hence sum - quotient*p is in (-p,p) and fits in int16. Narrowing the
 * operands first gives the same low 16 bits of this exact remainder. */
static inline int16x8_t leaf_reduce(int32x4_t lo, int32x4_t hi, int16x8_t p, int32x4_t reciprocal)
{
    int16x8_t quotient = vcombine_s16(vmovn_s32(vqrdmulhq_s32(lo, reciprocal)),
                                      vmovn_s32(vqrdmulhq_s32(hi, reciprocal)));
    int16x8_t sum = vcombine_s16(vmovn_s32(lo), vmovn_s32(hi));
    return neon_canon(vmlsq_s16(sum, quotient, p), p);
}

#if END_LEAF == 4
/* Four independent products per 32-bit accumulator vector. */
static inline int16x8_t leaf_dot4(int16x8x4_t a, int16x8_t b0, int16x8_t b1,
                                  int16x8_t b2, int16x8_t b3, int16x8_t p, int32x4_t reciprocal)
{
    int32x4_t lo = vmull_s16(vget_low_s16(a.val[0]), vget_low_s16(b0));
    int32x4_t hi = vmull_high_s16(a.val[0], b0);
    lo = vmlal_s16(lo, vget_low_s16(a.val[1]), vget_low_s16(b1));
    hi = vmlal_high_s16(hi, a.val[1], b1);
    lo = vmlal_s16(lo, vget_low_s16(a.val[2]), vget_low_s16(b2));
    hi = vmlal_high_s16(hi, a.val[2], b2);
    lo = vmlal_s16(lo, vget_low_s16(a.val[3]), vget_low_s16(b3));
    hi = vmlal_high_s16(hi, a.val[3], b3);
    return leaf_reduce(lo, hi, p, reciprocal);
}
#else
/* One vector holds a whole eight-coefficient leaf. For row j the high j
 * coefficients of wrapped_b precede the low 8-j coefficients of b. */
static inline int16x8_t leaf_dot8(int16x8_t a, int16x8_t b, int16x8_t wrapped_b,
                                  int16x8_t p, int32x4_t reciprocal)
{
    int32x4_t lo = vmull_laneq_s16(vget_low_s16(b), a, 0);
    int32x4_t hi = vmull_high_laneq_s16(b, a, 0);
#define ADD_ROTATED(lane) do {                                               \
        int16x8_t row = vextq_s16(wrapped_b, b, 8 - (lane));                 \
        lo = vmlal_laneq_s16(lo, vget_low_s16(row), a, (lane));              \
        hi = vmlal_high_laneq_s16(hi, row, a, (lane));                       \
    } while (0)
    ADD_ROTATED(1);
    ADD_ROTATED(2);
    ADD_ROTATED(3);
    ADD_ROTATED(4);
    ADD_ROTATED(5);
    ADD_ROTATED(6);
    ADD_ROTATED(7);
#undef ADD_ROTATED
    return leaf_reduce(lo, hi, p, reciprocal);
}
#endif

static inline void neon_basemul(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N],
                                const int16_t *twiddles, int16_t modulus, int16_t pinv16)
{
    const int16x8_t p = vdupq_n_s16(modulus), pinv = vdupq_n_s16(pinv16);
    const int32x4_t reciprocal = vdupq_n_s32((int32_t)(((UINT32_C(1) << 31) + modulus / 2) / modulus));
    unsigned i;

#if END_LEAF == 4
    for (i = 0; i < END_N; i += 32) {
        int16x8x4_t av = vld4q_s16(a + i), bv = vld4q_s16(b + i), out;
        int16x4_t z4 = vld1_s16(twiddles + 64 + i / 8);
        int16x8_t z = vcombine_s16(z4, z4);
        z = vzip1q_s16(z, vnegq_s16(z));
        int16x8_t zq = vmulq_s16(z, pinv);
        int16x8_t w1 = neon_montmul(bv.val[1], z, zq, p);
        int16x8_t w2 = neon_montmul(bv.val[2], z, zq, p);
        int16x8_t w3 = neon_montmul(bv.val[3], z, zq, p);
        out.val[0] = leaf_dot4(av, bv.val[0], w3, w2, w1, p, reciprocal);
        out.val[1] = leaf_dot4(av, bv.val[1], bv.val[0], w3, w2, p, reciprocal);
        out.val[2] = leaf_dot4(av, bv.val[2], bv.val[1], bv.val[0], w3, p, reciprocal);
        out.val[3] = leaf_dot4(av, bv.val[3], bv.val[2], bv.val[1], bv.val[0], p, reciprocal);
        vst4q_s16(r + i, out);
    }
#else
    /* Adjacent leaves have twiddles z and -z. */
    for (i = 0; i < END_N; i += 16) {
        int16x8_t z = vdupq_n_s16(twiddles[64 + i / 16]);
        int16x8_t zq = vmulq_s16(z, pinv);
        int16x8_t a0 = vld1q_s16(a + i), a1 = vld1q_s16(a + i + 8);
        int16x8_t b0 = vld1q_s16(b + i), b1 = vld1q_s16(b + i + 8);
        int16x8_t w0 = neon_montmul(b0, z, zq, p);
        int16x8_t w1 = vnegq_s16(neon_montmul(b1, z, zq, p));
        vst1q_s16(r + i, leaf_dot8(a0, b0, w0, p, reciprocal));
        vst1q_s16(r + i + 8, leaf_dot8(a1, b1, w1, p, reciprocal));
    }
#endif
}

#endif /* END_NEON_BASEMUL_H */
