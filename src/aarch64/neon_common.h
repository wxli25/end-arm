/*
 * END KEM -- AArch64 NEON helpers shared by the vector kernels.
 *
 * Montgomery multiplication (R = 2^16) in four instructions (Becker et al., "Neon NTT"):
 *   m  = a * zq            (zq = z * p^-1 mod 2^16, low half)
 *   h  = (2 a z) >> 16     (SQDMULH)
 *   l2 = (2 m p) >> 16     (SQDMULH)
 *   r  = (h - l2) / 2      (SHSUB)   = (a z - m p) / 2^16  exactly, in (-p, p)
 * z is the twiddle in Montgomery form (z R mod p, centred), so r = a * z_plain mod p.
 */
#ifndef END_NEON_COMMON_H
#define END_NEON_COMMON_H

#include <arm_neon.h>
#include "end_params.h"
#include "end_arith.h"

#define PINV16_257   (-255)      /* 257^-1   mod 2^16 */
#define PINV16_3329  (-3327)     /* 3329^-1  mod 2^16 */
#define PINV16_12289 (-12287)    /* 12289^-1 mod 2^16 */
#define BARRETT_V_3329  20159    /* round(2^26 / 3329)  */
#define BARRETT_V_12289 5461     /* round(2^26 / 12289) */

static inline int16x8_t neon_montmul(int16x8_t a, int16x8_t z, int16x8_t zq, int16x8_t pv)
{
    int16x8_t m = vmulq_s16(a, zq);
    int16x8_t h = vqdmulhq_s16(a, z);
    int16x8_t l2 = vqdmulhq_s16(m, pv);
    return vhsubq_s16(h, l2);
}

/* any int16 -> [0,257) (Fermat: x = hi*256 + lo = lo - hi mod 257) */
static inline int16x8_t neon_red257(int16x8_t a)
{
    const int16x8_t q = vdupq_n_s16(257), m255 = vdupq_n_s16(255);
    int16x8_t r = vsubq_s16(vandq_s16(a, m255), vshrq_n_s16(a, 8));        /* [-127, 383] */
    r = vsubq_s16(r, vandq_s16(q, vreinterpretq_s16_u16(vcgeq_s16(r, q)))); /* >= 257 -> -257 */
    r = vaddq_s16(r, vandq_s16(q, vshrq_n_s16(r, 15)));                     /* < 0 -> +257 */
    return r;
}

/* centred Barrett reduction, |r| <= p/2 + 1 */
static inline int16x8_t neon_barrett(int16x8_t a, int16x8_t v, int16x8_t pv)
{
    int16x8_t t = vqrdmulhq_s16(a, v);
    t = vrshrq_n_s16(t, 11);
    return vmlsq_s16(a, t, pv);
}

/* exact centring to [-(p-1)/2, (p-1)/2]: fixes the +-1 slack of neon_barrett at the boundary */
static inline int16x8_t neon_center_exact(int16x8_t a, int16x8_t pv, int16x8_t half)
{
    /* a > half -> a - p ; a < -half -> a + p */
    a = vsubq_s16(a, vandq_s16(pv, vreinterpretq_s16_u16(vcgtq_s16(a, half))));
    a = vaddq_s16(a, vandq_s16(pv, vreinterpretq_s16_u16(vcltq_s16(a, vnegq_s16(half)))));
    return a;
}

/* centred -> canonical [0,p) */
static inline int16x8_t neon_canon(int16x8_t a, int16x8_t pv)
{
    return vaddq_s16(a, vandq_s16(pv, vshrq_n_s16(a, 15)));
}

/* int32x4 -> int16x4 : x * R^-1 mod p, in (-p, p)  (|x| < 2^15 p) */
static inline int16x4_t neon_montred32(int32x4_t x, int32x4_t pinv32, int16x4_t p16)
{
    int16x4_t m = vmovn_s32(vmulq_s32(x, pinv32));          /* x * p^-1 mod 2^16 */
    int32x4_t t = vmlsl_s16(x, m, p16);                     /* x - m p  (= 0 mod 2^16) */
    return vshrn_n_s32(t, 16);
}

/* 64-bit lane interleave: [a.lo, b.lo] and [a.hi, b.hi] */
static inline int16x8_t neon_zip64_lo(int16x8_t a, int16x8_t b)
{
    return vreinterpretq_s16_s64(vzip1q_s64(vreinterpretq_s64_s16(a), vreinterpretq_s64_s16(b)));
}
static inline int16x8_t neon_zip64_hi(int16x8_t a, int16x8_t b)
{
    return vreinterpretq_s16_s64(vzip2q_s64(vreinterpretq_s64_s16(a), vreinterpretq_s64_s16(b)));
}
/* 32-bit lane unzip / zip */
static inline int16x8_t neon_uzp32_1(int16x8_t a, int16x8_t b)
{
    return vreinterpretq_s16_s32(vuzp1q_s32(vreinterpretq_s32_s16(a), vreinterpretq_s32_s16(b)));
}
static inline int16x8_t neon_uzp32_2(int16x8_t a, int16x8_t b)
{
    return vreinterpretq_s16_s32(vuzp2q_s32(vreinterpretq_s32_s16(a), vreinterpretq_s32_s16(b)));
}
static inline int16x8_t neon_zip32_1(int16x8_t a, int16x8_t b)
{
    return vreinterpretq_s16_s32(vzip1q_s32(vreinterpretq_s32_s16(a), vreinterpretq_s32_s16(b)));
}
static inline int16x8_t neon_zip32_2(int16x8_t a, int16x8_t b)
{
    return vreinterpretq_s16_s32(vzip2q_s32(vreinterpretq_s32_s16(a), vreinterpretq_s32_s16(b)));
}

#endif /* END_NEON_COMMON_H */
