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

/* any int16 -> [-127, 383], congruent mod 257 (one Fermat fold, no correction) */
static inline int16x8_t neon_fold257(int16x8_t a)
{
    return vsubq_s16(vandq_s16(a, vdupq_n_s16(255)), vshrq_n_s16(a, 8));
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

/*
 * Three merged layers on eight vectors v[0..7] that lie the same distance apart (neon_ntt7.h).
 * Cooley-Tukey: twiddles z[0] | z[1], z[2] | z[3..6];  Gentleman-Sande: z[0..3] | z[4], z[5] | z[6].
 */
#define NEON_CT(i, j, k) do { \
        int16x8_t t_ = neon_montmul(v[j], z[k], zq[k], pv); \
        v[j] = vsubq_s16(v[i], t_); v[i] = vaddq_s16(v[i], t_); \
    } while (0)
#define NEON_GS(i, j, k) do { \
        int16x8_t d_ = vsubq_s16(v[i], v[j]); \
        v[i] = vaddq_s16(v[i], v[j]); v[j] = neon_montmul(d_, z[k], zq[k], pv); \
    } while (0)
static inline void neon_ct3(int16x8_t v[8], const int16x8_t z[7], const int16x8_t zq[7], int16x8_t pv)
{
    NEON_CT(0, 4, 0); NEON_CT(1, 5, 0); NEON_CT(2, 6, 0); NEON_CT(3, 7, 0);
    NEON_CT(0, 2, 1); NEON_CT(1, 3, 1); NEON_CT(4, 6, 2); NEON_CT(5, 7, 2);
    NEON_CT(0, 1, 3); NEON_CT(2, 3, 4); NEON_CT(4, 5, 5); NEON_CT(6, 7, 6);
}
static inline void neon_gs3(int16x8_t v[8], const int16x8_t z[7], const int16x8_t zq[7], int16x8_t pv)
{
    NEON_GS(0, 1, 0); NEON_GS(2, 3, 1); NEON_GS(4, 5, 2); NEON_GS(6, 7, 3);
    NEON_GS(0, 2, 4); NEON_GS(1, 3, 4); NEON_GS(4, 6, 5); NEON_GS(5, 7, 5);
    NEON_GS(0, 4, 6); NEON_GS(1, 5, 6); NEON_GS(2, 6, 6); NEON_GS(3, 7, 6);
}
#define NTT7_CAT2(a, b) a##b
#define NTT7_CAT(a, b) NTT7_CAT2(a, b)

#endif /* END_NEON_COMMON_H */
