/*
 * END KEM -- AArch64 NEON kernel: constant-time vector helpers and the inner loops of the decoder.
 *   poly_rot16 : negacyclic rotation by a secret k. The part of k below 8 is one table lookup per vector
 *                (TBL over the previous and the current vector, indices derived from k: TBL has
 *                data-independent timing); the rest is a barrel shifter over whole vectors.
 *   fcl_select : the m0 smallest (score, index) keys with vector min/max, four keys per step
 *   sqnorm, cmov and the decoder loops: straightforward vector code
 */
#include <string.h>
#include "poly.h"
#include "neon_common.h"

#define NV (END_N / 8)
#define ROT_STAGES (END_LOGN - 3)       /* barrel-shifter stages over whole vectors */

void poly_rot16(int16_t r[END_N], const int16_t a[END_N], uint32_t k)
{
    int16_t tmp[END_N] END_ALIGN;
    /* the buffers alternate between r and tmp so that the last stage writes r */
    int16_t *dst = (ROT_STAGES & 1) ? tmp : r, *src;
    static const uint8_t base[16] = { 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31 };
    unsigned b, i;
    /* shift by s = k mod 8 coefficients: output byte j of a vector is byte 16 + j - 2s of (previous | current) */
    {
        const uint8x16_t sel = vsubq_u8(vld1q_u8(base), vdupq_n_u8((uint8_t)(2 * (k & 7u))));
        uint8x16x2_t tab;
        tab.val[0] = vreinterpretq_u8_s16(vnegq_s16(vld1q_s16(a + END_N - 8)));       /* the vector that wraps around */
        for (i = 0; i < NV; i++) {
            tab.val[1] = vreinterpretq_u8_s16(vld1q_s16(a + 8 * i));
            vst1q_s16(dst + 8 * i, vreinterpretq_s16_u8(vqtbl2q_u8(tab, sel)));
            tab.val[0] = tab.val[1];
        }
    }
    /* shift by k / 8 vectors, one masked stage per bit; the indices below are public */
    for (b = 0; b < ROT_STAGES; b++) {
        const unsigned sv = 1u << b;
        const uint16x8_t m = vdupq_n_u16((uint16_t)(0u - ((k >> (b + 3)) & 1u)));
        src = dst;
        dst = (src == r) ? tmp : r;
        for (i = 0; i < sv; i++)
            vst1q_s16(dst + 8 * i, vbslq_s16(m, vnegq_s16(vld1q_s16(src + 8 * (i + NV - sv))), vld1q_s16(src + 8 * i)));
        for (i = sv; i < NV; i++)
            vst1q_s16(dst + 8 * i, vbslq_s16(m, vld1q_s16(src + 8 * (i - sv)), vld1q_s16(src + 8 * i)));
    }
}

/* The pair (score, index) is one 32-bit key (score in the high bits), as in the portable version. Every lane keeps
 * its own sorted list of the m0 smallest keys; the m0 smallest keys overall are among these 4 m0 values. */
void fcl_select(uint32_t idx[], unsigned m0, const uint32_t scores[END_N])
{
    static const uint32_t lane[4] = { 0, 1, 2, 3 };
    uint32x4_t best[END_FCL_MAX], index = vld1q_u32(lane);
    uint32_t all[4 * END_FCL_MAX], out[END_FCL_MAX];
    unsigned i, j;
    for (j = 0; j < m0; j++) { best[j] = vdupq_n_u32(0xFFFFFFFFu); out[j] = 0xFFFFFFFFu; }
    for (i = 0; i < END_N; i += 4) {
        uint32x4_t cs = vorrq_u32(vshlq_n_u32(vld1q_u32(scores + i), END_LOGN), index);
        for (j = 0; j < m0; j++) {
            uint32x4_t lo = vminq_u32(best[j], cs);
            cs = vmaxq_u32(best[j], cs);
            best[j] = lo;
        }
        index = vaddq_u32(index, vdupq_n_u32(4));
    }
    for (j = 0; j < m0; j++) vst1q_u32(all + 4 * j, best[j]);
    for (i = 0; i < 4 * m0; i++) {
        uint32_t cs = all[i];
        for (j = 0; j < m0; j++) {
            uint32_t d = (cs - out[j]) & (uint32_t)(((uint64_t)cs - out[j]) >> 32);
            out[j] += d;
            cs -= d;
        }
    }
    for (j = 0; j < m0; j++) idx[j] = out[j] & (END_N - 1);
}

uint64_t poly_sqnorm16(const int16_t a[END_N])
{
    int64x2_t acc = vdupq_n_s64(0);
    unsigned i;
    for (i = 0; i < END_N; i += 8) {
        int16x8_t v = vld1q_s16(a + i);
        /* squares of int16 reach 2^30: widen to 64 bits after every 4 products (exact for any input) */
        acc = vaddq_s64(acc, vpaddlq_s32(vmull_s16(vget_low_s16(v), vget_low_s16(v))));
        acc = vaddq_s64(acc, vpaddlq_s32(vmull_s16(vget_high_s16(v), vget_high_s16(v))));
    }
    return (uint64_t)vaddvq_s64(acc);
}

void poly_cmov16(int16_t r[END_N], const int16_t a[END_N], uint16_t mask)
{
    uint16x8_t m = vdupq_n_u16(mask);
    unsigned i;
    for (i = 0; i < END_N; i += 8)
        vst1q_s16(r + i, vbslq_s16(m, vld1q_s16(a + i), vld1q_s16(r + i)));
}

/* ------------------------------------------------------------------------------------------------------------
 * Inner loops of the trapdoor decoder (poly.h). They compute exactly what the portable versions in
 * src/common/end_deckern.c compute; end_backend.h tells that file to leave these out.
 * ------------------------------------------------------------------------------------------------------------ */

/* x in [-256, 512] -> [0, 257): one conditional addition and one conditional subtraction */
static inline int16x8_t canon257_small(int16x8_t x, int16x8_t q)
{
    x = vaddq_s16(x, vandq_s16(q, vshrq_n_s16(x, 15)));                               /* [0, 513) */
    return vsubq_s16(x, vandq_s16(q, vreinterpretq_s16_u16(vcgeq_s16(x, q))));        /* [0, 257) */
}

void poly_lift257(int16_t out[END_N], const int16_t t[END_N], const int16_t cst[END_N], uint16_t parity)
{
    const int16x8_t q = vdupq_n_s16(END_Q), q2 = vdupq_n_s16(2 * END_Q), one = vdupq_n_s16(1);
    const int16x8_t par = vdupq_n_s16((int16_t)(parity & 1));
    unsigned i;
    for (i = 0; i < END_N; i += 8) {
        int16x8_t x = vsubq_s16(vshlq_n_s16(vld1q_s16(t + i), 1), vld1q_s16(cst + i));   /* 2t - cst in [-256, 512] */
        x = canon257_small(canon257_small(x, q), q);                                       /* (2t - cst) mod 257 */
        /* add 257 where the parity is wrong, then centre: [0, 514) -> [-257, 257) */
        x = vaddq_s16(x, vandq_s16(q, vreinterpretq_s16_u16(vtstq_s16(veorq_s16(x, par), one))));
        x = vsubq_s16(x, vandq_s16(q2, vreinterpretq_s16_u16(vcgeq_s16(x, q))));
        vst1q_s16(out + i, x);
    }
}

void poly_cand_round(int16_t r[END_N], int32_t c2t[END_N], const int32_t C2[END_N],
                     const int16_t w0[END_N], const int16_t w1[END_N], int32_t s0, int32_t s1)
{
    /* the arithmetic of round_div_2qqp_rem (end_arith.h) on four 32-bit lanes; the quotient estimate uses
     * SQDMULH with M = floor(2^31 / 2qQ) = 339: floor(x M / 2^31) is floor(x / 2qQ) or one less for x < 2^29 */
    const int16x4_t m0 = vdup_n_s16((int16_t)(-2 * END_Q * s0)), m1 = vdup_n_s16((int16_t)(-2 * END_Q * s1));
    const int32x4_t off = vdupq_n_s32(END_QQP + 32 * END_2QQP), k = vdupq_n_s32(END_2QQP);
    const int32x4_t m = vdupq_n_s32(339), half = vdupq_n_s32(END_QQP), c31 = vdupq_n_s32(31);
    unsigned i;
    for (i = 0; i < END_N; i += 4) {
        int32x4_t x = vld1q_s32(C2 + i), q, rem, neg;
        x = vmlal_s16(x, vld1_s16(w0 + i), m0);
        x = vmlal_s16(x, vld1_s16(w1 + i), m1);                  /* v = C2 - 2q (s0 w0 + s1 w1) */
        x = vaddq_s32(x, off);
        q = vqdmulhq_s32(x, m);
        rem = vsubq_s32(vmlsq_s32(x, q, k), k);                  /* [-2qQ, 2qQ) */
        neg = vshrq_n_s32(rem, 31);
        q = vsubq_s32(vaddq_s32(q, neg), c31);
        rem = vsubq_s32(vaddq_s32(rem, vandq_s32(k, neg)), half);
        vst1_s16(r + i, vmovn_s32(q));
        vst1q_s32(c2t + i, rem);
    }
}

uint64_t poly_cand_combine(int16_t t[END_N], const int16_t sp[END_N], const int16_t F0[END_N], const int16_t F1[END_N],
                           const int16_t rf[END_N], int32_t s0, int32_t s1)
{
    const int16x8_t c0 = vdupq_n_s16((int16_t)(2 * s0)), c1 = vdupq_n_s16((int16_t)(2 * s1));
    int64x2_t acc = vdupq_n_s64(0);
    unsigned i;
    for (i = 0; i < END_N; i += 8) {
        int16x8_t x = vmulq_s16(vld1q_s16(F0 + i), c0);
        x = vmlaq_s16(x, vld1q_s16(F1 + i), c1);
        x = vaddq_s16(x, vshlq_n_s16(vld1q_s16(rf + i), 1));
        x = vaddq_s16(x, vld1q_s16(sp + i));
        vst1q_s16(t + i, x);
        acc = vaddq_s64(acc, vpaddlq_s32(vmull_s16(vget_low_s16(x), vget_low_s16(x))));
        acc = vaddq_s64(acc, vpaddlq_s32(vmull_s16(vget_high_s16(x), vget_high_s16(x))));
    }
    return (uint64_t)vaddvq_s64(acc);
}

uint32_t poly_cand_check(int16_t t[END_N], const int16_t b[END_N], const int16_t f0[END_N], const int16_t f1[END_N],
                         int32_t s0, int32_t s1)
{
    const int16x8_t c0 = vdupq_n_s16((int16_t)(2 * s0)), c1 = vdupq_n_s16((int16_t)(2 * s1));
    const int16x8_t one = vdupq_n_s16(1), two = vdupq_n_s16(2);
    int16x8_t acc = vdupq_n_s16(0);
    unsigned i;
    for (i = 0; i < END_N; i += 8) {
        int16x8_t x = vmulq_s16(vld1q_s16(f0 + i), c0);
        x = vmlaq_s16(x, vld1q_s16(f1 + i), c1);
        x = vsubq_s16(vld1q_s16(b + i), x);
        vst1q_s16(t + i, x);
        acc = vorrq_s16(acc, vbicq_s16(vaddq_s16(x, one), two));     /* x = +-1  <=>  (x + 1) & ~2 == 0 */
    }
    return (uint32_t)vmaxvq_u16(vreinterpretq_u16_s16(acc));
}

void poly_cmov32(int32_t r[END_N], const int32_t a[END_N], uint32_t mask)
{
    const uint32x4_t m = vdupq_n_u32(mask);
    unsigned i;
    for (i = 0; i < END_N; i += 4)
        vst1q_s32(r + i, vbslq_s32(m, vld1q_s32(a + i), vld1q_s32(r + i)));
}
