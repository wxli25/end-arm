/*
 * END KEM -- AArch64 NEON kernel: constant-time vector helpers.
 *   poly_rot16 : negacyclic rotation by a secret k, log2(n) masked barrel-shifter stages
 *                (stages < 8 use EXT across vector boundaries, stages >= 8 move whole vectors)
 *   sqnorm, cmov: straightforward vector code;  fcl_select: scalar masked insertion (as in ref)
 */
#include <string.h>
#include "poly.h"
#include "neon_common.h"

#define NV (END_N / 8)

static inline void rot_stage_small(int16x8_t *dst, const int16x8_t *src, unsigned sh, uint16x8_t m)
{
    /* dst = m ? x^sh * src : src ; prev vector for i = 0 is -src[NV-1] (negacyclic wrap) */
    int16x8_t prev = vnegq_s16(src[NV - 1]);
    unsigned i;
    for (i = 0; i < NV; i++) {
        int16x8_t cur = src[i], sh_v;
        switch (sh) {
        case 1: sh_v = vextq_s16(prev, cur, 7); break;
        case 2: sh_v = vextq_s16(prev, cur, 6); break;
        default: sh_v = vextq_s16(prev, cur, 4); break;
        }
        dst[i] = vbslq_s16(m, sh_v, cur);
        prev = cur;
    }
}

void poly_rot16(int16_t r[END_N], const int16_t a[END_N], uint32_t k)
{
    int16x8_t t[NV], u[NV];
    unsigned b, i;
    for (i = 0; i < NV; i++) t[i] = vld1q_s16(a + 8 * i);
    for (b = 0; b < END_LOGN; b++) {
        unsigned sh = 1u << b;
        uint16x8_t m = vdupq_n_u16((uint16_t)(0u - ((k >> b) & 1u)));
        if (sh < 8) {
            rot_stage_small(u, t, sh, m);
        } else {
            unsigned sv = sh / 8;                       /* whole-vector shift */
            for (i = 0; i < NV; i++) {
                int16x8_t src = (i >= sv) ? t[i - sv] : vnegq_s16(t[i - sv + NV]);   /* public indices */
                u[i] = vbslq_s16(m, src, t[i]);
            }
        }
        memcpy(t, u, sizeof t);
    }
    for (i = 0; i < NV; i++) vst1q_s16(r + 8 * i, t[i]);
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

void poly_cmov32(int32_t r[END_N], const int32_t a[END_N], uint32_t mask)
{
    uint32x4_t m = vdupq_n_u32(mask);
    unsigned i;
    for (i = 0; i < END_N; i += 4)
        vst1q_s32(r + i, vbslq_s32(m, vld1q_s32(a + i), vld1q_s32(r + i)));
}

void fcl_select(uint32_t idx[], unsigned m0, const uint32_t scores[END_N])
{
    uint32_t best[8], bidx[8];
    unsigned i, j;
    for (j = 0; j < m0; j++) { best[j] = 0xFFFFFFFFu; bidx[j] = 0; }
    for (i = 0; i < END_N; i++) {
        uint32_t cs = scores[i], ci = i;
        for (j = 0; j < m0; j++) {
            uint32_t lt = ct_mask_lt_u32(cs, best[j]);
            uint32_t ns = (best[j] & lt) | (cs & ~lt);
            uint32_t ni = (bidx[j] & lt) | (ci & ~lt);
            best[j] = (cs & lt) | (best[j] & ~lt);
            bidx[j] = (ci & lt) | (bidx[j] & ~lt);
            cs = ns; ci = ni;
        }
    }
    for (j = 0; j < m0; j++) idx[j] = bidx[j];
}
