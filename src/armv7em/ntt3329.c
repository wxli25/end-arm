/*
 * END KEM -- ARMv7E-M kernel: NTT mod 3329 (7-layer tree) with packed 16-bit DSP arithmetic.
 * Forward: inputs Barrett-reduced (|a| <= 1665), 7 lazy layers (< 1665 + 7*3329), canonical output.
 * Inverse: GS sums Barrett-reduced after every third layer (bound 2^k * 3329 < 2^15), final scaling.
 */
#include "poly.h"
#include "end_arith.h"
#include "dsp.h"

#define L END_LEAF
#define P 3329
#define NEGP (-3329)
#define PINV (-3327)     /* 3329^-1 mod 2^16 */

static inline int32_t barrett_pair(int32_t x)
{
    int32_t a = barrett3329((int16_t)dsp_lo(x)), b = barrett3329((int16_t)dsp_hi(x));
    return dsp_pkhbt(a, b);
}
static inline int32_t canon_pair(int32_t x)
{
    int32_t a = canon3329((int16_t)dsp_lo(x)), b = canon3329((int16_t)dsp_hi(x));
    return dsp_pkhbt(a, b);
}

void poly_ntt3329(int16_t a[END_N])
{
    int32_t *w = (int32_t *)a;
    unsigned i, len, start, j, k = 1;
    for (i = 0; i < END_N / 2; i++)
        w[i] = barrett_pair(w[i]);
    for (len = END_N >> 1; len >= L; len >>= 1) {
        unsigned hl = len >> 1;
        for (start = 0; start < END_N / 2; start += len) {
            int32_t z = end_zetas3329_mont[k++];
            int32_t zq = (int16_t)(z * PINV);
            for (j = start; j < start + hl; j++) {
                int32_t x = w[j], y = w[j + hl];
                int32_t t = dsp_montmul_pair(y, z, zq, NEGP);
                w[j + hl] = dsp_ssub16(x, t);
                w[j] = dsp_sadd16(x, t);
            }
        }
    }
    for (i = 0; i < END_N / 2; i++)
        w[i] = canon_pair(w[i]);
}

void poly_intt3329(int16_t a[END_N])
{
    int32_t *w = (int32_t *)a;
    unsigned i, len, b, j, layer = 0;
    for (i = 0; i < END_N / 2; i++)
        w[i] = barrett_pair(w[i]);
    for (len = L; len <= END_N >> 1; len <<= 1) {
        unsigned nb = END_N / (2 * len), hl = len >> 1;
        int reduce = ((++layer) % 3) == 0;
        for (b = 0; b < nb; b++) {
            int32_t z = end_zetas3329_inv_mont[nb + b];
            int32_t zq = (int16_t)(z * PINV);
            unsigned start = len * b;
            for (j = start; j < start + hl; j++) {
                int32_t x = w[j], y = w[j + hl];
                int32_t s = dsp_sadd16(x, y);
                w[j] = reduce ? barrett_pair(s) : s;
                w[j + hl] = dsp_montmul_pair(dsp_ssub16(x, y), z, zq, NEGP);
            }
        }
    }
    {
        int32_t z = END_INV128_MONT_MOD3329, zq = (int16_t)(z * PINV);
        for (i = 0; i < END_N / 2; i++) {
            int32_t t = dsp_montmul_pair(w[i], z, zq, NEGP);     /* (-p, p) */
            int32_t lo = barrett3329((int16_t)dsp_lo(t)), hi = barrett3329((int16_t)dsp_hi(t));
            w[i] = dsp_pkhbt(lo, hi);
        }
    }
}

static inline int16_t leaf_zeta3329_mont(unsigned blk)
{
    int16_t z = end_zetas3329_mont[64 + (blk >> 1)];
    return (blk & 1) ? (int16_t)(-z) : z;
}

/* same algorithm as the portable kernel (the scalar part is not a hot spot on Cortex-M) */
void poly_basemul3329(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N])
{
    unsigned blk, i, k;
    for (blk = 0; blk < END_N / L; blk++) {
        const int16_t *x = a + blk * L, *y = b + blk * L;
        int16_t *o = r + blk * L;
        int16_t zR2 = fqmul3329(leaf_zeta3329_mont(blk), END_MONT_R2_MOD3329);
        int16_t tmp[L];
        for (k = 0; k < L; k++) {
            int32_t lo = 0, hi = 0;
            for (i = 0; i < L; i++) {
                if (i <= k) lo += (int32_t)x[i] * y[k - i];
                else        hi += (int32_t)x[i] * y[k + L - i];
            }
            {
                int16_t h = mont3329(hi);
                int16_t yk = mont3329((int32_t)h * zR2);
                tmp[k] = red3329_32(lo + yk);
            }
        }
        for (k = 0; k < L; k++) o[k] = tmp[k];
    }
}
