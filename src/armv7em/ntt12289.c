/*
 * END KEM -- ARMv7E-M kernel: complete NTT mod 12289 with packed 16-bit DSP arithmetic.
 * Centred values; Barrett every two layers (bound 6145 + 2*12289 < 2^15); the len = 1 layer works inside a word.
 */
#include "poly.h"
#include "end_arith.h"
#include "dsp.h"

#define P 12289
#define NEGP (-12289)
#define PINV (-12287)    /* 12289^-1 mod 2^16 */
#if END_N == 512
#define ZETAS     end_zetas12289_512
#define ZETAS_INV end_zetas12289_inv_512
#define NINV_MONT END_NINV_MONT12289_512
#else
#define ZETAS     end_zetas12289_1024
#define ZETAS_INV end_zetas12289_inv_1024
#define NINV_MONT END_NINV_MONT12289_1024
#endif

static inline int32_t barrett_pair(int32_t x)
{
    int32_t a = barrett12289((int16_t)dsp_lo(x)), b = barrett12289((int16_t)dsp_hi(x));
    return dsp_pkhbt(a, b);
}
static inline void reduce_all(int32_t *w)
{
    unsigned i;
    for (i = 0; i < END_N / 2; i++) w[i] = barrett_pair(w[i]);
}

void poly_ntt12289(int16_t a[END_N])
{
    int32_t *w = (int32_t *)a;
    unsigned i, len, start, j, k = 1, layer = 0;
    reduce_all(w);
    for (len = END_N >> 1; len >= 2; len >>= 1) {
        unsigned hl = len >> 1;
        for (start = 0; start < END_N / 2; start += len) {
            int32_t z = ZETAS[k++];
            int32_t zq = (int16_t)(z * PINV);
            for (j = start; j < start + hl; j++) {
                int32_t x = w[j], y = w[j + hl];
                int32_t t = dsp_montmul_pair(y, z, zq, NEGP);
                w[j + hl] = dsp_ssub16(x, t);
                w[j] = dsp_sadd16(x, t);
            }
        }
        if ((++layer & 1) == 0)
            reduce_all(w);
    }
    if (layer & 1)
        reduce_all(w);
    /* len = 1: pairs inside a word, one twiddle per word */
    for (j = 0; j < END_N / 2; j++) {
        int32_t x = w[j];
        int32_t z = ZETAS[k++];
        int32_t zq = (int16_t)(z * PINV);
        int32_t t = dsp_smultb(x, z), m = dsp_smultb(x, zq);
        int32_t r = dsp_smlabb(m, NEGP, t) >> 16;                 /* x.hi * z mod p */
        int32_t lo = (int16_t)x;
        w[j] = dsp_pkhbt(lo + r, lo - r);
    }
    for (i = 0; i < END_N / 2; i++) {
        int32_t t = barrett_pair(w[i]);
        w[i] = dsp_pkhbt(canon12289((int16_t)dsp_lo(t)), canon12289((int16_t)dsp_hi(t)));
    }
}

void poly_intt12289(int16_t a[END_N])
{
    int32_t *w = (int32_t *)a;
    unsigned i, len, b, j, idx = 0, layer = 0;
    reduce_all(w);
    /* len = 1 */
    for (j = 0; j < END_N / 2; j++) {
        int32_t x = w[j];
        int32_t z = ZETAS_INV[idx++];
        int32_t zq = (int16_t)(z * PINV);
        int32_t lo = (int16_t)x, hi = x >> 16;
        int32_t d = lo - hi, s = barrett12289((int16_t)(lo + hi));
        int32_t t = dsp_smulbb(d, z), m = dsp_smulbb(d, zq);
        int32_t r = dsp_smlabb(m, NEGP, t) >> 16;
        w[j] = dsp_pkhbt(s, r);
    }
    /* GS sums are reduced every layer: inputs of a layer are then bounded by p, differences by 2p < 2^15 */
    for (len = 2; len < END_N; len <<= 1) {
        unsigned nb = END_N / (2 * len), hl = len >> 1;
        for (b = 0; b < nb; b++) {
            int32_t z = ZETAS_INV[idx++];
            int32_t zq = (int16_t)(z * PINV);
            unsigned start = len * b;
            for (j = start; j < start + hl; j++) {
                int32_t x = w[j], y = w[j + hl];
                w[j] = barrett_pair(dsp_sadd16(x, y));
                w[j + hl] = dsp_montmul_pair(dsp_ssub16(x, y), z, zq, NEGP);
            }
        }
    }
    (void)layer;
    {
        int32_t z = NINV_MONT, zq = (int16_t)(z * PINV);
        for (i = 0; i < END_N / 2; i++) {
            int32_t t = dsp_montmul_pair(w[i], z, zq, NEGP);
            int32_t lo = barrett12289((int16_t)dsp_lo(t)), hi = barrett12289((int16_t)dsp_hi(t));
            w[i] = dsp_pkhbt(lo, hi);
        }
    }
}

void poly_pointwise12289(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N])
{
    const int32_t *wa = (const int32_t *)a, *wb = (const int32_t *)b;
    int32_t *wr = (int32_t *)r;
    unsigned i;
    const int32_t r2 = END_MONT_R2_MOD12289, r2q = (int16_t)(r2 * PINV);
    for (i = 0; i < END_N / 2; i++) {
        int32_t x = wa[i], y = wb[i];
        /* per half: t = x*y ; montgomery -> x y R^-1 ; then * R^2 montgomery -> x y */
        int32_t t0 = dsp_smulbb(x, y), t1 = dsp_smultt(x, y);
        int32_t m0 = (int16_t)(t0 * PINV), m1 = (int16_t)(t1 * PINV);
        int32_t u0 = dsp_smlabb(m0, NEGP, t0) >> 16, u1 = dsp_smlabb(m1, NEGP, t1) >> 16;
        int32_t u = dsp_pkhbt(u0, u1);
        int32_t v = dsp_montmul_pair(u, r2, r2q, NEGP);
        wr[i] = dsp_pkhbt(canon12289((int16_t)dsp_lo(v)), canon12289((int16_t)dsp_hi(v)));
    }
}
