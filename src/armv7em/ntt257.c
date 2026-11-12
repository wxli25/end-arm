/*
 * END KEM -- ARMv7E-M kernel: NTT mod 257 (reference layout) with packed 16-bit DSP arithmetic.
 * Two coefficients per 32-bit word; butterflies use SMULxx/SMLAxx Montgomery products and SADD16/SSUB16.
 * Range discipline as in the portable kernel (forward < 8*257 after reduction of the inputs).
 */
#include "poly.h"
#include "end_arith.h"
#include "dsp.h"

#define L END_LEAF
#define P 257
#define NEGP (-257)
#define PINV (-255)      /* 257^-1 mod 2^16 */

static inline int32_t red257_pair(int32_t x)
{
    int32_t a = red257((int16_t)dsp_lo(x)), b = red257((int16_t)dsp_hi(x));
    return dsp_pkhbt(a, b);
}

void poly_ntt257(int16_t a[END_N])
{
    unsigned i, len, start, j, k = 1;
    for (i = 0; i < END_N / 2; i++)
        dsp_store_pair(a, i, red257_pair(dsp_load_pair(a, i)));
    for (len = END_N >> 1; len >= L; len >>= 1) {
        unsigned hl = len >> 1;
        for (start = 0; start < END_N / 2; start += len) {
            int32_t z = end_zetas257_mont[k++];
            int32_t zq = (int16_t)(z * PINV);
            for (j = start; j < start + hl; j++) {
                int32_t x = dsp_load_pair(a, j), y = dsp_load_pair(a, j + hl);
                int32_t t = dsp_montmul_pair(y, z, zq, NEGP);
                dsp_store_pair(a, j + hl, dsp_ssub16(x, t));
                dsp_store_pair(a, j, dsp_sadd16(x, t));
            }
        }
    }
    for (i = 0; i < END_N / 2; i++)
        dsp_store_pair(a, i, red257_pair(dsp_load_pair(a, i)));
}

void poly_intt257(int16_t a[END_N])
{
    unsigned i, len, b, j;
    for (i = 0; i < END_N / 2; i++)
        dsp_store_pair(a, i, red257_pair(dsp_load_pair(a, i)));
    for (len = L; len <= END_N >> 1; len <<= 1) {
        unsigned nb = END_N / (2 * len), hl = len >> 1;
        for (b = 0; b < nb; b++) {
            int32_t z = end_zetas257_inv_mont[nb + b];
            int32_t zq = (int16_t)(z * PINV);
            unsigned start = len * b;
            for (j = start; j < start + hl; j++) {
                int32_t x = dsp_load_pair(a, j), y = dsp_load_pair(a, j + hl);
                dsp_store_pair(a, j, red257_pair(dsp_sadd16(x, y)));
                dsp_store_pair(a, j + hl, dsp_montmul_pair(dsp_ssub16(x, y), z, zq, NEGP));
            }
        }
    }
    {
        int32_t z = END_INV128_MONT_MOD257, zq = (int16_t)(z * PINV);
        for (i = 0; i < END_N / 2; i++)
            dsp_store_pair(a, i, red257_pair(dsp_montmul_pair(dsp_load_pair(a, i), z, zq, NEGP)));
    }
}

static inline int16_t leaf_zeta257(unsigned blk)
{
    int16_t z = end_zetas257[64 + (blk >> 1)];
    return (blk & 1) ? (int16_t)(END_Q - z) : z;
}

#if L == 4
/* leaf product mod (x^4 - z): dual 16-bit MACs on packed pairs A01 = [a1:a0], A23 = [a3:a2] */
void poly_basemul257(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N])
{
    unsigned blk;
    for (blk = 0; blk < END_N / 4; blk++) {
        int32_t A01 = dsp_load_pair(a, 2 * blk), A23 = dsp_load_pair(a, 2 * blk + 1), B01 = dsp_load_pair(b, 2 * blk), B23 = dsp_load_pair(b, 2 * blk + 1);
        int32_t z = leaf_zeta257(blk);
        int32_t h0 = dsp_smultt(A01, B23) + dsp_smultt(A23, B01) + dsp_smulbb(A23, B23);   /* a1b3 + a3b1 + a2b2 */
        int32_t r0 = dsp_smulbb(A01, B01) + z * h0;
        int32_t r1 = dsp_smuadx(A01, B01) + z * dsp_smuadx(A23, B23);                     /* a0b1+a1b0 + z(a2b3+a3b2) */
        int32_t r2 = dsp_smulbb(A01, B23) + dsp_smultt(A01, B01) + dsp_smulbb(A23, B01) + z * dsp_smultt(A23, B23);
        int32_t r3 = dsp_smuadx(A01, B23) + dsp_smuadx(A23, B01);
        dsp_store_pair(r, 2 * blk, dsp_pkhbt(red257_32(r0), red257_32(r1)));
        dsp_store_pair(r, 2 * blk + 1, dsp_pkhbt(red257_32(r2), red257_32(r3)));
    }
}
#else
void poly_basemul257(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N])
{
    unsigned blk, i, k;
    for (blk = 0; blk < END_N / L; blk++) {
        const int16_t *x = a + blk * L, *y = b + blk * L;
        int32_t z = leaf_zeta257(blk);
        int16_t tmp[L];
        for (k = 0; k < L; k++) {
            int32_t lo = 0, hi = 0;
            for (i = 0; i < L; i++) {
                if (i <= k) lo += (int32_t)x[i] * y[k - i];
                else        hi += (int32_t)x[i] * y[k + L - i];
            }
            tmp[k] = red257_32(lo + z * hi);
        }
        for (k = 0; k < L; k++) r[blk * L + k] = tmp[k];
    }
}
#endif
