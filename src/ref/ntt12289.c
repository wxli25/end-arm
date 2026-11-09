/*
 * END KEM -- portable C kernel: complete NTT mod 12289 (Montgomery, R = 2^16).
 *
 * Used once per decryption for w*c1 mod 12289 (CRT part of the mod-2qQ reduction).
 * Values are kept in [0, 12289) throughout (one conditional correction per add/sub).
 */
#include "poly.h"
#include "end_arith.h"

#define P END_QP

#if END_N == 512
#define ZETAS     end_zetas12289_512
#define ZETAS_INV end_zetas12289_inv_512
#define NINV_MONT END_NINV_MONT12289_512
#else
#define ZETAS     end_zetas12289_1024
#define ZETAS_INV end_zetas12289_inv_1024
#define NINV_MONT END_NINV_MONT12289_1024
#endif

static inline int16_t csub(int16_t x)   /* [0,2P) -> [0,P) */
{
    x = (int16_t)(x - P);
    return (int16_t)(x + ((x >> 15) & P));
}
static inline int16_t cadd(int16_t x)   /* (-P,P) -> [0,P) */
{
    return (int16_t)(x + ((x >> 15) & P));
}

void poly_ntt12289(int16_t a[END_N])
{
    unsigned len, start, j, k = 1;
    for (j = 0; j < END_N; j++)
        a[j] = canon12289(a[j]);
    for (len = END_N >> 1; len >= 1; len >>= 1) {
        for (start = 0; start < END_N; start += 2 * len) {
            int16_t z = ZETAS[k++];
            for (j = start; j < start + len; j++) {
                int16_t t = cadd(fqmul12289(z, a[j + len]));      /* [0,P) */
                a[j + len] = cadd((int16_t)(a[j] - t));
                a[j] = csub((int16_t)(a[j] + t));
            }
        }
    }
}

void poly_intt12289(int16_t a[END_N])
{
    unsigned len, b, j, idx = 0;
    for (j = 0; j < END_N; j++)
        a[j] = canon12289(a[j]);
    for (len = 1; len < END_N; len <<= 1) {
        unsigned nb = END_N / (2 * len);
        for (b = 0; b < nb; b++) {
            int16_t zi = ZETAS_INV[idx++];
            unsigned start = 2 * len * b;
            for (j = start; j < start + len; j++) {
                int16_t t = a[j];
                a[j] = csub((int16_t)(t + a[j + len]));
                a[j + len] = cadd(fqmul12289(zi, cadd((int16_t)(t - a[j + len]))));
            }
        }
    }
    for (j = 0; j < END_N; j++) {
        int16_t x = fqmul12289(a[j], NINV_MONT);     /* (-P, P) */
        x = (int16_t)(x - (P & ~((x - (P / 2 + 1)) >> 15)));   /* > 6144 -> subtract P */
        x = (int16_t)(x + (P & ((x + P / 2) >> 15)));          /* < -6144 -> add P */
        a[j] = x;
    }
}

void poly_pointwise12289(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N])
{
    unsigned j;
    for (j = 0; j < END_N; j++) {
        int16_t t = fqmul12289(a[j], b[j]);                    /* a b R^-1 */
        r[j] = cadd(fqmul12289(t, (int16_t)END_MONT_R2_MOD12289));   /* a b R^-1 R^2 R^-1 = a b */
    }
}
