/*
 * END KEM -- portable C kernel: NTT mod 3329 (same 7-layer tree as mod 257; Kyber-style Montgomery).
 *
 * Used by the decoder for the exact small quantities s~' (|.| <= 129) and r*f (|.| <= 864).
 */
#include "poly.h"
#include "end_arith.h"

#define L END_LEAF
#define P END_P3

void poly_ntt3329(int16_t a[END_N])
{
    unsigned len, start, j, k = 1;
    for (j = 0; j < END_N; j++)
        a[j] = barrett3329(a[j]);                 /* centred, |a| <= 1664 */
    for (len = END_N >> 1; len >= L; len >>= 1) {
        for (start = 0; start < END_N; start += 2 * len) {
            int16_t z = end_zetas3329_mont[k++];
            for (j = start; j < start + len; j++) {
                int16_t t = fqmul3329(z, a[j + len]);          /* exact zeta*a, in (-P, P) */
                a[j + len] = (int16_t)(a[j] - t);              /* growth < P per layer: 1664 + 7P < 2^15 */
                a[j] = (int16_t)(a[j] + t);
            }
        }
    }
    for (j = 0; j < END_N; j++)
        a[j] = canon3329(a[j]);
}

void poly_intt3329(int16_t a[END_N])
{
    unsigned len, b, j;
    for (j = 0; j < END_N; j++)
        a[j] = barrett3329(a[j]);
    for (len = L; len <= END_N >> 1; len <<= 1) {
        unsigned nb = END_N / (2 * len);
        for (b = 0; b < nb; b++) {
            int16_t zi = end_zetas3329_inv_mont[nb + b];
            unsigned start = 2 * len * b;
            for (j = start; j < start + len; j++) {
                int16_t t = a[j];
                a[j] = barrett3329((int16_t)(t + a[j + len]));
                a[j + len] = fqmul3329(zi, (int16_t)(t - a[j + len]));
            }
        }
    }
    for (j = 0; j < END_N; j++)
        a[j] = barrett3329(fqmul3329(a[j], END_INV128_MONT_MOD3329));
}

static inline int16_t leaf_zeta3329_mont(unsigned blk)
{
    int16_t z = end_zetas3329_mont[64 + (blk >> 1)];
    return (blk & 1) ? (int16_t)(-z) : z;
}

void poly_basemul3329(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N])
{
    unsigned blk, i, k;
    for (blk = 0; blk < END_N / L; blk++) {
        const int16_t *x = a + blk * L, *y = b + blk * L;
        int16_t *o = r + blk * L;
        /* zR2 = mont(zeta*R * R^2) = zeta*R^2 : cancels the two R^-1 factors of mont(hi) and mont(h*zR2) */
        int16_t zR2 = fqmul3329(leaf_zeta3329_mont(blk), END_MONT_R2_MOD3329);
        int16_t tmp[L];
        for (k = 0; k < L; k++) {
            int32_t lo = 0, hi = 0;
            for (i = 0; i < L; i++) {
                if (i <= k) lo += (int32_t)x[i] * y[k - i];
                else        hi += (int32_t)x[i] * y[k + L - i];
            }
            /* hi < L * P^2 < 2^27; mont(hi) = hi/R; mont(hi/R * zeta R^2) = hi*zeta */
            int16_t h = mont3329(hi);
            int16_t yk = mont3329((int32_t)h * zR2);
            tmp[k] = red3329_32(lo + yk);
        }
        for (k = 0; k < L; k++) o[k] = tmp[k];
    }
}
