/*
 * END KEM -- portable C kernel: NTT mod 257 (reference layout).
 *
 * 257 is a Fermat prime: 2^8 = -1 and 2^16 = 1 (mod 257). Consequences used here:
 *   * Montgomery reduction with R = 2^16 returns a*R^-1 = a (mod 257): no Montgomery factors anywhere;
 *   * the canonicalisation red257 is a shift/subtract (x = hi*256 + lo = lo - hi mod 257).
 * Range discipline (int16, no overflow):
 *   forward : inputs are first reduced to [0,257); every CT layer adds at most 257 in magnitude,
 *             so after 7 layers |a| < 8*257 = 2056.
 *   inverse : the GS sum is re-reduced every layer, the GS difference goes through mont257.
 * NOTE: int16 wrap-around is NOT harmless mod 257 (65536 = 1, not 0); the bounds above are real bounds.
 */
#include "poly.h"
#include "end_arith.h"

#define L END_LEAF

void poly_ntt257(int16_t a[END_N])
{
    unsigned len, start, j, k = 1;
    for (j = 0; j < END_N; j++)
        a[j] = red257(a[j]);                      /* [0,257) */
    for (len = END_N >> 1; len >= L; len >>= 1) {
        for (start = 0; start < END_N; start += 2 * len) {
            int16_t z = end_zetas257[k++];
            for (j = start; j < start + len; j++) {
                int16_t t = mont257((int32_t)z * a[j + len]);   /* (-257, 257) */
                a[j + len] = (int16_t)(a[j] - t);               /* |.| < 2056 after 7 layers */
                a[j] = (int16_t)(a[j] + t);
            }
        }
    }
    for (j = 0; j < END_N; j++)
        a[j] = red257(a[j]);
}

void poly_intt257(int16_t a[END_N])
{
    unsigned len, b, j;
    for (j = 0; j < END_N; j++)
        a[j] = red257(a[j]);
    for (len = L; len <= END_N >> 1; len <<= 1) {
        unsigned nb = END_N / (2 * len);
        for (b = 0; b < nb; b++) {
            int16_t zi = end_zetas257_inv[nb + b];
            unsigned start = 2 * len * b;
            for (j = start; j < start + len; j++) {
                int16_t t = a[j];
                a[j] = red257((int16_t)(t + a[j + len]));       /* keep the sum in [0,257) */
                a[j + len] = mont257((int32_t)zi * (t - a[j + len]));
            }
        }
    }
    for (j = 0; j < END_N; j++)
        a[j] = red257(mont257((int32_t)a[j] * END_INV128_MOD257));
}

static inline int16_t leaf_zeta257(unsigned blk)
{
    int16_t z = end_zetas257[64 + (blk >> 1)];
    return (blk & 1) ? (int16_t)(END_Q - z) : z;
}

/* o = x * y mod (X^L - z), inputs in [0,257), output in [0,257) */
static void leaf_mul257(int16_t *o, const int16_t *x, const int16_t *y, int16_t z)
{
    unsigned i, k;
    int16_t tmp[L];
    for (k = 0; k < L; k++) {
        int32_t lo = 0, hi = 0;
        for (i = 0; i < L; i++) {
            if (i <= k) lo += (int32_t)x[i] * y[k - i];
            else        hi += (int32_t)x[i] * y[k + L - i];
        }
        tmp[k] = red257_32(lo + (int32_t)z * hi);   /* |lo + z*hi| < 2^20 + 2^28 */
    }
    for (k = 0; k < L; k++) o[k] = tmp[k];
}

void poly_basemul257(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N])
{
    unsigned blk;
    for (blk = 0; blk < END_N / L; blk++)
        leaf_mul257(r + blk * L, a + blk * L, b + blk * L, leaf_zeta257(blk));
}

