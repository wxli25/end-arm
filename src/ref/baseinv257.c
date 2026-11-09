/*
 * END KEM -- portable kernel: inverse in the NTT-257 domain (key generation only).
 * Norm-descent inversion of each leaf Z_257[x]/(x^L - z): a^-1 = a(-x) * (a(x)a(-x))^-1, recursively
 * down to a scalar inverse computed by exponentiation (constant control flow).
 */
#include "poly.h"
#include "end_arith.h"

#define L END_LEAF

static inline int16_t leaf_zeta257(unsigned blk)
{
    int16_t z = end_zetas257[64 + (blk >> 1)];
    return (blk & 1) ? (int16_t)(END_Q - z) : z;
}

/* a^255 = a^-1 mod 257 (constant time; 0 -> 0) */
static int16_t inv257(int16_t a)
{
    int16_t r = a;
    unsigned i;
    for (i = 0; i < 7; i++) {
        r = red257(fqmul257(r, r));
        r = red257(fqmul257(r, a));
    }
    return r;
}

/*
 * Inverse in Z_257[x]/(x^m - z) by norm descent:  a^-1 = a(-x) * (a(x) a(-x))^-1, where
 * a(x)a(-x) only has even powers and therefore lives in Z_257[y]/(y^(m/2) - z), y = x^2.
 * Fixed control flow (constant time); 'fail' accumulates a non-zero value if some scalar is 0.
 */
static void leaf_inv257(int16_t *r, const int16_t *a, unsigned m, int16_t z, uint32_t *fail)
{
    int16_t c[L] = {0}, prod[L] = {0}, bb[L] = {0}, binv[L] = {0}, d[L] = {0};
    unsigned i, k;
    if (m == 1) {
        *fail |= ct_is_zero_u32((uint32_t)a[0]);
        r[0] = inv257(a[0]);
        return;
    }
    for (i = 0; i < m; i++)
        c[i] = (i & 1) ? red257((int16_t)(END_Q - a[i])) : a[i];
    /* prod = a * c mod (x^m - z) */
    for (k = 0; k < m; k++) {
        int32_t lo = 0, hi = 0;
        for (i = 0; i < m; i++) {
            if (i <= k) lo += (int32_t)a[i] * c[k - i];
            else        hi += (int32_t)a[i] * c[k + m - i];
        }
        prod[k] = red257_32(lo + (int32_t)z * hi);
    }
    for (i = 0; i < m / 2; i++) bb[i] = prod[2 * i];
    leaf_inv257(binv, bb, m / 2, z, fail);
    for (i = 0; i < m; i++) d[i] = (i & 1) ? 0 : binv[i / 2];
    /* r = c * d mod (x^m - z) */
    for (k = 0; k < m; k++) {
        int32_t lo = 0, hi = 0;
        for (i = 0; i < m; i++) {
            if (i <= k) lo += (int32_t)c[i] * d[k - i];
            else        hi += (int32_t)c[i] * d[k + m - i];
        }
        r[k] = red257_32(lo + (int32_t)z * hi);
    }
}

int poly_baseinv257(int16_t r[END_N], const int16_t a[END_N])
{
    unsigned blk;
    uint32_t fail = 0;
    for (blk = 0; blk < END_N / L; blk++)
        leaf_inv257(r + blk * L, a + blk * L, L, leaf_zeta257(blk), &fail);
    return (int)fail;
}
