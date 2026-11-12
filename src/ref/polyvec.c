/*
 * END KEM -- portable C kernel: constant-time vector helpers
 *   poly_rot16   : negacyclic rotation by a secret amount (log2(n)-stage barrel shifter, public addresses)
 *   poly_sqnorm16: exact squared Euclidean norm
 *   poly_cmov*   : masked copy
 *   fcl_select   : constant-time selection of the m0 smallest scores (Free Candidate Localization)
 */
#include <string.h>
#include "poly.h"
#include "end_arith.h"

void poly_rot16(int16_t r[END_N], const int16_t a[END_N], uint32_t k)
{
    int16_t t[END_N], u[END_N];
    unsigned b, i;
    memcpy(t, a, sizeof t);
    for (b = 0; b < END_LOGN; b++) {
        unsigned sh = 1u << b;
        uint16_t m = (uint16_t)(0u - ((k >> b) & 1u));
        /* u = x^sh * t */
        for (i = 0; i < sh; i++)
            u[i] = (int16_t)(-t[i + END_N - sh]);
        for (i = sh; i < END_N; i++)
            u[i] = t[i - sh];
        for (i = 0; i < END_N; i++)
            t[i] = (int16_t)((u[i] & m) | (t[i] & (uint16_t)~m));
    }
    memcpy(r, t, sizeof t);
}

uint64_t poly_sqnorm16(const int16_t a[END_N])
{
    uint64_t s = 0;
    unsigned i;
    for (i = 0; i < END_N; i++)
        s += (uint64_t)((int32_t)a[i] * a[i]);
    return s;
}

void poly_cmov16(int16_t r[END_N], const int16_t a[END_N], uint16_t mask)
{
    unsigned i;
    for (i = 0; i < END_N; i++)
        r[i] = (int16_t)((a[i] & mask) | (r[i] & (uint16_t)~mask));
}

void poly_cmov32(int32_t r[END_N], const int32_t a[END_N], uint32_t mask)
{
    unsigned i;
    for (i = 0; i < END_N; i++)
        r[i] = (int32_t)(((uint32_t)a[i] & mask) | ((uint32_t)r[i] & ~mask));
}

/*
 * Keep a sorted list (ascending score) of the m0 best entries; every input element is
 * inserted with a fixed chain of masked swaps, so the running time and the memory access
 * pattern are independent of the scores.  Ties resolve towards the smaller index.
 */
void fcl_select(uint32_t idx[], unsigned m0, const uint32_t scores[END_N])
{
    uint32_t best[END_FCL_MAX], bidx[END_FCL_MAX];
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
