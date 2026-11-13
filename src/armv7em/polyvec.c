/*
 * END KEM -- ARMv7E-M kernel: constant-time vector helpers on packed 32-bit words.
 *   poly_rot16 : barrel shifter; even shifts move whole words, the odd stage (shift 1) re-packs halves.
 *   sqnorm     : SMLALD (two squares per instruction into a 64-bit accumulator)
 */
#include <string.h>
#include "poly.h"
#include "end_arith.h"
#include "dsp.h"

#define NW (END_N / 2)

void poly_rot16(int16_t r[END_N], const int16_t a[END_N], uint32_t k)
{
    int32_t t[NW], u[NW];
    unsigned b, i;
    memcpy(t, a, sizeof t);
    /* stage 0: shift by one coefficient: u[i] = [t[i].lo : t[i-1].hi], u[0].lo = -t[NW-1].hi */
    {
        uint32_t m = 0u - (k & 1u);
        int32_t prev = -(t[NW - 1] >> 16);                  /* negacyclic wrap of the last coefficient */
        for (i = 0; i < NW; i++) {
            int32_t cur = t[i];
            int32_t sh = dsp_pkhbt(prev, cur);              /* [cur.lo : prev] */
            u[i] = (int32_t)(((uint32_t)sh & m) | ((uint32_t)cur & ~m));
            prev = cur >> 16;
        }
        memcpy(t, u, sizeof t);
    }
    for (b = 1; b < END_LOGN; b++) {
        unsigned sw = 1u << (b - 1);                        /* shift in words */
        uint32_t m = 0u - ((k >> b) & 1u);
        for (i = 0; i < sw; i++)   /* wrapped words are negated half-wise (packed negation) */
            u[i] = (int32_t)(((uint32_t)dsp_ssub16(0, t[i + NW - sw]) & m) | ((uint32_t)t[i] & ~m));
        for (i = sw; i < NW; i++)
            u[i] = (int32_t)(((uint32_t)t[i - sw] & m) | ((uint32_t)t[i] & ~m));
        memcpy(t, u, sizeof t);
    }
    memcpy(r, t, sizeof t);
}

uint64_t poly_sqnorm16(const int16_t a[END_N])
{
    int64_t acc = 0;
    unsigned i;
    for (i = 0; i < NW; i++) {
        int32_t pair = dsp_load_pair(a, i);
        acc = dsp_smlald(pair, pair, acc);
    }
    return (uint64_t)acc;
}

void poly_cmov16(int16_t r[END_N], const int16_t a[END_N], uint16_t mask)
{
    uint32_t m = (uint32_t)mask | ((uint32_t)mask << 16);
    unsigned i;
    for (i = 0; i < NW; i++)
        dsp_store_pair(r, i, (int32_t)(((uint32_t)dsp_load_pair(a, i) & m) | ((uint32_t)dsp_load_pair(r, i) & ~m)));
}

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
