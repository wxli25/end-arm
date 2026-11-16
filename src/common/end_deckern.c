/*
 * END KEM -- portable versions of the inner loops of the trapdoor decoder and of the FCL selection (poly.h).
 *
 * A backend that provides its own version of one of these functions defines the matching END_ASM_* macro in
 * its end_backend.h (src/<backend>/end_backend.h; src/ref/end_backend.h is the empty default).
 * All loops have public bounds and addresses and no data-dependent branches.
 */
#include "end_params.h"
#include "end_arith.h"
#include "poly.h"
#include "end_backend.h"

#ifndef END_ASM_FCL_SELECT
/* The pair (score, index) is packed into one 32-bit key, score in the high bits: keys are distinct, their order is
 * the lexicographic order, and a compare-exchange is one unsigned comparison. best[] is kept sorted by insertion. */
void fcl_select(uint32_t idx[], unsigned m0, const uint32_t scores[END_N])
{
    uint32_t best[END_FCL_MAX];
    unsigned i, j;
    for (j = 0; j < m0; j++) best[j] = 0xFFFFFFFFu;
    for (i = 0; i < END_N; i++) {
        uint32_t cs = (scores[i] << END_LOGN) | i;
        for (j = 0; j < m0; j++) {
            /* the top word of the 64-bit difference is all-ones exactly when cs < best[j] (a borrow) */
            uint32_t d = (cs - best[j]) & (uint32_t)(((uint64_t)cs - best[j]) >> 32);
            best[j] += d;                       /* min(best[j], cs) */
            cs -= d;                            /* max(best[j], cs) moves on */
        }
    }
    for (j = 0; j < m0; j++) idx[j] = best[j] & (END_N - 1);
}
#endif

#ifndef END_ASM_POLY_LIFT257
/* r in [0,257) -> representative mod 514 with the given parity, centred to [-257, 256] */
static inline int16_t lift_center_parity(int16_t r, uint16_t parity)
{
    int32_t x = r;
    x += END_Q & -(int32_t)(((uint32_t)x ^ parity) & 1u);     /* fix parity: x or x + 257       */
    x -= (2 * END_Q) & ~((x - END_Q) >> 31);                   /* x >= 257 -> x - 514           */
    return (int16_t)x;
}

void poly_lift257(int16_t out[END_N], const int16_t t[END_N], const int16_t cst[END_N], uint16_t parity)
{
    unsigned i;
    for (i = 0; i < END_N; i++)
        out[i] = lift_center_parity(red257((int16_t)(2 * t[i] - cst[i])), parity);
}
#endif

#ifndef END_ASM_POLY_CRT_C2
void poly_crt_c2(int32_t C2[END_N], const int16_t t257[END_N], const int16_t t12289[END_N], const int16_t c2[END_N],
                 const int16_t wconst[END_N], uint16_t pwc1)
{
    unsigned i;
    for (i = 0; i < END_N; i++) {
        int32_t u257 = red257((int16_t)(2 * t257[i] - wconst[i]));              /* w c1 mod 257 */
        int32_t u12 = t12289[i] + ((t12289[i] >> 15) & END_QP);                 /* [0, 12289)   */
        int32_t d = u12 - u257;                                                  /* (-257, 12289) */
        int32_t tt = mont12289(d * END_INV257_MONT_MOD12289);                   /* d * 257^-1 mod 12289, (-Q',Q') */
        int32_t u, v;
        tt += (tt >> 31) & END_QP;                                               /* [0, 12289) */
        u = u257 + END_Q * tt;                                                   /* [0, qQ): w c1 mod qQ */
        u += END_QQP & -(int32_t)(((uint32_t)u ^ pwc1) & 1u);                   /* parity lift to mod 2qQ */
        v = END_QP * (int32_t)c2[i] - u;                                         /* (-3qQ, qQ) */
        C2[i] = center_2qqp(v);                                                  /* [-qQ, qQ) */
    }
}
#endif

#ifndef END_ASM_POLY_CAND_ROUND
void poly_cand_round(int16_t r[END_N], int32_t c2t[END_N], const int32_t C2[END_N],
                     const int16_t w0[END_N], const int16_t w1[END_N], int32_t s0, int32_t s1)
{
    unsigned i;
    for (i = 0; i < END_N; i++) {
        int32_t we = s0 * w0[i] + s1 * w1[i];
        int32_t c;
        r[i] = (int16_t)round_div_2qqp_rem(C2[i] - 2 * END_Q * we, &c);          /* |r| <= 6 */
        c2t[i] = c;
    }
}
#endif

#ifndef END_ASM_POLY_CAND_COMBINE
uint64_t poly_cand_combine(int16_t t[END_N], const int16_t sp[END_N], const int16_t F0[END_N], const int16_t F1[END_N],
                           const int16_t rf[END_N], int32_t s0, int32_t s1)
{
    uint64_t acc = 0;
    unsigned i;
    for (i = 0; i < END_N; i++) {
        int32_t x = (int16_t)(sp[i] + 2 * (s0 * F0[i] + s1 * F1[i]) + 2 * rf[i]);
        t[i] = (int16_t)x;
        acc += (uint32_t)(x * x);
    }
    return acc;
}
#endif

#ifndef END_ASM_POLY_CAND_CHECK
uint32_t poly_cand_check(int16_t t[END_N], const int16_t b[END_N], const int16_t f0[END_N], const int16_t f1[END_N],
                         int32_t s0, int32_t s1)
{
    uint32_t acc = 0;
    unsigned i;
    for (i = 0; i < END_N; i++) {
        int32_t x = b[i] - 2 * (s0 * f0[i] + s1 * f1[i]);
        t[i] = (int16_t)x;
        acc |= (uint32_t)(x * x - 1);                       /* 0 iff x = +-1 */
    }
    return acc;
}
#endif

#ifndef END_ASM_POLY_CMOV32
void poly_cmov32(int32_t r[END_N], const int32_t a[END_N], uint32_t mask)
{
    unsigned i;
    for (i = 0; i < END_N; i++)
        r[i] = ct_sel_i32(a[i], r[i], mask);
}
#endif
