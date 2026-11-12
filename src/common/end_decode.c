/*
 * END KEM -- trapdoor decoder (L6): Alg. 6 (END.DecodeTrap) in the doubled / centred domain.
 *
 * Notation (design brief V5):  s~ = 2s - 1,  mu~ = 2 mu' - delta,  c1 = g s~ + f mu~ mod 2q (centred),
 *   C2 = (Q c2 - w c1) mod 2qQ (centred),  F_Q = Q F - w f,  G_Q = Q G - w g,
 *   s~' = (c1 F_Q - C2 f) / (qQ)  =  s~ - 2 e1 F + 2 e2 f           (exact division; computed mod 3329)
 *   gf stage : for a candidate e~ (sparse, +-1):  r = round((C2 - 2 q w e~)/(2qQ)),  s~'' = s~' + 2 e~ F + 2 r f
 *              C2~ = C2 - 2 q w e~ - 2qQ r ;  score = |s~''|^2 (+ |mu~''|^2 with END_FULL_SCORE)
 *   GF stage : for a candidate e~ (sparse, +-1):  s~ = s~'' - 2 e~ f ;  accept iff all entries are +-1.
 *
 * Every loop bound, memory address and candidate count is public; secret-dependent choices are made
 * with masks (poly_cmov*, fcl_select, poly_rot16).
 */
#include <string.h>
#include "end_internal.h"
#include "end_arith.h"

#define M0A END_M0_GF
#define M0B END_M0_GF2

/* r in [0,257) -> representative mod 514 with the given parity, centred to [-257, 256] */
static inline int16_t lift_center_parity(int16_t r, uint16_t parity)
{
    int32_t x = r;
    x += END_Q & -(int32_t)(((uint32_t)x ^ parity) & 1u);     /* fix parity: x or x + 257       */
    x -= (2 * END_Q) & ~((x - END_Q) >> 31);                   /* x >= 257 -> x - 514           */
    return (int16_t)x;
}

/* NTT-domain helper: t = INTT257(cn o k), then v = lift(2 t - cst mod 257) */
static void residue_lift(int16_t out[END_N], const int16_t cn[END_N], const int16_t kn[END_N],
                         const int16_t cst[END_N], uint16_t parity, int16_t tmp[END_N])
{
    unsigned i;
    poly_basemul257(tmp, cn, kn);
    poly_intt257(tmp);
    for (i = 0; i < END_N; i++)
        out[i] = lift_center_parity(red257((int16_t)(2 * tmp[i] - cst[i])), parity);
}

/* candidate enumeration order (must match golden/end_golden.py candidates()):
 *   zero; singles: for k in S: sigma = -1, +1 ; pairs: for k0 < k1: (-1,-1), (-1,+1), (+1,-1), (+1,+1) */
typedef struct { int8_t s0, s1; uint8_t k0, k1; } cand_t;   /* s0 = coefficient of x^S[k0], s1 of x^S[k1] (0 = absent) */

static unsigned enumerate(cand_t *c, unsigned m0)
{
    unsigned n = 0, k0, k1;
    int sg0, sg1;
    c[n].s0 = 0; c[n].s1 = 0; c[n].k0 = 0; c[n].k1 = 0; n++;
    for (k0 = 0; k0 < m0; k0++)
        for (sg0 = -1; sg0 <= 1; sg0 += 2) { c[n].s0 = (int8_t)sg0; c[n].s1 = 0; c[n].k0 = (uint8_t)k0; c[n].k1 = 0; n++; }
    for (k0 = 0; k0 < m0; k0++)
        for (k1 = k0 + 1; k1 < m0; k1++)
            for (sg0 = -1; sg0 <= 1; sg0 += 2)
                for (sg1 = -1; sg1 <= 1; sg1 += 2) {
                    c[n].s0 = (int8_t)sg0; c[n].s1 = (int8_t)sg1; c[n].k0 = (uint8_t)k0; c[n].k1 = (uint8_t)k1; n++;
                }
    return n;
}

int decode_trap(int16_t sbar[END_N], const end_expkey *ek, const int16_t cprime[END_N], uint32_t diag[8])
{
    int16_t cn[END_N] END_ALIGN, t1[END_N] END_ALIGN, t2[END_N] END_ALIGN;
    int16_t c1[END_N] END_ALIGN, c2[END_N] END_ALIGN;
    int32_t C2[END_N] END_ALIGN, C2b[END_N] END_ALIGN;
    int16_t sp[END_N] END_ALIGN;       /* s~'                              */
    int16_t best[END_N] END_ALIGN;     /* best s~'' of the gf stage        */
    int16_t rv[END_N] END_ALIGN;       /* rounding vector r of a candidate */
    int16_t wrot[M0A][END_N] END_ALIGN, Frot[M0A][END_N] END_ALIGN;
    int16_t frot[M0B][END_N] END_ALIGN;
    uint32_t sc[END_N], idx[END_FCL_MAX];
    cand_t cand[END_NCAND_GF > END_NCAND_GF2 ? END_NCAND_GF : END_NCAND_GF2];
    uint64_t best_score;
    uint32_t found, best_i, i, ci, nc;

    /* ---- step 1-3: residues --------------------------------------------- */
    memcpy(cn, cprime, sizeof cn);
    poly_ntt257(cn);
    residue_lift(c1, cn, ek->fn257, ek->const1, ek->p1, t1);     /* c1 = g s~ + f mu~ (mod 514, centred) */
    residue_lift(c2, cn, ek->Fn257, ek->const2, ek->p2, t1);     /* c2 = G s~ + F mu~ (mod 514, centred) */
    /* u = w c1 mod qQ via CRT(257, 12289); mod 257 from the NTT-257 domain, mod 12289 from NTT-12289 */
    poly_basemul257(t1, cn, ek->wfn257);
    poly_intt257(t1);                                            /* t1 = w f c' mod 257 */
    memcpy(t2, c1, sizeof t2);
    poly_ntt12289(t2);                                           /* t2 = NTT12289(c1) */
    poly_pointwise12289(t2, t2, ek->wn12289);
    poly_intt12289(t2);                                          /* t2 = w c1 mod 12289, centred */
    for (i = 0; i < END_N; i++) {
        int32_t u257 = red257((int16_t)(2 * t1[i] - ek->wconst[i]));           /* w c1 mod 257 */
        int32_t u12 = t2[i] + ((t2[i] >> 15) & END_QP);                         /* [0, 12289)   */
        int32_t d = u12 - u257;                                                  /* (-257, 12289) */
        int32_t tt = mont12289(d * END_INV257_MONT_MOD12289);                   /* d * 257^-1 mod 12289, (-Q',Q') */
        int32_t u, v;
        tt += (tt >> 31) & END_QP;                                               /* [0, 12289) */
        u = u257 + END_Q * tt;                                                   /* [0, qQ): w c1 mod qQ */
        u += END_QQP & -(int32_t)(((uint32_t)u ^ ek->pwc1) & 1u);               /* parity lift to mod 2qQ */
        v = END_QP * (int32_t)c2[i] - u;                                         /* (-3qQ, qQ) */
        C2[i] = center_2qqp(v);                                                  /* [-qQ, qQ) */
    }
    /* ---- step 4: s~' = (c1 F_Q - C2 f)/(qQ)  (mod 3329, exact) ------------- */
    memcpy(t1, c1, sizeof t1);
    poly_ntt3329(t1);
    for (i = 0; i < END_N; i++) t2[i] = red3329_32(C2[i]);
    poly_ntt3329(t2);
    poly_basemul3329(t1, t1, ek->FQn3329);
    poly_basemul3329(t2, t2, ek->fn3329);
    for (i = 0; i < END_N; i++)
        sp[i] = fqmul3329(canon3329((int16_t)(t1[i] - t2[i])), END_INV_QQP_MONT_MOD3329);
    poly_intt3329(sp);                                           /* s~' centred, exact (|s~'| <= 129) */

    /* ---- gf stage: FCL on c1 (boundary q), candidates, norm score ---------- */
    for (i = 0; i < END_N; i++) {
        int32_t a = c1[i];
        a = (a ^ (a >> 31)) - (a >> 31);                         /* |c1| */
        sc[i] = (uint32_t)(END_Q - a);
    }
    fcl_select(idx, M0A, sc);
    for (ci = 0; ci < M0A; ci++) {
        poly_rot16(wrot[ci], ek->w, idx[ci]);
        poly_rot16(Frot[ci], ek->F16, idx[ci]);
    }
    nc = enumerate(cand, M0A);
    memcpy(best, sp, sizeof best);
    memcpy(C2b, C2, sizeof C2b);
    best_score = poly_sqnorm16(sp);
    best_i = 0;
    for (ci = 1; ci < nc; ci++) {
        const cand_t *c = &cand[ci];
        const int16_t *w0 = wrot[c->k0], *w1 = wrot[c->k1], *F0 = Frot[c->k0], *F1 = Frot[c->k1];
        int32_t s0 = c->s0, s1 = c->s1;
        uint64_t score;
        uint32_t better;
        /* r = round((C2 - 2 q w e~) / (2qQ)), |r| <= 6, kept in rv; t2 = r f (exact, |r f| <= 864) */
        for (i = 0; i < END_N; i++) {
            int32_t we = s0 * w0[i] + s1 * w1[i];
            rv[i] = (int16_t)round_div_2qqp(C2[i] - 2 * END_Q * we);
            t2[i] = rv[i];
        }
        poly_ntt3329(t2);
        poly_basemul3329(t2, t2, ek->fn3329);
        poly_intt3329(t2);
        for (i = 0; i < END_N; i++)
            t1[i] = (int16_t)(sp[i] + 2 * (s0 * F0[i] + s1 * F1[i]) + 2 * t2[i]);   /* s~'' */
        score = poly_sqnorm16(t1);
        better = (uint32_t)ct_lt_u64(score, best_score);
        best_score = (score & (0ull - (uint64_t)better)) | (best_score & ~(0ull - (uint64_t)better));
        best_i = (ci & (0u - better)) | (best_i & ~(0u - better));
        poly_cmov16(best, t1, (uint16_t)(0u - better));
        /* C2~ = C2 - 2 q w e~ - 2qQ r, selected with the same mask */
        for (i = 0; i < END_N; i++) {
            int32_t we = s0 * w0[i] + s1 * w1[i];
            int32_t v = C2[i] - 2 * END_Q * we - (int32_t)rv[i] * END_2QQP;
            C2b[i] = ct_sel_i32(v, C2b[i], 0u - better);
        }
    }
    if (diag) { diag[0] = best_i; diag[2] = idx[0]; diag[3] = (M0A > 1) ? idx[1] : 0; }

    /* ---- GF stage: FCL on C2~ (boundary qQ), candidates, binary check ------- */
    for (i = 0; i < END_N; i++) {
        int32_t a = C2b[i];
        a = (a ^ (a >> 31)) - (a >> 31);
        sc[i] = (uint32_t)(END_QQP - a);
    }
    fcl_select(idx, M0B, sc);
    for (ci = 0; ci < M0B; ci++)
        poly_rot16(frot[ci], ek->f16, idx[ci]);
    nc = enumerate(cand, M0B);
    memcpy(sbar, best, sizeof best);     /* fallback: uncorrected s~'' (decryption failure -> implicit rejection) */
    found = 0; best_i = 0;
    for (ci = 0; ci < nc; ci++) {
        const cand_t *c = &cand[ci];
        const int16_t *f0 = frot[c->k0], *f1 = frot[c->k1];
        int32_t s0 = c->s0, s1 = c->s1;
        uint32_t acc = 0, pass, take;
        for (i = 0; i < END_N; i++) {
            int32_t x = best[i] - 2 * (s0 * f0[i] + s1 * f1[i]);
            t1[i] = (int16_t)x;
            acc |= (uint32_t)(x * x - 1);                       /* 0 iff x = +-1 */
        }
        pass = ct_is_zero_u32(acc);
        take = pass & ~found;
        poly_cmov16(sbar, t1, (uint16_t)(0u - take));
        best_i = (ci & (0u - take)) | (best_i & ~(0u - take));
        found |= pass;
    }
    if (diag) { diag[1] = best_i; diag[4] = idx[0]; diag[5] = (M0B > 1) ? idx[1] : 0; diag[6] = found; }
    return (int)found;
}
