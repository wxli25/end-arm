/*
 * END KEM -- polynomial kernel interface (L2).
 *
 * Every backend (src/ref, src/aarch64, src/armv7em) implements exactly these functions with
 * the same input/output contracts. The NTT-257 layout and its canonical [0, 257) values are fixed
 * because they are the public-key wire format (reference layout: 7 Cooley-Tukey layers, twiddle
 * k = 3^brv7(k), leaves of END_LEAF coefficients, leaf blk is modulo x^LEAF - (+/-)zetas257[64 + blk/2]).
 *
 * The NTT-3329 and NTT-12289 domains are internal. Their values are plain residues (no hidden
 * factor), but a backend chooses which representative it stores: the assembly backend returns
 * unreduced int16 values from the forward transforms. A caller may pass such data to the product
 * and inverse-transform functions of the same backend; before relying on a range it must reduce
 * (canon3329 / barrett12289). The input bounds below are what every backend accepts.
 *
 * All functions are constant time with respect to the coefficient values.
 */
#ifndef END_POLY_H
#define END_POLY_H

#include <stdint.h>
#include "end_params.h"

/* convenient short names inside the library */
#define poly_ntt257 END_NS(poly_ntt257)
#define poly_intt257 END_NS(poly_intt257)
#define poly_basemul257 END_NS(poly_basemul257)
#define poly_baseinv257 END_NS(poly_baseinv257)
#define poly_ntt3329 END_NS(poly_ntt3329)
#define poly_intt3329 END_NS(poly_intt3329)
#define poly_basemul3329 END_NS(poly_basemul3329)
#define poly_ntt12289 END_NS(poly_ntt12289)
#define poly_intt12289 END_NS(poly_intt12289)
#define poly_pointwise12289 END_NS(poly_pointwise12289)
#define poly_rot16 END_NS(poly_rot16)
#define poly_sqnorm16 END_NS(poly_sqnorm16)
#define poly_cmov16 END_NS(poly_cmov16)
#define fcl_select END_NS(fcl_select)
#define poly_lift257 END_NS(poly_lift257)
#define poly_crt_c2 END_NS(poly_crt_c2)
#define poly_cand_round END_NS(poly_cand_round)
#define poly_cand_combine END_NS(poly_cand_combine)
#define poly_cand_check END_NS(poly_cand_check)
#define poly_cmov32 END_NS(poly_cmov32)


/* ---------------- NTT mod 257 (incomplete, 7 layers) ---------------- */
/* in: |a| <= 16384;                     out: canonical [0,257), NTT layout                    */
void poly_ntt257(int16_t a[END_N]);
/* in: |a| <= 4095 in NTT layout;        out: canonical [0,257), coefficient order             */
void poly_intt257(int16_t a[END_N]);
/* r = a * b in the NTT domain (leaf-wise products).  in: [0,257);  out: [0,257)               */
void poly_basemul257(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N]);
/* r = a^-1 in the NTT domain. returns 0 if invertible, non-zero otherwise. in/out: [0,257)   */
int  poly_baseinv257(int16_t r[END_N], const int16_t a[END_N]);

/* ---------------- NTT mod 3329 (incomplete, 7 layers, same tree as mod 257) ---------------- */
/* in: |a| <= 16384;                     out: NTT domain (backend-defined representatives)     */
void poly_ntt3329(int16_t a[END_N]);
/* in: output of poly_basemul3329, or NTT-domain values with |a| <= 4095;
 * out: centred (-1665, 1665) -- exact if |true| <= 1664                                      */
void poly_intt3329(int16_t a[END_N]);
/* r = a * b leaf-wise.  in: outputs of poly_ntt3329;  out: a representative with |r| < 3329   */
void poly_basemul3329(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N]);

/* ---------------- NTT mod 12289 (complete) ---------------- */
/* in: |a| <= 8191;                      out: NTT domain (backend-defined representatives)     */
void poly_ntt12289(int16_t a[END_N]);
/* in: output of poly_pointwise12289, or NTT-domain values with |a| <= 8191;
 * out: centred (-6145, 6145)                                                                 */
void poly_intt12289(int16_t a[END_N]);
/* r = a o b (true product mod 12289).  in: outputs of poly_ntt12289;  out: |r| < 12289        */
void poly_pointwise12289(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N]);

/* ---------------- vector helpers (constant time in the secret index k) ---------------- */
/* r = x^k * a (negacyclic rotation), 0 <= k < END_N, k secret; r must not overlap a          */
void poly_rot16(int16_t r[END_N], const int16_t a[END_N], uint32_t k);
/* sum of squares (exact) */
uint64_t poly_sqnorm16(const int16_t a[END_N]);
/* r = mask ? a : r   (mask is 0 or all-ones) */
void poly_cmov16(int16_t r[END_N], const int16_t a[END_N], uint16_t mask);
/* Free Candidate Localization: indices of the m0 smallest scores (ties: smaller index first); scores < 2^22 */
void fcl_select(uint32_t idx[], unsigned m0, const uint32_t scores[END_N]);

/* ---------------- inner loops of the trapdoor decoder (portable versions in end_deckern.c) ----------------
 * A candidate is e~ = s0 x^k0 + s1 x^k1 with s0, s1 in {-1, 0, 1}. The candidate list is public, so s0 and s1
 * are public; every coefficient value is secret and all functions are constant time in them.               */
/* out = representative mod 2q, centred in [-q, q), of (2 t - cst) mod q with the parity bit `parity`; t, cst in [0, q) */
void poly_lift257(int16_t out[END_N], const int16_t t[END_N], const int16_t cst[END_N], uint16_t parity);
/* C2 = (Q c2 - u) centred mod 2qQ, where u = w c1 mod 2qQ is rebuilt by CRT from u mod q = (2 t257 - wconst) mod q,
 * u mod Q = t12289 (centred) and its parity pwc1.  t257, wconst in [0, q); c2 in [-q, q)                         */
void poly_crt_c2(int32_t C2[END_N], const int16_t t257[END_N], const int16_t t12289[END_N], const int16_t c2[END_N],
                 const int16_t wconst[END_N], uint16_t pwc1);
/* v = C2 - 2 q (s0 w0 + s1 w1);  r = round(v / 2qQ) (ties up);  c2t = v - 2qQ r in [-qQ, qQ).   |C2| <= qQ      */
void poly_cand_round(int16_t r[END_N], int32_t c2t[END_N], const int32_t C2[END_N],
                     const int16_t w0[END_N], const int16_t w1[END_N], int32_t s0, int32_t s1);
/* t = sp + 2 (s0 F0 + s1 F1) + 2 rf (no overflow: |sp|, |rf| <= 2047, |F| <= 127); returns sum t_i^2            */
uint64_t poly_cand_combine(int16_t t[END_N], const int16_t sp[END_N], const int16_t F0[END_N], const int16_t F1[END_N],
                           const int16_t rf[END_N], int32_t s0, int32_t s1);
/* t = b - 2 (s0 f0 + s1 f1) (|b| <= 16383, |f| <= 127); returns 0 iff every t_i is +1 or -1                     */
uint32_t poly_cand_check(int16_t t[END_N], const int16_t b[END_N], const int16_t f0[END_N], const int16_t f1[END_N],
                         int32_t s0, int32_t s1);
/* r = mask ? a : r   (mask is 0 or all-ones) */
void poly_cmov32(int32_t r[END_N], const int32_t a[END_N], uint32_t mask);

#endif /* END_POLY_H */
