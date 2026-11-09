/*
 * END KEM -- polynomial kernel interface (L2).
 *
 * Every backend (src/ref, src/aarch64, src/armv7em) implements exactly these functions with
 * the same input/output contracts. Only the *internal* ordering of NTT-3329 / NTT-12289 data
 * may differ between backends; the NTT-257 layout is fixed because it is the public-key wire
 * layout (reference layout: 7 Cooley-Tukey layers, twiddle k = 3^brv7(k), leaves of
 * END_LEAF coefficients, leaf blk is modulo x^LEAF - (+/-)zetas257[64 + blk/2]).
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
#define poly_cmov32 END_NS(poly_cmov32)
#define fcl_select END_NS(fcl_select)


/* ---------------- NTT mod 257 (incomplete, 7 layers) ---------------- */
/* in: any int16 (reduced internally);  out: canonical [0,257), NTT layout                     */
void poly_ntt257(int16_t a[END_N]);
/* in: any int16 in NTT layout;        out: canonical [0,257), coefficient order              */
void poly_intt257(int16_t a[END_N]);
/* r = a * b in the NTT domain (leaf-wise products).  in: [0,257);  out: [0,257)               */
void poly_basemul257(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N]);
/* r = a^-1 in the NTT domain. returns 0 if invertible, non-zero otherwise. in/out: [0,257)   */
int  poly_baseinv257(int16_t r[END_N], const int16_t a[END_N]);

/* ---------------- NTT mod 3329 (incomplete, 7 layers, same tree as mod 257) ---------------- */
/* in: any int16;                      out: canonical [0,3329)                                */
void poly_ntt3329(int16_t a[END_N]);
/* in: any int16 (NTT domain);         out: centred (-1665, 1665) -- exact if |true| <= 1664   */
void poly_intt3329(int16_t a[END_N]);
/* r = a * b leaf-wise.  in: [0,3329);  out: [0,3329)                                          */
void poly_basemul3329(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N]);

/* ---------------- NTT mod 12289 (complete) ---------------- */
/* in: any int16;                      out: canonical [0,12289)                               */
void poly_ntt12289(int16_t a[END_N]);
/* in: any int16 (NTT domain);         out: centred (-6145, 6145)                              */
void poly_intt12289(int16_t a[END_N]);
/* r = a o b (true product mod 12289, no Montgomery factor).  in: [0,12289);  out: [0,12289)   */
void poly_pointwise12289(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N]);

/* ---------------- vector helpers (constant time in the secret index k) ---------------- */
/* r = x^k * a (negacyclic rotation), 0 <= k < END_N, k secret                               */
void poly_rot16(int16_t r[END_N], const int16_t a[END_N], uint32_t k);
/* sum of squares (exact) */
uint64_t poly_sqnorm16(const int16_t a[END_N]);
/* r = mask ? a : r   (mask is 0 or all-ones) */
void poly_cmov16(int16_t r[END_N], const int16_t a[END_N], uint16_t mask);
void poly_cmov32(int32_t r[END_N], const int32_t a[END_N], uint32_t mask);
/* Free Candidate Localization: indices of the m0 smallest scores (ties: smaller index first) */
void fcl_select(uint32_t idx[], unsigned m0, const uint32_t scores[END_N]);

#endif /* END_POLY_H */
