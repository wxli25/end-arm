/*
 * END KEM -- scalar modular arithmetic and constant-time helpers (L1).
 *
 * Reduction specifications (each backend must implement the same mathematical maps):
 *   mont257(a)    : a * 2^-16 mod 257 == a mod 257 (2^16 = 1 mod 257); |a| < 2^28 -> result fits int16
 *   red257(x)     : any int16 -> [0,257)                (Fermat: x = hi*256 + lo = lo - hi mod 257)
 *   mont3329(a)   : a * 2^-16 mod 3329, |a| < 2^15*3329 -> (-3329, 3329)
 *   barrett3329(x): any int16 -> centred (-1665, 1665)
 *   mont12289(a)  : a * 2^-16 mod 12289, |a| < 2^15*12289 -> (-12289, 12289)
 *   barrett12289(x): any int16 -> centred, abs(r) <= 6145 (the value -6145 = 6144 mod 12289 can occur)
 */
#ifndef END_ARITH_H
#define END_ARITH_H

#include <stdint.h>
#include "end_params.h"
#include "end_tables.h"

/* ---- constant-time primitives --------------------------------------- */
static inline uint32_t ct_lt_u32(uint32_t a, uint32_t b)   /* 1 if a < b else 0 */
{
    return (a ^ ((a ^ b) | ((a - b) ^ b))) >> 31;
}
static inline uint32_t ct_mask_lt_u32(uint32_t a, uint32_t b) { return 0u - ct_lt_u32(a, b); }
static inline uint64_t ct_lt_u64(uint64_t a, uint64_t b)
{
    return (a ^ ((a ^ b) | ((a - b) ^ b))) >> 63;
}
static inline uint32_t ct_eq_u32(uint32_t a, uint32_t b)   /* 1 if a == b */
{
    uint32_t x = a ^ b;
    return 1u ^ ((x | (0u - x)) >> 31);
}
static inline uint32_t ct_is_zero_u32(uint32_t x) { return 1u ^ ((x | (0u - x)) >> 31); }
static inline int32_t ct_sel_i32(int32_t a, int32_t b, uint32_t mask) /* mask ? a : b (mask all-ones or 0) */
{
    return (int32_t)(((uint32_t)a & mask) | ((uint32_t)b & ~mask));
}
static inline int16_t wrap16(int32_t x) { return (int16_t)(uint16_t)((uint32_t)x); }

/* ---- generic Montgomery reduction, R = 2^16 ---------------------------- */
static inline int16_t mont_reduce16(int32_t a, int32_t p, uint16_t pinv_neg)
{
    uint32_t lo = (uint32_t)a & 0xFFFFu;                    /* unsigned arithmetic only: no UB */
    int16_t u = (int16_t)(uint16_t)(lo * (uint32_t)pinv_neg); /* u = a * (-p^-1) mod 2^16 */
    int32_t t = a + (int32_t)u * p;                          /* = 0 mod 2^16 */
    return (int16_t)(t >> 16);
}

/* ---- mod 257 ---------------------------------------------------------- */
static inline int16_t mont257(int32_t a) { return mont_reduce16(a, END_Q, (uint16_t)END_QINV16_MOD257); }
static inline int16_t red257(int16_t x)
{
    int32_t r = (x & 255) - (x >> 8);              /* lo - hi, in [-127, 383] */
    r -= END_Q & ~((r - END_Q) >> 31);             /* r >= 257 -> r - 257      */
    r += END_Q & (r >> 31);                        /* r <  0   -> r + 257      */
    return (int16_t)r;                             /* [0, 257)                 */
}
static inline int16_t fqmul257(int16_t a, int16_t b) { return mont257((int32_t)a * b); }
static inline int16_t red257_32(int32_t a)          /* |a| < 2^28 -> [0,257) */
{
    return red257(mont257(a));
}

/* ---- mod 3329 --------------------------------------------------------- */
static inline int16_t mont3329(int32_t a) { return mont_reduce16(a, END_P3, (uint16_t)END_QINV16_MOD3329); }
static inline int16_t barrett3329(int16_t x)         /* any int16 -> (-1665, 1665) */
{
    const int32_t v = ((1 << 26) + END_P3 / 2) / END_P3;
    int32_t t = ((int32_t)v * x + (1 << 25)) >> 26;
    return (int16_t)(x - t * END_P3);
}
static inline int16_t canon3329(int16_t x)           /* any int16 -> [0, 3329) */
{
    int16_t r = barrett3329(x);
    return (int16_t)(r + ((r >> 15) & END_P3));
}
static inline int16_t fqmul3329(int16_t a, int16_t b) { return mont3329((int32_t)a * b); }
static inline int16_t red3329_32(int32_t a)          /* |a| < 2^27 -> [0, 3329) */
{
    /* Barrett with M = ceil(2^39 / 3329): quotient error at most 1 in either direction */
    const int64_t M = 165140983; /* ceil(2^39/3329) */
    int32_t q = (int32_t)(((int64_t)a * M) >> 39);
    int32_t r = a - q * END_P3;              /* in (-3329, 2*3329) */
    r += (r >> 31) & END_P3;
    r -= END_P3 & ~((r - END_P3) >> 31);
    return (int16_t)r;
}

/* ---- mod 12289 -------------------------------------------------------- */
static inline int16_t mont12289(int32_t a) { return mont_reduce16(a, END_QP, (uint16_t)END_QINV16_MOD12289); }
static inline int16_t barrett12289(int16_t x)        /* any int16 -> (-6145, 6145) */
{
    const int32_t v = ((1 << 26) + END_QP / 2) / END_QP;
    int32_t t = ((int32_t)v * x + (1 << 25)) >> 26;
    return (int16_t)(x - t * END_QP);
}
static inline int16_t canon12289(int16_t x)          /* any int16 -> [0, 12289) */
{
    int16_t r = barrett12289(x);
    return (int16_t)(r + ((r >> 15) & END_QP));
}
static inline int16_t fqmul12289(int16_t a, int16_t b) { return mont12289((int32_t)a * b); }

/* ---- rounding helpers used by the decoder ----------------------------- */
/* round((v)/(2qQ)) with ties towards +inf, for |v| < 2^27, constant time (Barrett, M = ceil(2^50/2qQ)) */
static inline int32_t round_div_2qqp(int32_t v)
{
    const int64_t M = 178198389; /* ceil(2^50 / 6316546) */
    int64_t x = (int64_t)v + END_QQP + (int64_t)8 * END_2QQP;   /* positive */
    int32_t q = (int32_t)((x * M) >> 50);
    int32_t r = (int32_t)(x - (int64_t)q * END_2QQP);
    /* correct a possible +1 overestimate (r < 0) or -1 underestimate (r >= 2qQ) */
    q += (r >> 31);
    r += (r >> 31) & END_2QQP;
    q += ~((r - END_2QQP) >> 31) & 1;
    return q - 8;
}

/* centre v mod 2qQ into [-qQ, qQ) for |v| < 2*2qQ */
static inline int32_t center_2qqp(int32_t v)
{
    v += END_2QQP & ((v + END_QQP) >> 31);          /* if v < -qQ: add 2qQ */
    v -= END_2QQP & ~((v - END_QQP) >> 31);         /* if v >= qQ: subtract 2qQ */
    return v;
}

#endif /* END_ARITH_H */
