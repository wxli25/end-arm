/*
 * END KEM -- ARMv7E-M / ARMv8-M Mainline DSP helpers (inline asm wrappers; no dependency on arm_acle.h).
 *
 * Packed representation: one 32-bit word holds two consecutive int16 coefficients
 * (low half = even index, high half = odd index); all buffers are 4-byte aligned (END_ALIGN).
 *
 * Montgomery multiplication of a packed pair y by a twiddle z (R = 2^16), zq = z * p^-1 mod 2^16:
 *   t = y_i * z ; m = low16(y_i * zq) ; r_i = (t - m p) >> 16      (SMULxx, SMLAxx, PKHTB)
 */
#ifndef END_DSP_H
#define END_DSP_H

#include <stdint.h>
#include <string.h>

/* Copy the representation of a coefficient pair without violating C aliasing rules.
 * __builtin_memcpy (not memcpy): pqm4 compiles with -ffreestanding, where a plain memcpy
 * becomes a library call per coefficient pair; the builtin always folds to one LDR/STR. */
static inline int32_t dsp_load_pair(const int16_t *a, unsigned word)
{
    int32_t value;
    __builtin_memcpy(&value, a + 2 * word, sizeof value);
    return value;
}

static inline void dsp_store_pair(int16_t *a, unsigned word, int32_t value)
{
    __builtin_memcpy(a + 2 * word, &value, sizeof value);
}

#define DSP2(name, insn) \
    static inline int32_t name(int32_t a, int32_t b) { int32_t r; __asm__(insn " %0, %1, %2" : "=r"(r) : "r"(a), "r"(b)); return r; }
#define DSP3(name, insn) \
    static inline int32_t name(int32_t a, int32_t b, int32_t c) { int32_t r; __asm__(insn " %0, %1, %2, %3" : "=r"(r) : "r"(a), "r"(b), "r"(c)); return r; }

DSP2(dsp_smulbb, "smulbb")   /* a.lo * b.lo */
DSP2(dsp_smulbt, "smulbt")   /* a.lo * b.hi */
DSP2(dsp_smultb, "smultb")   /* a.hi * b.lo */
DSP2(dsp_smultt, "smultt")   /* a.hi * b.hi */
DSP3(dsp_smlabb, "smlabb")   /* c + a.lo * b.lo */
DSP3(dsp_smlabt, "smlabt")
DSP3(dsp_smlatb, "smlatb")
DSP3(dsp_smlatt, "smlatt")
DSP2(dsp_smuad,  "smuad")    /* a.lo*b.lo + a.hi*b.hi */
DSP2(dsp_smuadx, "smuadx")   /* a.lo*b.hi + a.hi*b.lo */
DSP3(dsp_smlad,  "smlad")    /* c + a.lo*b.lo + a.hi*b.hi */
DSP3(dsp_smladx, "smladx")
DSP2(dsp_sadd16, "sadd16")   /* packed 16-bit add (wrapping) */
DSP2(dsp_ssub16, "ssub16")
DSP2(dsp_uadd16, "uadd16")

static inline int64_t dsp_smlald(int32_t a, int32_t b, int64_t acc)   /* acc + a.lo*b.lo + a.hi*b.hi (64-bit) */
{
    __asm__("smlald %Q0, %R0, %1, %2" : "+r"(acc) : "r"(a), "r"(b));
    return acc;
}
static inline int32_t dsp_pkhbt(int32_t lo, int32_t hi)   /* [hi.lo : lo.lo] */
{
    int32_t r; __asm__("pkhbt %0, %1, %2, lsl #16" : "=r"(r) : "r"(lo), "r"(hi)); return r;
}
static inline int32_t dsp_pkhtb(int32_t hi, int32_t lo)   /* [hi.hi : lo.hi]  (lo arithmetically shifted right by 16) */
{
    int32_t r; __asm__("pkhtb %0, %1, %2, asr #16" : "=r"(r) : "r"(hi), "r"(lo)); return r;
}
static inline int32_t dsp_lo(int32_t x) { return (int16_t)x; }
static inline int32_t dsp_hi(int32_t x) { return x >> 16; }

/* Montgomery product of a packed pair with the twiddle z (zq = z*p^-1 mod 2^16, negp = -p), result in (-p,p) */
static inline int32_t dsp_montmul_pair(int32_t y, int32_t z, int32_t zq, int32_t negp)
{
    int32_t t0 = dsp_smulbb(y, z), t1 = dsp_smultb(y, z);
    int32_t m0 = dsp_smulbb(y, zq), m1 = dsp_smultb(y, zq);
    int32_t r0 = dsp_smlabb(m0, negp, t0);          /* = (t0 - m0 p), low 16 bits zero */
    int32_t r1 = dsp_smlabb(m1, negp, t1);
    return dsp_pkhtb(r1, r0);
}

#endif /* END_DSP_H */
