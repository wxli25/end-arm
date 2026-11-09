/*
 * END KEM -- byte-level codecs (L3), bit-exact with the reference implementation (dde0a77).
 *
 * Public key : NTT-257 domain coefficients taken modulo 258 = 86 * 3 (CRT), the residues mod 86
 *              packed 7 per 45-bit word, the residues mod 3 packed 17 per 27-bit word, with the
 *              exact bit layout of the reference code.c  (514 / 1027 bytes).
 * Ciphertext : END-512: 4 x 6-bit coefficients in 3 bytes.  END-1024: 8 base-90 digits per 52-bit word,
 *              128 words packed into 104 little-endian uint64 (832 bytes).
 * Message    : s[i*MB + j] = bit i of m'[j] ; s[n/2 + i*MB + j] = bit i of r[j].
 *
 * All decoders are total functions (any byte string decodes to in-range coefficients) and run in
 * constant time; divisions use reciprocal multiplication with limb-wise correction (no 128-bit types).
 */
#include <string.h>
#include "end_internal.h"
#include "end_arith.h"

static inline uint32_t ct_ge_mask_u32(uint32_t x, uint32_t y)   /* all-ones if x >= y (x, y < 2^31) */
{
    return 0u - (1u ^ ((x - y) >> 31));
}
static inline uint32_t ct_sub_if_ge_u32(uint32_t x, uint32_t y)
{
    return x - (y & ct_ge_mask_u32(x, y));
}
static inline uint16_t ct_mod258_i16(int16_t a)   /* canonical residue mod 258 of a signed int16 */
{
    uint32_t x = (uint32_t)((int32_t)a + 33024);
    x = ct_sub_if_ge_u32(x, 33024u); x = ct_sub_if_ge_u32(x, 16512u); x = ct_sub_if_ge_u32(x, 8256u);
    x = ct_sub_if_ge_u32(x, 4128u);  x = ct_sub_if_ge_u32(x, 2064u);  x = ct_sub_if_ge_u32(x, 1032u);
    x = ct_sub_if_ge_u32(x, 516u);   x = ct_sub_if_ge_u32(x, 258u);
    return (uint16_t)x;
}
static inline uint16_t ct_mod258_u32_small(uint32_t x)
{
    x = ct_sub_if_ge_u32(x, 4128u); x = ct_sub_if_ge_u32(x, 2064u); x = ct_sub_if_ge_u32(x, 1032u);
    x = ct_sub_if_ge_u32(x, 516u);  x = ct_sub_if_ge_u32(x, 258u);
    return (uint16_t)x;
}
static inline uint16_t ct_mod86_u16(uint16_t a)   /* a < 258 */
{
    uint32_t x = a;
    x = ct_sub_if_ge_u32(x, 172u); x = ct_sub_if_ge_u32(x, 86u);
    return (uint16_t)x;
}
static inline uint16_t ct_mod3_u16(uint16_t a)    /* a < 258 */
{
    uint32_t x = a;
    x = ct_sub_if_ge_u32(x, 192u); x = ct_sub_if_ge_u32(x, 96u); x = ct_sub_if_ge_u32(x, 48u);
    x = ct_sub_if_ge_u32(x, 24u);  x = ct_sub_if_ge_u32(x, 12u); x = ct_sub_if_ge_u32(x, 6u);
    x = ct_sub_if_ge_u32(x, 3u);
    return (uint16_t)x;
}

/* x / d and x mod d for a 64-bit x and a small constant d, limb-wise (16 bits per step), constant time.
 * mfloor = floor(2^32 / d). Each step's quotient estimate is floor(t/d) or floor(t/d)-1 (t < d*2^16). */
static inline void ct_divmod_small_u64(uint64_t x, uint32_t d, uint32_t mfloor, uint64_t *q_out, uint16_t *r_out)
{
    uint64_t q = 0;
    uint32_t r = 0;
    int s;
    for (s = 48; s >= 0; s -= 16) {
        uint32_t t = (r << 16) | (uint32_t)((x >> s) & 0xFFFFu);
        uint32_t qq = (uint32_t)(((uint64_t)t * mfloor) >> 32);
        uint32_t m;
        r = t - qq * d;
        m = ct_ge_mask_u32(r, d);
        r -= d & m;
        qq += m & 1u;
        q = (q << 16) | qq;
    }
    *q_out = q;
    *r_out = (uint16_t)r;
}
#define MFLOOR86 49941480u   /* floor(2^32/86) */
#define MFLOOR90 47721858u   /* floor(2^32/90) */
#define MFLOOR3  1431655765u /* floor(2^32/3)  */

static inline void store64_le(uint8_t *out, uint64_t x) { int i; for (i = 0; i < 8; i++) out[i] = (uint8_t)(x >> (8 * i)); }
static inline uint64_t load64_le(const uint8_t *in) { uint64_t x = 0; int i; for (i = 7; i >= 0; i--) x = (x << 8) | in[i]; return x; }
static inline void store32_le(uint8_t *out, uint32_t x) { int i; for (i = 0; i < 4; i++) out[i] = (uint8_t)(x >> (8 * i)); }
static inline uint32_t load32_le(const uint8_t *in) { return (uint32_t)in[0] | ((uint32_t)in[1] << 8) | ((uint32_t)in[2] << 16) | ((uint32_t)in[3] << 24); }

static const uint64_t pk_table1[7] = { 1ULL, 86ULL, 7396ULL, 636056ULL, 54700816ULL, 4704270176ULL, 404567235136ULL };
static const uint32_t pk_table2[17] = { 1U, 3U, 9U, 27U, 81U, 243U, 729U, 2187U, 6561U, 19683U, 59049U, 177147U,
                                        531441U, 1594323U, 4782969U, 14348907U, 43046721U };

/* ------------------------------------------------------------------------- */
/* public key                                                                 */
/* ------------------------------------------------------------------------- */
#if END_N == 512
#define PK_NLOW   72      /* base-86 words stored in the 64-bit low section */
#define PK_NW86   73
#define PK_HB_OFF 384
#define PK_M3_OFF 412
#define PK_NW3    30
#else
#define PK_NLOW   144
#define PK_NW86   147
#define PK_HB_OFF 768
#define PK_M3_OFF 823
#define PK_NW3    60
#endif

void pk_encode(uint8_t out[END_PKBYTES], const int16_t h[END_N])
{
    uint64_t tmp[PK_NW86];
    uint32_t tmp32[PK_NW3 + 1];
    uint32_t t;
    size_t i, j, idx;

    /* residues mod 86, 7 per word */
    for (i = 0; i < PK_NLOW + (END_N == 512 ? 1 : 2); i++) {
        uint64_t acc = 0;
        for (j = 0; j < 7; j++) {
            size_t c = 7 * i + j;
            if (c < END_N - (END_N == 512 ? 1 : 2) || END_N == 512) {
                if (c < (size_t)(END_N - (END_N == 512 ? 1 : 2)))
                    acc += (uint64_t)ct_mod86_u16(ct_mod258_i16(h[c])) * pk_table1[j];
            }
        }
        tmp[i] = acc;
    }
#if END_N == 1024
    tmp[146] = (uint64_t)ct_mod86_u16(ct_mod258_i16(h[END_N - 2])) * pk_table1[0]
             + (uint64_t)ct_mod86_u16(ct_mod258_i16(h[END_N - 1])) * pk_table1[1];
#endif
    /* low 38 bits of words 0..PK_NLOW-1 (words 1,2 mod 3 are stored fully) */
    j = 0;
    for (i = 0; i < PK_NLOW; i += 3) {
        uint64_t r0 = (tmp[i + 1] << 19) | (tmp[i] & 0x7FFFFULL);
        uint64_t r1 = (tmp[i + 2] << 19) | ((tmp[i] >> 19) & 0x7FFFFULL);
        store64_le(out + 8 * j, r0);
        store64_le(out + 8 * (j + 1), r1);
        j += 2;
    }
    /* high 7 bits of words 0,3,6,... */
    idx = PK_HB_OFF;
    for (i = 0; i < PK_NLOW; i += 24) {
        for (j = 0; j < 7; j++)
            out[idx + j] = (uint8_t)(((tmp[i + 3 * j] >> 38) << 1) | ((tmp[i + 21] >> (38 + j)) & 1ULL));
        idx += 7;
    }
#if END_N == 512
    for (i = 0; i < 6; i++) out[idx++] = (uint8_t)(tmp[72] >> (8 * i));
    out[idx++] = (uint8_t)ct_mod86_u16(ct_mod258_i16(h[END_N - 1]));
#else
    for (i = 0; i < 5; i++) out[idx++] = (uint8_t)(tmp[144] >> (8 * i));
    for (i = 0; i < 5; i++) out[idx++] = (uint8_t)(tmp[145] >> (8 * i));
    out[idx++] = (uint8_t)((((tmp[144] >> 40) & 31ULL) << 3) | (tmp[146] & 7ULL));
    out[idx++] = (uint8_t)((((tmp[145] >> 40) & 31ULL) << 3) | ((tmp[146] >> 3) & 7ULL));
    out[idx++] = (uint8_t)(tmp[146] >> 6);
#endif
    /* residues mod 3, 17 per word */
    for (i = 0; i < PK_NW3; i++) {
        uint32_t acc = 0;
        for (j = 0; j < 17; j++)
            acc += (uint32_t)ct_mod3_u16(ct_mod258_i16(h[17 * i + j])) * pk_table2[j];
        tmp32[i] = acc;
    }
    tmp32[PK_NW3] = 0;
    for (i = 0; i < (END_N == 512 ? 2u : 4u); i++)
        tmp32[PK_NW3] += (uint32_t)ct_mod3_u16(ct_mod258_i16(h[17 * PK_NW3 + i])) * pk_table2[i];
    j = 0;
    t = (END_N == 512) ? 0u : tmp32[PK_NW3];
    for (i = 0; i < PK_NW3; i += 6) {
        size_t k;
        for (k = 0; k < 5; k++)
            store32_le(out + PK_M3_OFF + 4 * (j + k), tmp32[i + k] | (((tmp32[i + 5] >> (5 * k)) & 31U) << 27));
        j += 5;
        t = (t << 2) | (tmp32[i + 5] >> 25);
    }
#if END_N == 512
    out[512] = (uint8_t)t;
    out[513] = (uint8_t)((t >> 8) | (tmp32[30] << 2));
#else
    out[1023] = (uint8_t)t; out[1024] = (uint8_t)(t >> 8); out[1025] = (uint8_t)(t >> 16); out[1026] = (uint8_t)(t >> 24);
#endif
}

void pk_decode(int16_t h[END_N], const uint8_t in[END_PKBYTES])
{
    uint64_t tmp[PK_NW86];
    uint32_t tmp32[PK_NW3 + 1];
    uint16_t r86[END_N], r3[END_N];
    uint32_t t;
    size_t i, j, k, idx;

    memset(tmp, 0, sizeof tmp);
    idx = PK_HB_OFF;
    for (i = 0; i < PK_NLOW; i += 24) {
        uint64_t w21 = 0;
        for (j = 0; j < 7; j++) {
            tmp[i + 3 * j] = ((uint64_t)(in[idx + j] >> 1)) << 38;
            w21 |= ((uint64_t)(in[idx + j] & 1u)) << (38 + j);
        }
        tmp[i + 21] = w21;
        idx += 7;
    }
    j = 0;
    for (i = 0; i < PK_NLOW; i += 3) {
        uint64_t r0 = load64_le(in + 8 * j), r1 = load64_le(in + 8 * (j + 1));
        tmp[i + 0] += ((r1 & 0x7FFFFULL) << 19) | (r0 & 0x7FFFFULL);
        tmp[i + 1] = r0 >> 19;
        tmp[i + 2] = r1 >> 19;
        j += 2;
    }
#if END_N == 512
    tmp[72] = 0;
    for (i = 0; i < 6; i++) tmp[72] |= ((uint64_t)in[405 + i]) << (8 * i);
    r86[END_N - 1] = in[411];
    for (i = 0; i < 73; i++) {
        uint64_t x = tmp[i];
        for (k = 0; k < 7; k++) { uint64_t q; uint16_t r; ct_divmod_small_u64(x, 86, MFLOOR86, &q, &r); r86[7 * i + k] = r; x = q; }
    }
#else
    tmp[144] = 0; tmp[145] = 0;
    for (i = 0; i < 5; i++) { tmp[144] |= ((uint64_t)in[810 + i]) << (8 * i); tmp[145] |= ((uint64_t)in[815 + i]) << (8 * i); }
    tmp[144] |= ((uint64_t)((in[820] >> 3) & 31u)) << 40;
    tmp[145] |= ((uint64_t)((in[821] >> 3) & 31u)) << 40;
    tmp[146] = (uint64_t)(in[820] & 7u) | ((uint64_t)(in[821] & 7u) << 3) | ((uint64_t)in[822] << 6);
    for (i = 0; i < 146; i++) {
        uint64_t x = tmp[i];
        for (k = 0; k < 7; k++) { uint64_t q; uint16_t r; ct_divmod_small_u64(x, 86, MFLOOR86, &q, &r); r86[7 * i + k] = r; x = q; }
    }
    {
        uint64_t x = tmp[146], q; uint16_t r;
        ct_divmod_small_u64(x, 86, MFLOOR86, &q, &r); r86[END_N - 2] = r;
        ct_divmod_small_u64(q, 86, MFLOOR86, &q, &r); r86[END_N - 1] = r;
    }
#endif
    /* base-3 part */
#if END_N == 512
    tmp32[30] = in[513] >> 2;
    t = (((uint32_t)in[513] & 3u) << 8) | (uint32_t)in[512];
#else
    t = load32_le(in + 1023);
    tmp32[60] = t >> 20;
#endif
    j = 0;
    for (i = 0; i < PK_NW3; i += 6) {
        uint32_t rs[5];
        for (k = 0; k < 5; k++) rs[k] = load32_le(in + PK_M3_OFF + 4 * (j + k));
        for (k = 0; k < 5; k++) tmp32[i + k] = rs[k] & 0x07FFFFFFu;
        tmp32[i + 5] = 0;
        for (k = 0; k < 5; k++) tmp32[i + 5] |= (rs[k] >> 27) << (5 * k);
        j += 5;
    }
    {
        size_t ngroups = PK_NW3 / 6, gi;
        for (gi = 0; gi < ngroups; gi++)
            tmp32[6 * gi + 5] += ((t >> (2 * (ngroups - 1 - gi))) & 3u) << 25;
    }
    for (i = 0; i < PK_NW3; i++) {
        uint32_t x = tmp32[i];
        for (k = 0; k < 17; k++) { uint64_t q; uint16_t r; ct_divmod_small_u64(x, 3, MFLOOR3, &q, &r); r3[17 * i + k] = r; x = (uint32_t)q; }
    }
    {
        uint32_t x = tmp32[PK_NW3];
        for (k = 0; k < (END_N == 512 ? 2u : 4u); k++) { uint64_t q; uint16_t r; ct_divmod_small_u64(x, 3, MFLOOR3, &q, &r); r3[17 * PK_NW3 + k] = r; x = (uint32_t)q; }
    }
    /* CRT mod 258: x = r86 * 87 + r3 * 172 (mod 258); a value of 257 is folded to 0 by red257 */
    for (i = 0; i < END_N; i++)
        h[i] = red257((int16_t)ct_mod258_u32_small((uint32_t)r86[i] * 87u + (uint32_t)r3[i] * 172u));
}

/* ------------------------------------------------------------------------- */
/* ciphertext                                                                 */
/* ------------------------------------------------------------------------- */
#if END_N == 512
void ct_encode(uint8_t out[END_CTBYTES], const int16_t c[END_N])
{
    size_t i, j = 0;
    for (i = 0; i < END_N; i += 4) {
        out[j]     = (uint8_t)(((c[i + 3] & 3) << 6) | (c[i] & 63));
        out[j + 1] = (uint8_t)((((c[i + 3] >> 2) & 3) << 6) | (c[i + 1] & 63));
        out[j + 2] = (uint8_t)((((c[i + 3] >> 4) & 3) << 6) | (c[i + 2] & 63));
        j += 3;
    }
}
void ct_decode(int16_t c[END_N], const uint8_t in[END_CTBYTES])
{
    size_t i, j = 0;
    for (i = 0; i < END_CTBYTES; i += 3) {
        c[j] = in[i] & 0x3f; c[j + 1] = in[i + 1] & 0x3f; c[j + 2] = in[i + 2] & 0x3f;
        c[j + 3] = (int16_t)((((in[i + 2] >> 6) & 3) << 4) | (((in[i + 1] >> 6) & 3) << 2) | ((in[i] >> 6) & 3));
        j += 4;
    }
}
#else
static const uint64_t c_table90[8] = { 1ULL, 90ULL, 8100ULL, 729000ULL, 65610000ULL, 5904900000ULL, 531441000000ULL, 47829690000000ULL };
void ct_encode(uint8_t out[END_CTBYTES], const int16_t c[END_N])
{
    uint64_t tmp[128], res[104];
    size_t i, index;
    for (i = 0; i < 128; i++) {
        uint64_t acc = 0; size_t k;
        for (k = 0; k < 8; k++) acc += (uint64_t)(uint16_t)c[8 * i + k] * c_table90[k];
        tmp[i] = acc;
    }
    memset(res, 0, sizeof res);
    index = 0;
    for (i = 0; i < 80; i += 5) {
        res[index + 0] = (tmp[i + 0] << 12) | (tmp[i + 4] & 0xFFFULL);
        res[index + 1] = (tmp[i + 1] << 12) | ((tmp[i + 4] >> 12) & 0xFFFULL);
        res[index + 2] = (tmp[i + 2] << 12) | ((tmp[i + 4] >> 24) & 0xFFFULL);
        res[index + 3] = (tmp[i + 3] << 12) | ((tmp[i + 4] >> 36) & 0xFFFULL);
        index += 4;
        res[64] = (res[64] << 4) | (tmp[i + 4] >> 48);
    }
    index = 65;
    for (i = 80; i < 125; i += 5) {
        res[index + 0] = (tmp[i + 0] << 12) | (tmp[i + 4] & 0xFFFULL);
        res[index + 1] = (tmp[i + 1] << 12) | ((tmp[i + 4] >> 12) & 0xFFFULL);
        res[index + 2] = (tmp[i + 2] << 12) | ((tmp[i + 4] >> 24) & 0xFFFULL);
        res[index + 3] = (tmp[i + 3] << 12) | ((tmp[i + 4] >> 36) & 0xFFFULL);
        index += 4;
        res[101] = (res[101] << 4) | (tmp[i + 4] >> 48);
    }
    res[101] = (res[101] << 28) | (tmp[125] & 0x0FFFFFFFULL);
    res[102] = (tmp[126] << 12) | ((tmp[125] >> 28) & 0xFFFULL);
    res[103] = (tmp[127] << 12) | ((tmp[125] >> 40) & 0xFFFULL);
    for (i = 0; i < 104; i++) store64_le(out + 8 * i, res[i]);
}
void ct_decode(int16_t c[END_N], const uint8_t in[END_CTBYTES])
{
    uint64_t res[104], tmp[128], t;
    size_t i, j, index;
    for (i = 0; i < 104; i++) res[i] = load64_le(in + 8 * i);
    tmp[125] = ((res[103] & 0xFFFULL) << 40) | ((res[102] & 0xFFFULL) << 28) | (res[101] & 0x0FFFFFFFULL);
    tmp[126] = (res[102] >> 12) & 0xFFFFFFFFFFFFFULL;
    tmp[127] = (res[103] >> 12) & 0xFFFFFFFFFFFFFULL;
    index = 65; t = res[101] >> 28; j = 0;
    for (i = 80; i < 125; i += 5) {
        tmp[i + 0] = res[index + 0] >> 12; tmp[i + 1] = res[index + 1] >> 12;
        tmp[i + 2] = res[index + 2] >> 12; tmp[i + 3] = res[index + 3] >> 12;
        tmp[i + 4] = (((t >> (32 - j * 4)) & 15ULL) << 48) | ((res[index + 3] & 0xFFFULL) << 36) |
                     ((res[index + 2] & 0xFFFULL) << 24) | ((res[index + 1] & 0xFFFULL) << 12) | (res[index + 0] & 0xFFFULL);
        index += 4; j++;
    }
    index = 0; t = res[64]; j = 0;
    for (i = 0; i < 80; i += 5) {
        tmp[i + 0] = res[index + 0] >> 12; tmp[i + 1] = res[index + 1] >> 12;
        tmp[i + 2] = res[index + 2] >> 12; tmp[i + 3] = res[index + 3] >> 12;
        tmp[i + 4] = (((t >> (60 - j * 4)) & 15ULL) << 48) | ((res[index + 3] & 0xFFFULL) << 36) |
                     ((res[index + 2] & 0xFFFULL) << 24) | ((res[index + 1] & 0xFFFULL) << 12) | (res[index + 0] & 0xFFFULL);
        index += 4; j++;
    }
    for (i = 0; i < 128; i++) {
        uint64_t x = tmp[i];
        for (j = 0; j < 8; j++) { uint64_t q; uint16_t r; ct_divmod_small_u64(x, 90, MFLOOR90, &q, &r); c[8 * i + j] = (int16_t)r; x = q; }
    }
}
#endif

/* ------------------------------------------------------------------------- */
/* message layout                                                             */
/* ------------------------------------------------------------------------- */
void s_from_mr(int16_t s[END_N], const uint8_t mprime[END_MSGBYTES], const uint8_t r[END_MSGBYTES])
{
    size_t i, j;
    for (i = 0; i < 8; i++)
        for (j = 0; j < END_MSGBYTES; j++) {
            s[i * END_MSGBYTES + j] = (int16_t)((mprime[j] >> i) & 1);
            s[END_N / 2 + i * END_MSGBYTES + j] = (int16_t)((r[j] >> i) & 1);
        }
}

void mr_from_s(uint8_t mprime[END_MSGBYTES], uint8_t r[END_MSGBYTES], const int16_t s[END_N])
{
    size_t i, j;
    memset(mprime, 0, END_MSGBYTES); memset(r, 0, END_MSGBYTES);
    for (i = 0; i < 8; i++)
        for (j = 0; j < END_MSGBYTES; j++) {
            mprime[j] |= (uint8_t)((s[i * END_MSGBYTES + j] & 1) << i);
            r[j] |= (uint8_t)((s[END_N / 2 + i * END_MSGBYTES + j] & 1) << i);
        }
}
