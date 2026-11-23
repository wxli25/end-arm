/* Cortex-M public-key decoder. Decode the base-86 words directly into the
 * destination; the base-3 stream then completes those same coefficients. */
#include "end_internal.h"
#include "end_arith.h"

#ifdef END_M_PORTABLE_PK_DECODE
void END_NS(pk_decode_portable)(int16_t h[END_N], const uint8_t in[END_PKBYTES]);

void pk_decode(int16_t h[END_N], const uint8_t in[END_PKBYTES])
{
    END_NS(pk_decode_portable)(h, in);
}
#else

#define LOW_WORDS (END_N == 512 ? 72u : 144u)
#define TRIT_OFFSET (END_N == 512 ? 412u : 823u)
#define TRIT_BLOCKS (END_N == 512 ? 5u : 10u)

static inline uint32_t load32le(const uint8_t *p)
{
    return (uint32_t)p[0] | ((uint32_t)p[1] << 8)
         | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}

static inline uint32_t mulhi(uint32_t a, uint32_t b)
{
    return (uint32_t)(((uint64_t)a * b) >> 32);
}

/* d*ceil(2^32/d) - 2^32 is 70 for d=86 and 2 for d=3. Thus these
 * quotients are exact for x < 2^32/70 and x < 2^31 respectively. */
static inline uint32_t div86(uint32_t x) { return mulhi(x, 49941481u); }
static inline uint32_t div3(uint32_t x) { return mulhi(x, 1431655766u); }

/* Split a word below 2^48 at 86^4. Its two halves are below 86^4
 * and 2^48/86^4, so div86 covers both, including noncanonical words. */
static void unpack86(int16_t *dst, uint32_t low, uint32_t high)
{
    const uint32_t power = 54700816u;
    uint32_t b = mulhi((high << 16) | (low >> 16), 5145718u);
    uint32_t a = low - b * power;
    unsigned i;
    /* The estimate is at most two short. a is the true remainder before
     * correction modulo 2^32 and is below 3*power. */
    for (i = 0; i < 2; i++) {
        uint32_t ge = 1u ^ ((a - power) >> 31);
        a -= (0u - ge) & power;
        b += ge;
    }
    for (i = 0; i < 3; i++) {
        uint32_t qa = div86(a), qb = div86(b);
        dst[i] = (int16_t)(a - qa * 86u);
        dst[4 + i] = (int16_t)(b - qb * 86u);
        a = qa;
        b = qb;
    }
    dst[3] = (int16_t)a;
}

/* For each triple of 45-bit words, the two low-section records contain
 * 19 low bits of the first word each and the other two words in full. */
static void unpack_triple(int16_t *dst, const uint8_t *p, uint32_t top7)
{
    uint32_t a = load32le(p), b = load32le(p + 4);
    uint32_t c = load32le(p + 8), d = load32le(p + 12);
    unpack86(dst, (a & 0x7ffffu) | (c << 19), ((c & 0x7ffffu) >> 13) | (top7 << 6));
    unpack86(dst + 7, (a >> 19) | (b << 13), b >> 19);
    unpack86(dst + 14, (c >> 19) | (d << 13), d >> 19);
}

static inline int16_t crt_digit(uint32_t a, uint32_t b)
{
    /* CRT is a + 86*((a + 2*b) mod 3). For t=3q+r <= 89,
     * (171*t mod 512) = q + 171*r; bits 7..8 therefore equal r. */
    uint32_t t = a + 2u * b;
    uint32_t r = ((171u * t) >> 7) & 3u;
    uint32_t x = a + 86u * r;
    return (int16_t)(x - (257u & (0u - ((x + 255u) >> 9))));
}

/* Interleave two independent quotient chains, reusing the two input
 * coefficient streams as the output. Packed ternary words are < 2^27. */
static void merge_trit_pair(int16_t *dst, uint32_t x, uint32_t y)
{
    unsigned i;
    for (i = 0; i < 17; i++) {
        uint32_t qx = div3(x), qy = div3(y);
        dst[i] = crt_digit((uint32_t)dst[i], x - 3u * qx);
        dst[17 + i] = crt_digit((uint32_t)dst[17 + i], y - 3u * qy);
        x = qx;
        y = qy;
    }
}

void pk_decode(int16_t h[END_N], const uint8_t in[END_PKBYTES])
{
    const uint8_t *p = in, *top = in + LOW_WORDS / 3 * 16;
    int16_t *out = h;
    uint32_t tail;
    unsigned i, j;
    for (i = 0; i < LOW_WORDS / 24; i++) {
        uint32_t last = 0;
        for (j = 0; j < 7; j++) {
            uint32_t bits = *top++;
            last |= (bits & 1u) << j;
            unpack_triple(out, p, bits >> 1);
            p += 16;
            out += 21;
        }
        unpack_triple(out, p, last);
        p += 16;
        out += 21;
    }
#if END_N == 512
    unpack86(out, load32le(top), (uint32_t)top[4] | ((uint32_t)top[5] << 8));
    out[7] = top[6];                              /* arbitrary raw byte */
    tail = (uint32_t)in[512] | ((uint32_t)in[513] << 8);
#else
    unpack86(out, load32le(top), (uint32_t)top[4] | ((uint32_t)(top[10] >> 3) << 8));
    unpack86(out + 7, load32le(top + 5), (uint32_t)top[9] | ((uint32_t)(top[11] >> 3) << 8));
    {
        uint32_t v = (top[10] & 7u) | ((uint32_t)(top[11] & 7u) << 3) | ((uint32_t)top[12] << 6);
        uint32_t q = div86(v);
        out[14] = (int16_t)(v - q * 86u);
        out[15] = (int16_t)(q - div86(q) * 86u);
    }
    tail = load32le(in + 1023);
#endif
    p = in + TRIT_OFFSET;
    out = h;
    for (i = 0; i < TRIT_BLOCKS; i++) {
        const uint32_t mask = 0x07ffffffu;
        uint32_t a = load32le(p), b = load32le(p + 4), c = load32le(p + 8);
        uint32_t d = load32le(p + 12), e = load32le(p + 16);
        uint32_t f = (a >> 27) | ((b >> 27) << 5) | ((c >> 27) << 10)
                   | ((d >> 27) << 15) | ((e >> 27) << 20);
        f |= ((tail >> (2 * (TRIT_BLOCKS - 1 - i))) & 3u) << 25;
        merge_trit_pair(out, a & mask, b & mask);
        merge_trit_pair(out + 34, c & mask, d & mask);
        merge_trit_pair(out + 68, e & mask, f);
        out += 102;
        p += 20;
    }
    tail >>= 2 * TRIT_BLOCKS;
    for (i = 0; i < END_N - TRIT_BLOCKS * 102; i++) {
        uint32_t q = div3(tail);
        uint32_t x = (uint32_t)out[i] * 87u + (tail - 3u * q) * 172u;
        /* Preserve the portable decoder's finite reduction on the raw
         * END-512 tail, including bytes above the normal digit range. */
        uint32_t d;
        for (d = 4128; d >= 258; d >>= 1)
            x -= d & (0u - (1u ^ ((x - d) >> 31)));
        out[i] = red257((int16_t)x);
        tail = q;
    }
}

#endif /* END_M_PORTABLE_PK_DECODE */
