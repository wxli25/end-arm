/* Differential test for the Cortex-M decoder, including total decoding of
 * noncanonical public keys and callers supplying unaligned byte buffers. */
#include <stdio.h>
#include <string.h>
#include "end_internal.h"

void END_NS(pk_decode_portable)(int16_t out[END_N], const uint8_t in[END_PKBYTES]);

static uint32_t rng_state = 0xb93d1e7au;

static uint32_t random_word(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return rng_state;
}

static int check(const uint8_t *bytes, const int16_t *expected, unsigned trial)
{
    static int16_t reference[END_N], actual[END_N];
    unsigned i;
    END_NS(pk_decode_portable)(reference, bytes);
    pk_decode(actual, bytes);
    for (i = 0; i < END_N; i++) {
        if (actual[i] != reference[i] || (expected && actual[i] != expected[i])) {
            printf("pk_decode case %u coefficient %u: %d (ref %d)\n",
                   trial, i, actual[i], reference[i]);
            return 1;
        }
    }
    return 0;
}

int main(void)
{
    static uint8_t storage[END_PKBYTES + 3];
    static int16_t poly[END_N];
    unsigned i, trial;
    for (trial = 0; trial < 32; trial++) {
        uint8_t *pk = storage + (trial & 3);
        for (i = 0; i < END_N; i++) poly[i] = (int16_t)(random_word() % 257);
        pk_encode(pk, poly);
        if (check(pk, poly, trial)) return 1;
    }
    for (trial = 0; trial < 256; trial++) {
        uint8_t *pk = storage + (trial & 3);
        for (i = 0; i < END_PKBYTES; i++) pk[i] = (uint8_t)random_word();
#if END_N == 512
        pk[411] = (uint8_t)trial;           /* every possible raw last digit */
#else
        pk[820] = (uint8_t)trial;
        pk[821] = (uint8_t)~trial;
        pk[822] = (uint8_t)(trial * 19u);
#endif
        if (check(pk, NULL, 32 + trial)) return 1;
    }
    memset(storage, 0, sizeof storage);
    if (check(storage, NULL, 288)) return 1;
    memset(storage, 0xff, sizeof storage);
    if (check(storage, NULL, 289)) return 1;
    printf("END-%d Cortex-M pk_decode: 290 differential cases passed\n", END_N);
    return 0;
}
