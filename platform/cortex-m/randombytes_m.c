/*
 * randombytes() for the bare-metal Cortex-M builds.
 *
 * Board builds must supply platform_entropy() backed by a TRNG or seeded DRBG.
 * The fixed-seed stream is compiled only into the QEMU test build.
 */
#include <stdint.h>
#include <string.h>
#include "end_api.h"
#ifdef END_TEST_DRBG
#include "fips202.h"

static const uint8_t drbg_seed[32] = { 0x45, 0x4e, 0x44, 0x2d, 0x41, 0x52, 0x4d, 0x2d, 0x54, 0x45, 0x53, 0x54, 0x2d, 0x4f, 0x4e, 0x4c, 0x59 };
static uint32_t drbg_ctr = 0;
#endif

/* Return 1 after filling the entire buffer. Return 0 on failure. */
__attribute__((weak)) int platform_entropy(uint8_t *out, size_t len)
{
    (void)out; (void)len;
    return 0;
}

int END_NS(randombytes)(uint8_t *out, size_t outlen)
{
    if (platform_entropy(out, outlen) == 1)
        return 0;
#ifdef END_TEST_DRBG
    uint8_t in[36];
    if (drbg_ctr == UINT32_MAX)
        return -1;
    memcpy(in, drbg_seed, 32);
    in[32] = (uint8_t)drbg_ctr; in[33] = (uint8_t)(drbg_ctr >> 8); in[34] = (uint8_t)(drbg_ctr >> 16); in[35] = (uint8_t)(drbg_ctr >> 24);
    drbg_ctr++;
    shake256(out, outlen, in, 36);
    return 0;
#else
    return -1;
#endif
}
