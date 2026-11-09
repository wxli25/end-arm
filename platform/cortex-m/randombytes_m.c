/*
 * randombytes() for the bare-metal Cortex-M builds.
 *
 * !!! TEST-ONLY DRBG: SHAKE256(seed || counter) with a fixed seed. Replace platform_entropy() with the
 * board's TRNG (e.g. STM32 RNG peripheral) before any real use. The seed can also be overridden from
 * the host through semihosting arguments in the bench/test harness.
 */
#include <stdint.h>
#include <string.h>
#include "end_api.h"
#include "fips202.h"

static uint8_t drbg_seed[32] = { 0x45, 0x4e, 0x44, 0x2d, 0x41, 0x52, 0x4d, 0x2d, 0x54, 0x45, 0x53, 0x54, 0x2d, 0x4f, 0x4e, 0x4c, 0x59 };
static uint32_t drbg_ctr = 0;

__attribute__((weak)) int platform_entropy(uint8_t *out, size_t len)   /* return 0 if no TRNG available */
{
    (void)out; (void)len;
    return 0;
}

void END_NS(randombytes)(uint8_t *out, size_t outlen)
{
    uint8_t in[36];
    if (platform_entropy(out, outlen))
        return;
    memcpy(in, drbg_seed, 32);
    in[32] = (uint8_t)drbg_ctr; in[33] = (uint8_t)(drbg_ctr >> 8); in[34] = (uint8_t)(drbg_ctr >> 16); in[35] = (uint8_t)(drbg_ctr >> 24);
    drbg_ctr++;
    shake256(out, outlen, in, 36);
}
