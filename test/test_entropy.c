/* Exercise the board RNG adapter on the host without enabling its test DRBG. */
#include <stdio.h>
#include <string.h>
#include "end_api.h"

static int entropy_ready;

int platform_entropy(uint8_t *out, size_t len)
{
    memset(out, 0x42, len); /* Also simulate an entropy source that fails after writing. */
    return entropy_ready;
}

static int untouched(const uint8_t *buf, size_t len)
{
    size_t i;
    for (i = 0; i < len; i++)
        if (buf[i] != 0xa5)
            return 0;
    return 1;
}

int main(void)
{
    static uint8_t pk[END_PKBYTES], sk[END_SKBYTES], ct[END_CTBYTES];
    uint8_t ss[END_SSBYTES], recovered[END_SSBYTES];
    memset(pk, 0xa5, sizeof pk);
    memset(sk, 0xa5, sizeof sk);
    memset(ct, 0xa5, sizeof ct);
    memset(ss, 0xa5, sizeof ss);

    if (END_NS(crypto_kem_keypair)(pk, sk) != -1 ||
        END_NS(crypto_kem_enc)(ct, ss, pk) != -1 ||
        !untouched(pk, sizeof pk) || !untouched(sk, sizeof sk) ||
        !untouched(ct, sizeof ct) || !untouched(ss, sizeof ss)) {
        fprintf(stderr, "Board entropy failure did not leave KEM outputs untouched\n");
        return 1;
    }

    entropy_ready = 1;
    if (END_NS(crypto_kem_keypair)(pk, sk) ||
        END_NS(crypto_kem_enc)(ct, ss, pk) ||
        END_NS(crypto_kem_dec)(recovered, ct, sk) ||
        memcmp(ss, recovered, sizeof ss) != 0) {
        fprintf(stderr, "Board entropy success did not complete a KEM round trip\n");
        return 1;
    }
    puts("Board entropy adapter: failure and success paths passed");
    return 0;
}
