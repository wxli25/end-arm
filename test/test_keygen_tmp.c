/* measure how much of the key-generation scratch buffer is actually written (canary scan) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "end_api.h"
int main(int argc, char **argv)
{
    static uint8_t pk[END_PKBYTES], sk[END_SKBYTES];
    size_t big = 64 * END_N;
    uint8_t *tmp = malloc(big), seed[32];
    int s, nseeds = argc > 1 ? atoi(argv[1]) : 40;
    size_t maxused = 0;
    for (s = 0; s < nseeds; s++) {
        size_t i, used = 0;
        memset(tmp, 0xA5, big);
        memset(seed, 0, 32); seed[0] = (uint8_t)s; seed[1] = (uint8_t)(s >> 8); seed[5] = 0x77;
        if (s == 0) for (i = 0; i < 32; i++) seed[i] = (uint8_t)i;     /* the bench seed */
        END_NS(crypto_kem_keypair_derand)(pk, sk, seed, tmp, big);
        for (i = 0; i < big; i++) if (tmp[i] != 0xA5) used = i + 1;
        if (used > maxused) maxused = used;
        printf("seed %2d: scratch bytes touched = %6u  (%5.2f * n)\n", s, (unsigned)used, (double)used / END_N);
    }
    printf("max touched = %u bytes = %.2f * n (configured END_KEYGEN_TMPBYTES = %d = %.2f * n)\n",
           (unsigned)maxused, (double)maxused / END_N, END_KEYGEN_TMPBYTES, (double)END_KEYGEN_TMPBYTES / END_N);
    return 0;
}
