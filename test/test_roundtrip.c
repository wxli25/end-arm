/*
 * Random round trips: keypair -> enc -> dec, plus implicit-rejection checks on tampered ciphertexts
 * (the rejected shared secret must equal SHA3-256(k || ct')) and sk/expanded-key consistency.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "end_api.h"
#include "fips202.h"

int main(int argc, char **argv)
{
    int n = argc > 1 ? atoi(argv[1]) : 20;
    static uint8_t pk[END_PKBYTES], sk[END_SKBYTES], ct[END_CTBYTES], ct2[END_CTBYTES];
    static uint8_t ss[END_SSBYTES], ss2[END_SSBYTES], ss3[END_SSBYTES], buf[END_MSGBYTES + END_CTBYTES];
    static end_expkey ek;
    int i, fails = 0;
    if (n < 1) {
        fprintf(stderr, "Expected a positive iteration count\n");
        return 2;
    }
    printf("test_roundtrip %s: %d iterations\n", END_NAME, n);
    for (i = 0; i < n; i++) {
        if (END_NS(crypto_kem_keypair)(pk, sk)) { printf("keypair failed\n"); return 1; }
        if (END_NS(crypto_kem_enc)(ct, ss, pk)) { printf("encaps failed\n"); return 1; }
        END_NS(crypto_kem_dec)(ss2, ct, sk);
        if (memcmp(ss, ss2, END_SSBYTES)) { printf("iter %d: shared secret mismatch\n", i); fails++; }
        END_NS(expkey_from_sk)(&ek, sk);
        END_NS(crypto_kem_dec_exp)(ss3, ct, &ek);
        if (memcmp(ss, ss3, END_SSBYTES)) { printf("iter %d: dec_exp mismatch\n", i); fails++; }
        /* tamper: flip one bit -> implicit rejection with K = SHA3-256(k || ct') */
        memcpy(ct2, ct, END_CTBYTES);
        ct2[i % END_CTBYTES] ^= (uint8_t)(1u << (i % 8));
        END_NS(crypto_kem_dec)(ss2, ct2, sk);
        memcpy(buf, sk + END_SK_OFF_K, END_MSGBYTES);
        memcpy(buf + END_MSGBYTES, ct2, END_CTBYTES);
        sha3_256(ss3, buf, sizeof buf);
        if (!memcmp(ss, ss2, END_SSBYTES)) { printf("iter %d: tampered ct accepted\n", i); fails++; }
        if (memcmp(ss2, ss3, END_SSBYTES)) { printf("iter %d: rejection key != H(k||ct')\n", i); fails++; }
    }
    printf("test_roundtrip %s: %d failures -> %s\n", END_NAME, fails, fails ? "FAIL" : "ALL OK");
    return fails != 0;
}
