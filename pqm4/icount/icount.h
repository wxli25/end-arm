/*
 * Cumulative KEM workload for dynamic instruction counting under QEMU (TCG plugin libinsn).
 *
 * pqm4 builds mupq/crypto_kem/<name>.c into elf/<impl>_<name>.elf for every scheme, so the four
 * files icount_{none,keypair,enc,dec}.c (ICOUNT_OPS = 0..3) give, for any KEM in pqm4:
 *   keypair = I(keypair) - I(none),  encaps = I(enc) - I(keypair),  decaps = I(dec) - I(enc)
 * measured with pqm4's own harness, HAL, randombytes and Keccak assembly.
 * icount_kg10.c runs ten key generations (ICOUNT_KEYPAIR_REPS) for an average over seeds:
 *   keypair(avg) = (I(kg10) - I(none)) / 10.
 */
#include "api.h"
#include "hal.h"

#include <stdint.h>
#include <string.h>

#define PASTER(x, y) x##y
#define EVALUATOR(x, y) PASTER(x, y)
#define NAMESPACE(fun) EVALUATOR(MUPQ_NAMESPACE, fun)

#define MUPQ_CRYPTO_BYTES           NAMESPACE(CRYPTO_BYTES)
#define MUPQ_CRYPTO_PUBLICKEYBYTES  NAMESPACE(CRYPTO_PUBLICKEYBYTES)
#define MUPQ_CRYPTO_SECRETKEYBYTES  NAMESPACE(CRYPTO_SECRETKEYBYTES)
#define MUPQ_CRYPTO_CIPHERTEXTBYTES NAMESPACE(CRYPTO_CIPHERTEXTBYTES)
#define MUPQ_crypto_kem_keypair NAMESPACE(crypto_kem_keypair)
#define MUPQ_crypto_kem_enc NAMESPACE(crypto_kem_enc)
#define MUPQ_crypto_kem_dec NAMESPACE(crypto_kem_dec)

int main(void)
{
    unsigned char key_a[MUPQ_CRYPTO_BYTES], key_b[MUPQ_CRYPTO_BYTES];
    unsigned char sk[MUPQ_CRYPTO_SECRETKEYBYTES];
    unsigned char pk[MUPQ_CRYPTO_PUBLICKEYBYTES];
    unsigned char ct[MUPQ_CRYPTO_CIPHERTEXTBYTES];

    hal_setup(CLOCK_BENCHMARK);
    memset(key_a, 0, sizeof key_a);
    memset(key_b, 0, sizeof key_b);
#ifndef ICOUNT_KEYPAIR_REPS
#define ICOUNT_KEYPAIR_REPS 1
#endif
#if ICOUNT_OPS >= 1
    for (int rep = 0; rep < ICOUNT_KEYPAIR_REPS; rep++)  /* averages variable-time key generation */
        MUPQ_crypto_kem_keypair(pk, sk);
#endif
#if ICOUNT_OPS >= 2
    MUPQ_crypto_kem_enc(ct, key_a, pk);
#endif
#if ICOUNT_OPS >= 3
    MUPQ_crypto_kem_dec(key_b, ct, sk);
#endif
#if ICOUNT_OPS >= 3
    hal_send_str(memcmp(key_a, key_b, sizeof key_a) == 0 ? "icount OK" : "icount KEY MISMATCH");
#else
    hal_send_str("icount OK");
#endif
    hal_send_str("#");
    (void)pk; (void)sk; (void)ct;
    return 0;
}
