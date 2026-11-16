/*
 * Cumulative signature workload for the Cortex-M4 cycle model (pqm4/m4model), the signature
 * counterpart of pqm4/icount/icount.h.  The file names mirror the KEM harness so that
 * run_m4model.sh and m4model_report.py treat both alike:
 *   icount_keypair.c : keypair            icount_kg10.c : ten key generations
 *   icount_enc.c     : keypair + sign     icount_dec.c  : keypair + sign + verify (open)
 * so keygen = I(keypair) - I(none), sign = I(enc) - I(keypair), verify = I(dec) - I(enc).
 * The message is 59 bytes from randombytes, as in pqm4's mupq/crypto_sign/speed.c.
 */
#include "api.h"
#include "hal.h"
#include "randombytes.h"

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#define MLEN 59

#define PASTER(x, y) x##y
#define EVALUATOR(x, y) PASTER(x, y)
#define NAMESPACE(fun) EVALUATOR(MUPQ_NAMESPACE, fun)

#define MUPQ_CRYPTO_PUBLICKEYBYTES NAMESPACE(CRYPTO_PUBLICKEYBYTES)
#define MUPQ_CRYPTO_SECRETKEYBYTES NAMESPACE(CRYPTO_SECRETKEYBYTES)
#define MUPQ_CRYPTO_BYTES          NAMESPACE(CRYPTO_BYTES)
#define MUPQ_crypto_sign_keypair   NAMESPACE(crypto_sign_keypair)
#define MUPQ_crypto_sign           NAMESPACE(crypto_sign)
#define MUPQ_crypto_sign_open      NAMESPACE(crypto_sign_open)

#ifndef ICOUNT_KEYPAIR_REPS
#define ICOUNT_KEYPAIR_REPS 1
#endif

int main(void)
{
    unsigned char sk[MUPQ_CRYPTO_SECRETKEYBYTES];
    unsigned char pk[MUPQ_CRYPTO_PUBLICKEYBYTES];
    unsigned char m[MLEN];
    unsigned char sm[MLEN + MUPQ_CRYPTO_BYTES];
    size_t smlen = 0;
    int rc = 0;

    hal_setup(CLOCK_BENCHMARK);
#if ICOUNT_OPS >= 1
    for (int rep = 0; rep < ICOUNT_KEYPAIR_REPS; rep++)
        MUPQ_crypto_sign_keypair(pk, sk);
#endif
#if ICOUNT_OPS >= 2
    randombytes(m, MLEN);
    MUPQ_crypto_sign(sm, &smlen, m, MLEN, sk);
#endif
#if ICOUNT_OPS >= 3
    rc = MUPQ_crypto_sign_open(sm, &smlen, sm, smlen, pk);
#endif
    hal_send_str(rc == 0 ? "icount OK" : "icount VERIFY FAILED");
    hal_send_str("#");
    (void)pk; (void)sk; (void)m; (void)sm; (void)smlen;
    return 0;
}
