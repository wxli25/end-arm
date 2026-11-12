/*
 * END KEM -- glue between the namespaced END API and pqm4's harness.
 *
 * pqm4 links one scheme per binary and calls the un-namespaced crypto_kem_* functions
 * declared in api.h. Randomness comes from pqm4's randombytes() (the board's TRNG, or the
 * deterministic generator of the testvectors build).
 */
#include "api.h"
#include "end_api.h"

#if CRYPTO_PUBLICKEYBYTES != END_PKBYTES || CRYPTO_SECRETKEYBYTES != END_SKBYTES || \
    CRYPTO_CIPHERTEXTBYTES != END_CTBYTES || CRYPTO_BYTES != END_SSBYTES
#error "api.h sizes do not match end_params.h"
#endif

static int pqm4_randombytes(uint8_t *out, size_t outlen);

/* Defined before randombytes.h is included: that header turns "randombytes" into a macro,
 * which would otherwise rename END_NS(randombytes). */
int END_NS(randombytes)(uint8_t *out, size_t outlen)
{
    return pqm4_randombytes(out, outlen);
}

int crypto_kem_keypair(uint8_t *pk, uint8_t *sk)
{
    return END_NS(crypto_kem_keypair)(pk, sk);
}

int crypto_kem_enc(uint8_t *ct, uint8_t *ss, const uint8_t *pk)
{
    return END_NS(crypto_kem_enc)(ct, ss, pk);
}

int crypto_kem_dec(uint8_t *ss, const uint8_t *ct, const uint8_t *sk)
{
    return END_NS(crypto_kem_dec)(ss, ct, sk);
}

#include "randombytes.h"

static int pqm4_randombytes(uint8_t *out, size_t outlen)
{
    return randombytes(out, outlen) == 0 ? 0 : -1;
}
