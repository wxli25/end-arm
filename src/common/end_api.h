/*
 * END KEM -- public API (L8).
 *
 * Two parameter sets are selected at compile time (-DEND_N=512 | 1024); all symbols are prefixed
 * end512_ / end1024_ so that both can be linked into one binary.
 *
 *   keypair / keypair_derand  : NTRU trapdoor generation (slow: ~10^7 cycles), needs a scratch buffer
 *   enc / enc_derand          : KEM encapsulation (FO transform, implicit rejection)
 *   dec                       : decapsulation from the compact secret key
 *   expkey_from_sk + dec_exp  : decapsulation from a pre-expanded key (saves ~15 NTTs per call)
 *
 * All functions return 0 on success. Randomized key generation and encapsulation
 * return -1 on entropy failure, leaving their output buffers untouched.
 * Decapsulation uses implicit rejection.
 */
#ifndef END_API_H
#define END_API_H

#include <stdint.h>
#include <stddef.h>
#include "end_params.h"

#define END_CRYPTO_PUBLICKEYBYTES  END_PKBYTES
#define END_CRYPTO_SECRETKEYBYTES  END_SKBYTES
#define END_CRYPTO_CIPHERTEXTBYTES END_CTBYTES
#define END_CRYPTO_BYTES           END_SSBYTES

/* expanded secret key: everything decapsulation needs, with the NTT forms precomputed */
typedef struct END_ALIGN {
    int16_t h_ntt[END_N];          /* public key h, NTT-257 layout, [0,257)                     */
    uint8_t k[END_MSGBYTES];       /* implicit-rejection secret                                 */
    int16_t f16[END_N];            /* f (ternary) as int16                                      */
    int16_t F16[END_N];            /* F (|.| <= 31) as int16                                    */
    int16_t w[END_N];              /* w, |w| <= 2^15                                            */
    int16_t fn257[END_N];          /* NTT257(f)                                                 */
    int16_t Fn257[END_N];          /* NTT257(F)                                                 */
    int16_t wfn257[END_N];         /* NTT257(w*f)                                               */
    int16_t const1[END_N];         /* (g + delta f) * 1  mod 257                                */
    int16_t const2[END_N];         /* (G + delta F) * 1  mod 257                                */
    int16_t wconst[END_N];         /* w * (g + delta f) * 1  mod 257                            */
    int16_t wn12289[END_N];        /* NTT12289(w)                                               */
    int16_t fn3329[END_N];         /* NTT3329(f)                                                */
    int16_t FQn3329[END_N];        /* NTT3329(F_Q), F_Q = Q F - w f                             */
    uint16_t p1;                   /* parity of every coefficient of c1 = g s~ + f mu~          */
    uint16_t p2;                   /* parity of every coefficient of G s~ + F mu~               */
    uint16_t pwc1;                 /* parity of every coefficient of w*c1  (= p1 & parity(sum w))*/
} end_expkey;

int END_NS(crypto_kem_keypair)(uint8_t *pk, uint8_t *sk);
int END_NS(crypto_kem_keypair_derand)(uint8_t *pk, uint8_t *sk, const uint8_t seed[END_SEEDBYTES],
                                      void *tmp, size_t tmplen);   /* tmp: >= END_KEYGEN_TMPBYTES bytes */
int END_NS(crypto_kem_enc)(uint8_t *ct, uint8_t *ss, const uint8_t *pk);
int END_NS(crypto_kem_enc_derand)(uint8_t *ct, uint8_t *ss, const uint8_t *pk, const uint8_t m[END_MSGBYTES]);
int END_NS(crypto_kem_dec)(uint8_t *ss, const uint8_t *ct, const uint8_t *sk);
int END_NS(expkey_from_sk)(end_expkey *ek, const uint8_t *sk);
int END_NS(crypto_kem_dec_exp)(uint8_t *ss, const uint8_t *ct, const end_expkey *ek);

/* secret-key field accessors (compact format: k || pk || f || g || F || G || w) */
#define END_SK_K(sk)  ((sk) + END_SK_OFF_K)
#define END_SK_PK(sk) ((sk) + END_SK_OFF_PK)
#define END_SK_F(sk)  ((const int8_t *)((sk) + END_SK_OFF_F))
#define END_SK_G(sk)  ((const int8_t *)((sk) + END_SK_OFF_G))
#define END_SK_FF(sk) ((const int8_t *)((sk) + END_SK_OFF_FF))
#define END_SK_GG(sk) ((const int8_t *)((sk) + END_SK_OFF_GG))

/* Platform entropy provider. Return 0 after filling out, or -1 on failure. */
int END_NS(randombytes)(uint8_t *out, size_t outlen);

#endif /* END_API_H */
