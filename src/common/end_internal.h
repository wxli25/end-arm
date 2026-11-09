/*
 * END KEM -- internal interfaces shared by the scheme-level modules (L3-L7).
 */
#ifndef END_INTERNAL_H
#define END_INTERNAL_H

#include <stdint.h>
#include <stddef.h>
#include "end_params.h"
#include "end_api.h"
#include "poly.h"
#include "fips202.h"

/* ---- codec (L3) ----------------------------------------------------- */
void END_NS(pk_encode)(uint8_t out[END_PKBYTES], const int16_t h_ntt[END_N]);  /* h in [0,257), NTT-257 layout */
void END_NS(pk_decode)(int16_t h_ntt[END_N], const uint8_t in[END_PKBYTES]);    /* out in [0,257)               */
void END_NS(ct_encode)(uint8_t out[END_CTBYTES], const int16_t c[END_N]);       /* c in [0, END_DE)             */
void END_NS(ct_decode)(int16_t c[END_N], const uint8_t in[END_CTBYTES]);         /* out in [0, END_DE)           */
void END_NS(s_from_mr)(int16_t s[END_N], const uint8_t mprime[END_MSGBYTES], const uint8_t r[END_MSGBYTES]);
void END_NS(mr_from_s)(uint8_t mprime[END_MSGBYTES], uint8_t r[END_MSGBYTES], const int16_t s[END_N]);
#define pk_encode END_NS(pk_encode)
#define pk_decode END_NS(pk_decode)
#define ct_encode END_NS(ct_encode)
#define ct_decode END_NS(ct_decode)
#define s_from_mr END_NS(s_from_mr)
#define mr_from_s END_NS(mr_from_s)

static inline int16_t end_compress(int16_t x)      /* x in [0,256] -> [0, END_DE) */
{
#if END_N == 512
    return (int16_t)(((uint32_t)x * 255u) >> 10);
#else
    return (int16_t)(((uint32_t)x * 717u) >> 11);
#endif
}
static inline int16_t end_decompress(int16_t c)    /* c in [0, END_DE) -> [0,256] */
{
#if END_N == 512
    return (int16_t)(4 * c + 3);
#else
    return (int16_t)((((uint32_t)c * 731u) >> 8) + 2);
#endif
}

/* ---- hashes ----------------------------------------------------------- */
static inline void end_hash1(uint8_t out[END_HASH1_BYTES], const uint8_t *in, size_t inlen)
{
#if END_N == 512
    sha3_256(out, in, inlen);
#else
    sha3_512(out, in, inlen);
#endif
}

/* ---- trapdoor decoder (L6) ------------------------------------------ */
/*
 * Alg. 6 (END.DecodeTrap) in the doubled domain with the corrected signs.
 * cprime: decompressed ciphertext coefficients in [0,256].
 * sbar  : receives s~ (entries +-1 on success).
 * returns 1 if a GF-stage candidate passed the distributional (binary) check, 0 otherwise.
 * diag (may be NULL): [0] = gf candidate index chosen, [1] = GF candidate index chosen,
 *                      [2..3] = gf FCL indices, [4..5] = GF FCL indices.
 */
int END_NS(decode_trap)(int16_t sbar[END_N], const end_expkey *ek, const int16_t cprime[END_N], uint32_t diag[8]);
#define decode_trap END_NS(decode_trap)

/* ---- key material --------------------------------------------------- */
void END_NS(expkey_from_parts)(end_expkey *ek, const uint8_t k[END_MSGBYTES], const int16_t h_ntt[END_N],
                               const int8_t f[END_N], const int8_t g[END_N], const int8_t F[END_N],
                               const int8_t G[END_N], const int16_t w[END_N]);
#define expkey_from_parts END_NS(expkey_from_parts)

/* ---- DPKE / PKE (L7) ------------------------------------------------- */
void END_NS(dpke_encrypt)(int16_t c[END_N], const int16_t h_ntt[END_N], const int16_t s[END_N]);
void END_NS(pke_encrypt)(uint8_t ct[END_CTBYTES], const int16_t h_ntt[END_N], const uint8_t m[END_MSGBYTES],
                         const uint8_t r[END_MSGBYTES]);
int  END_NS(pke_decrypt)(uint8_t m[END_MSGBYTES], const end_expkey *ek, const uint8_t ct[END_CTBYTES]);
#define dpke_encrypt END_NS(dpke_encrypt)
#define pke_encrypt  END_NS(pke_encrypt)
#define pke_decrypt  END_NS(pke_decrypt)

/* ---- key generation internals (L4) ----------------------------------- */
void END_NS(ternary_sample)(int8_t a[END_N], keccak_state *state);
#define ternary_sample END_NS(ternary_sample)

#endif /* END_INTERNAL_H */
