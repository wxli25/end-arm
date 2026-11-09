/*
 * END KEM -- scheme layer (L7/L8): DPKE (Alg. 8-9), PKE = ACWC (Alg. 10-11), KEM = FO (Alg. 12-14),
 * secret-key (de)serialisation and expanded-key precomputation.
 */
#include <string.h>
#include "end_internal.h"
#include "end_arith.h"

/* ------------------------------------------------------------------------- */
/* expanded key                                                               */
/* ------------------------------------------------------------------------- */
static void ones_product_mod257(int16_t out[END_N], const int16_t an[END_N], int16_t tmp[END_N])
{
    /* out = a * (1 + x + ... + x^(n-1)) mod 257, a given in NTT-257 form */
    unsigned i;
    for (i = 0; i < END_N; i++) tmp[i] = 1;
    poly_ntt257(tmp);
    poly_basemul257(out, an, tmp);
    poly_intt257(out);
}

void expkey_from_parts(end_expkey *ek, const uint8_t k[END_MSGBYTES], const int16_t h_ntt[END_N],
                               const int8_t f[END_N], const int8_t g[END_N], const int8_t F[END_N],
                               const int8_t G[END_N], const int16_t w[END_N])
{
    int16_t t1[END_N] END_ALIGN, t2[END_N] END_ALIGN, t3[END_N] END_ALIGN;
    unsigned i;
    int32_t sw = 0, sf = 0, sg = 0, sF = 0, sG = 0;

    memcpy(ek->h_ntt, h_ntt, sizeof ek->h_ntt);
    memcpy(ek->k, k, END_MSGBYTES);
    memcpy(ek->w, w, sizeof ek->w);
    for (i = 0; i < END_N; i++) {
        ek->f16[i] = f[i]; ek->F16[i] = F[i];
        sw += w[i]; sf += f[i]; sg += g[i]; sF += F[i]; sG += G[i];
    }
    ek->p1 = (uint16_t)((sg + END_DELTA * sf) & 1);
    ek->p2 = (uint16_t)((sG + END_DELTA * sF) & 1);
    ek->pwc1 = (uint16_t)(ek->p1 & (sw & 1));

    /* NTT-257 forms */
    for (i = 0; i < END_N; i++) ek->fn257[i] = f[i];
    poly_ntt257(ek->fn257);
    for (i = 0; i < END_N; i++) ek->Fn257[i] = F[i];
    poly_ntt257(ek->Fn257);
    for (i = 0; i < END_N; i++) t1[i] = red257(w[i]);
    poly_ntt257(t1);                                                 /* NTT257(w) */
    poly_basemul257(ek->wfn257, t1, ek->fn257);                     /* NTT257(w f) */
    /* const1 = (g + delta f) * 1, const2 = (G + delta F) * 1, wconst = w * const1  (all mod 257) */
    for (i = 0; i < END_N; i++) t2[i] = red257((int16_t)(g[i] + END_DELTA * f[i]));
    poly_ntt257(t2);
    ones_product_mod257(ek->const1, t2, t3);
    poly_basemul257(t2, t2, t1);                                     /* NTT257(w (g + delta f)) */
    ones_product_mod257(ek->wconst, t2, t3);
    for (i = 0; i < END_N; i++) t2[i] = red257((int16_t)(G[i] + END_DELTA * F[i]));
    poly_ntt257(t2);
    ones_product_mod257(ek->const2, t2, t3);

    /* NTT-12289 form of w */
    memcpy(ek->wn12289, w, sizeof ek->wn12289);
    poly_ntt12289(ek->wn12289);

    /* NTT-3329 forms of f and F_Q = Q F - w f  (exact integers reduced mod 3329) */
    for (i = 0; i < END_N; i++) ek->fn3329[i] = f[i];
    poly_ntt3329(ek->fn3329);
    for (i = 0; i < END_N; i++) t1[i] = barrett3329(w[i]);
    poly_ntt3329(t1);
    poly_basemul3329(t2, t1, ek->fn3329);                            /* NTT3329(w f) */
    poly_intt3329(t2);                                               /* w f mod 3329, centred */
    for (i = 0; i < END_N; i++)
        ek->FQn3329[i] = red3329_32((int32_t)(END_QP % END_P3) * F[i] - t2[i]);   /* int32: 2302*31 overflows int16 */
    poly_ntt3329(ek->FQn3329);
#ifdef END_FULL_SCORE
    for (i = 0; i < END_N; i++) { ek->g16[i] = g[i]; ek->G16[i] = G[i]; ek->gn3329[i] = g[i]; }
    poly_ntt3329(ek->gn3329);
    poly_basemul3329(t2, t1, ek->gn3329);
    poly_intt3329(t2);
    for (i = 0; i < END_N; i++)
        ek->GQn3329[i] = red3329_32((int32_t)(END_QP % END_P3) * G[i] - t2[i]);
    poly_ntt3329(ek->GQn3329);
#endif
}

int END_NS(expkey_from_sk)(end_expkey *ek, const uint8_t *sk)
{
    int16_t h_ntt[END_N] END_ALIGN, w[END_N] END_ALIGN;
    unsigned i;
    pk_decode(h_ntt, END_SK_PK(sk));
    for (i = 0; i < END_N; i++)
        w[i] = (int16_t)((uint16_t)sk[END_SK_OFF_W + 2 * i] | ((uint16_t)sk[END_SK_OFF_W + 2 * i + 1] << 8));
    expkey_from_parts(ek, END_SK_K(sk), h_ntt, END_SK_F(sk), END_SK_G(sk), END_SK_FF(sk), END_SK_GG(sk), w);
    return 0;
}

/* ------------------------------------------------------------------------- */
/* DPKE / PKE                                                                 */
/* ------------------------------------------------------------------------- */
void dpke_encrypt(int16_t c[END_N], const int16_t h_ntt[END_N], const int16_t s[END_N])
{
    int16_t t[END_N] END_ALIGN;
    unsigned i;
    memcpy(t, s, sizeof t);
    poly_ntt257(t);
    poly_basemul257(t, h_ntt, t);
    poly_intt257(t);                       /* h s mod 257 in [0,257) */
    for (i = 0; i < END_N; i++)
        c[i] = end_compress(t[i]);
}

void pke_encrypt(uint8_t ct[END_CTBYTES], const int16_t h_ntt[END_N], const uint8_t m[END_MSGBYTES],
                         const uint8_t r[END_MSGBYTES])
{
    uint8_t mp[END_MSGBYTES];
    int16_t s[END_N] END_ALIGN, c[END_N] END_ALIGN;
    unsigned i;
    end_hash1(mp, r, END_MSGBYTES);
    for (i = 0; i < END_MSGBYTES; i++) mp[i] ^= m[i];
    s_from_mr(s, mp, r);
    dpke_encrypt(c, h_ntt, s);
    ct_encode(ct, c);
}

int pke_decrypt(uint8_t m[END_MSGBYTES], const end_expkey *ek, const uint8_t ct[END_CTBYTES])
{
    int16_t c[END_N] END_ALIGN, sbar[END_N] END_ALIGN;
    uint8_t mp[END_MSGBYTES], r[END_MSGBYTES], h[END_MSGBYTES];
    unsigned i;
    int ok;
    ct_decode(c, ct);
    for (i = 0; i < END_N; i++) c[i] = end_decompress(c[i]);
    ok = decode_trap(sbar, ek, c, NULL);
    for (i = 0; i < END_N; i++) sbar[i] = (int16_t)((sbar[i] + 1) >> 1);      /* s = (s~ + 1)/2 */
    mr_from_s(mp, r, sbar);
    end_hash1(h, r, END_MSGBYTES);
    for (i = 0; i < END_MSGBYTES; i++) m[i] = mp[i] ^ h[i];
    return ok;
}

/* ------------------------------------------------------------------------- */
/* KEM                                                                        */
/* ------------------------------------------------------------------------- */
int END_NS(crypto_kem_enc_derand)(uint8_t *ct, uint8_t *ss, const uint8_t *pk, const uint8_t m[END_MSGBYTES])
{
    int16_t h_ntt[END_N] END_ALIGN;
    uint8_t r[END_MSGBYTES];
    uint8_t buf[END_MSGBYTES + END_CTBYTES];
    pk_decode(h_ntt, pk);
    end_hash1(r, m, END_MSGBYTES);                 /* H2 = H1 (reference) */
    pke_encrypt(ct, h_ntt, m, r);
    memcpy(buf, m, END_MSGBYTES);
    memcpy(buf + END_MSGBYTES, ct, END_CTBYTES);
    sha3_256(ss, buf, sizeof buf);                 /* K = SHA3-256(m || ct) */
    return 0;
}

int END_NS(crypto_kem_enc)(uint8_t *ct, uint8_t *ss, const uint8_t *pk)
{
    uint8_t m[END_MSGBYTES];
    END_NS(randombytes)(m, END_MSGBYTES);
    return END_NS(crypto_kem_enc_derand)(ct, ss, pk, m);
}

static uint8_t ct_memneq(const uint8_t *a, const uint8_t *b, size_t len)   /* 1 if different */
{
    uint8_t d = 0;
    size_t i;
    for (i = 0; i < len; i++) d |= a[i] ^ b[i];
    return (uint8_t)((0u - (uint32_t)d) >> 31);
}

int END_NS(crypto_kem_dec_exp)(uint8_t *ss, const uint8_t *ct, const end_expkey *ek)
{
    uint8_t buf[END_MSGBYTES + END_CTBYTES];
    uint8_t r[END_MSGBYTES], ct2[END_CTBYTES];
    uint8_t fail, mask;
    unsigned i;
    pke_decrypt(buf, ek, ct);                      /* m' (whatever the decoder produced) */
    end_hash1(r, buf, END_MSGBYTES);
    pke_encrypt(ct2, ek->h_ntt, buf, r);           /* re-encryption */
    fail = ct_memneq(ct, ct2, END_CTBYTES);
    mask = (uint8_t)(0u - fail);
    for (i = 0; i < END_MSGBYTES; i++)             /* implicit rejection: m' <- k */
        buf[i] ^= mask & (buf[i] ^ ek->k[i]);
    memcpy(buf + END_MSGBYTES, ct, END_CTBYTES);
    sha3_256(ss, buf, sizeof buf);
    return 0;
}

int END_NS(crypto_kem_dec)(uint8_t *ss, const uint8_t *ct, const uint8_t *sk)
{
#ifdef END_STATIC_EXPKEY
    static end_expkey ek;
#else
    end_expkey ek;
#endif
    END_NS(expkey_from_sk)(&ek, sk);
    return END_NS(crypto_kem_dec_exp)(ss, ct, &ek);
}
