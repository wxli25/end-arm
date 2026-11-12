/*
 * END KEM -- key generation (L4 + Alg. 7 / Alg. 12):
 *   seed --SHAKE256--> k, then f, g (ternary, invertible mod 257), h = g f^-1 (NTT-257 domain),
 *   (F, G) from the NTRU solver (ntrugen profile SOLVE_END_257_*), w from compute_w (fixed-point FFT).
 *
 * compute_w() and the solver profiles are taken verbatim from the reference keygen.c (itself derived
 * from Thomas Pornin's ntrugen). END-512 consumes a continuous SHAKE256 stream.
 */
#include <string.h>
#include "end_internal.h"
#include "end_arith.h"
#include "module_p.h"
#include "big_integer.h"
#include "ntru.h"
#include "ntrugen_poly.h"
#ifdef END_KEYGEN_TRACE
#include <stdio.h>
#define KTRACE(...) do { printf(__VA_ARGS__); putchar(10); } while (0)
#else
#define KTRACE(...) do { } while (0)
#endif

#if END_N == 512
static const ntru_profile SOLVE_END_PROFILE = {
    257,
    9, 9,
    { 1, 1, 1, 2, 2, 4,  8,  15, 28, 56, 0 },
    { 1, 1, 2, 3, 6, 11, 21, 42, 82, 0 },
    { 1, 1, 2, 2, 3, 3, 3, 4, 4, 0 },
    4,
    { 0, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31 },
    { 0, 0, 1, 2, 2, 2, 2, 2, 2, 3, 3 }
};
#else
static const ntru_profile SOLVE_END_PROFILE = {
    257,
    10, 10,
    { 1, 1, 1, 2, 3, 4,  8,  16, 30, 59, 117 },
    { 1, 1, 2, 3, 6, 12, 23, 44, 87, 172 },
    { 1, 1, 1, 2, 3, 3,  3,  3,  4,  6 },
    4,
    { 0, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31 },
    { 0, 0, 1, 2, 2, 2, 2, 2, 2, 3, 3 }
};
#endif

/*
 * Compute the vector w (reference keygen.c, compute_w):
 *   w = round(q'*(gamma2*F*adj(f) + G*adj(g))/(gamma2*f*adj(f) + g*adj(g))),  q' = 12289
 * Returned value: 1 on success, 0 on error (overflow -> reject the key).
 * Value of w*R (R = 2^32 mod p = 20478) is written at the start of tmp[] (plain 32-bit).
 * RAM USAGE: 4.5*n words
 */
static int compute_w(unsigned logn, const int8_t *restrict f, const int8_t *restrict g,
                     const int8_t *restrict F, const int8_t *restrict G, uint32_t *restrict tmp)
{
    uint32_t qp = 12289;
    int32_t gamma2 = END_GAMMA2;
    size_t n = (size_t)1 << logn;
    size_t hn = n >> 1;
    uint32_t *t1 = tmp;
    uint32_t *t2 = t1 + n;
    uint32_t *t3 = t2 + n;
    uint32_t *t4 = t3 + n;
    uint32_t *t5 = t4 + n;   /* half-size */
    uint32_t p = PRIMES[0].p;
    uint32_t p0i = PRIMES[0].p0i;
    uint32_t R2 = PRIMES[0].R2;
    uint32_t *gm = t1;
    uint32_t gmv;
    uint32_t *igm;
    fxr *rt1, *rt2;
    int32_t lim1;
    fxr lim2;
    uint64_t Fdnorm = 0, Gdnorm = 0;
    size_t u;

    mp_mkgm(logn, gm, PRIMES[0].g, p, p0i);
    gmv = mp_montymul(R2, mp_set(gamma2, p), p, p0i);
    poly_mp_set_small(logn, t2, f, p);
    poly_mp_set_small(logn, t3, g, p);
    poly_mp_set_small(logn, t4, F, p);
    mp_NTT(logn, t2, gm, p, p0i);
    mp_NTT(logn, t3, gm, p, p0i);
    mp_NTT(logn, t4, gm, p, p0i);
    for (u = 0; u < hn; u ++) {
        uint32_t xf = t2[u];
        uint32_t xfa = t2[n - 1 - u];
        uint32_t xg = t3[u];
        uint32_t xga = t3[n - 1 - u];
        uint32_t xF = t4[u];
        uint32_t xFa = t4[n - 1 - u];
        uint32_t gmvf = mp_montymul(gmv, xf, p, p0i);
        uint32_t gmvfa = mp_montymul(gmv, xfa, p, p0i);
        t5[u] = mp_montymul(R2, mp_add(
            mp_montymul(xf, gmvfa, p, p0i),
            mp_montymul(xg, xga, p, p0i), p), p, p0i);
        t4[u] = mp_montymul(R2, mp_montymul(xF, gmvfa, p, p0i), p, p0i);
        t4[n - 1 - u] = mp_montymul(R2, mp_montymul(xFa, gmvf, p, p0i), p, p0i);
    }
    poly_mp_set_small(logn, t2, G, p);
    mp_NTT(logn, t2, gm, p, p0i);
    for (u = 0; u < n; u ++) {
        uint32_t xG = t2[u];
        uint32_t xga = t3[n - 1 - u];
        uint32_t xz = t4[u];
        uint32_t mga = mp_montymul(xga, R2, p, p0i);
        t1[u] = mp_add(mp_montymul(mga, xG, p, p0i), xz, p);
    }
    for (u = 0; u < hn; u ++) {
        t2[u] = t2[n - 1 - u] = t5[u];
    }
    igm = t3;
    mp_mkigm(logn, igm, PRIMES[0].ig, p, p0i);
    mp_iNTT(logn, t1, igm, p, p0i);
    mp_iNTT(logn, t2, igm, p, p0i);
    for (u = 0; u < n; u ++) {
        t1[u] = (uint32_t)mp_norm(t1[u], p);
        t2[u] = (uint32_t)mp_norm(t2[u], p);
    }
    rt1 = (fxr *)(t2 + n);
    for (u = 0; u < n; u ++) {
        rt1[u] = fxr_of(*(int32_t *)&t2[u]);
    }
    vect_FFT(logn, rt1);
    rt2 = (fxr *)t2;
    memmove(rt2, rt1, hn * sizeof *rt1);
    lim1 = (int32_t)(((uint64_t)1 << (41 - logn)) / qp);
    for (u = 0; u < n; u ++) {
        int32_t x = *(int32_t *)&t1[u];
        if (x <= -lim1 || x >= +lim1) {
            return 0;
        }
        rt1[u] = fxr_of_scaled32(((uint64_t)x * qp) << 20);
    }
    vect_FFT(logn, rt1);
    for (u = 0; u < hn; u ++) {
        fxr z1r = rt1[u];
        fxr z1i = rt1[u + hn];
        fxr z2 = rt2[u];
        if (!fxr_lt(fxr_div2e(fxr_abs(z1r), 30 - logn), z2)
            || !fxr_lt(fxr_div2e(fxr_abs(z1i), 30 - logn), z2))
        {
            return 0;
        }
        rt1[u] = fxr_div(z1r, z2);
        rt1[u + hn] = fxr_div(z1i, z2);
    }
    vect_iFFT(logn, rt1);
    lim2 = fxr_of(1 << 3);
    for (u = 0; u < n; u ++) {
        if (fxr_lt(lim2, fxr_abs(rt1[u]))) {
            return 0;
        }
        t1[u] = (uint32_t)fxr_round(fxr_mul2e(rt1[u], 12));
    }
    /* norm check of (gamma*Fd, Gd), Fd = q'*F - f*w, Gd = q'*G - g*w (values computed, not used for rejection) */
    gm = igm = t4;
    mp_mkgm(logn, gm, PRIMES[0].g, p, p0i);
    for (u = 0; u < n; u ++) {
        t1[u] = mp_montymul(R2, mp_set(*(int32_t *)&t1[u], p), p, p0i);
    }
    mp_NTT(logn, t1, gm, p, p0i);
    for (u = 0; u < n; u ++) {
        t2[u] = mp_set(f[u], p);
        t3[u] = mp_set((int32_t)qp * F[u], p);
    }
    mp_NTT(logn, t2, gm, p, p0i);
    mp_NTT(logn, t3, gm, p, p0i);
    for (u = 0; u < n; u ++) {
        t2[u] = mp_sub(t3[u], mp_montymul(t2[u], t1[u], p, p0i), p);
    }
    mp_mkigm(logn, igm, PRIMES[0].ig, p, p0i);
    mp_iNTT(logn, t2, igm, p, p0i);
    for (u = 0; u < n; u ++) {
        int32_t x = mp_norm(t2[u], p);
        Fdnorm += (uint64_t)((int64_t)x * (int64_t)x);
    }
    for (u = 0; u < n; u ++) {
        t2[u] = mp_set(g[u], p);
        t3[u] = mp_set((int32_t)qp * G[u], p);
    }
    mp_mkgm(logn, gm, PRIMES[0].g, p, p0i);
    mp_NTT(logn, t2, gm, p, p0i);
    mp_NTT(logn, t3, gm, p, p0i);
    for (u = 0; u < n; u ++) {
        t2[u] = mp_sub(t3[u], mp_montymul(t2[u], t1[u], p, p0i), p);
    }
    mp_mkigm(logn, igm, PRIMES[0].ig, p, p0i);
    mp_iNTT(logn, t2, igm, p, p0i);
    for (u = 0; u < n; u ++) {
        int32_t x = mp_norm(t2[u], p);
        Gdnorm += (uint64_t)((int64_t)x * (int64_t)x);
    }
    (void)Fdnorm; (void)Gdnorm;
    mp_iNTT(logn, t1, igm, p, p0i);
    for (u = 0; u < n; u ++) {
        t1[u] = (uint32_t)mp_norm(t1[u], p);
    }
    return 1;
}

/* ------------------------------------------------------------------------- */
/* ternary sampling (reference poly.c / keygen.c)                             */
/* ------------------------------------------------------------------------- */
#if END_N == 512
/* 384 consecutive SHAKE256 bytes -> 512 samples; k < 9 gives -1, k > 54 gives +1. */
void ternary_sample(int8_t a[END_N], end_xof_state *state)
{
    uint8_t buf[384];
    unsigned i = 0, j;
    end_xof_squeeze(buf, sizeof buf, state);
    for (j = 0; j < END_N; j += 4) {
        uint8_t x = buf[i++], y = buf[i++], z = buf[i++], k;
        int8_t m_neg, m_pos;
        k = x & 63; m_neg = (int8_t)-(k < 9); m_pos = (int8_t)-(k > 54); a[j]     = (int8_t)((m_neg & -1) | (m_pos & 1));
        k = y & 63; m_neg = (int8_t)-(k < 9); m_pos = (int8_t)-(k > 54); a[j + 1] = (int8_t)((m_neg & -1) | (m_pos & 1));
        k = z & 63; m_neg = (int8_t)-(k < 9); m_pos = (int8_t)-(k > 54); a[j + 2] = (int8_t)((m_neg & -1) | (m_pos & 1));
        k = (uint8_t)((((x >> 6) & 3) << 4) | (((y >> 6) & 3) << 2) | ((z >> 6) & 3));
        m_neg = (int8_t)-(k < 9); m_pos = (int8_t)-(k > 54); a[j + 3] = (int8_t)((m_neg & -1) | (m_pos & 1));
    }
}
#else
/* 680 bytes (5 SHAKE256 blocks); f = (b1 & b2 - b3 & b4) * b5 with independent bits (P(+-1) = 3/32) */
void ternary_sample(int8_t f[END_N], end_xof_state *state)
{
    uint8_t buf[680];
    int8_t tmp1[2 * END_N], tmp2[2 * END_N];
    unsigned i, j, t = 0;
    end_xof_squeezeblocks(buf, 5, state);
    for (i = 0; i < 8; i++)
        for (j = 0; j < 256; j++)
            tmp1[i * 256 + j] = (int8_t)((buf[t + j] >> i) & 1);
    t += 256;
    for (i = 0; i < END_N; i++) tmp1[i] = (int8_t)(tmp1[i] & tmp1[i + END_N]);
    for (i = 0; i < 8; i++)
        for (j = 0; j < 256; j++)
            tmp2[i * 256 + j] = (int8_t)((buf[t + j] >> i) & 1);
    t += 256;
    for (i = 0; i < END_N; i++) { tmp2[i] = (int8_t)(tmp2[i] & tmp2[i + END_N]); tmp1[i] = (int8_t)(tmp1[i] - tmp2[i]); }
    for (i = 0; i < 8; i++)
        for (j = 0; j < 128; j++)
            tmp2[i * 128 + j] = (int8_t)((buf[t + j] >> i) & 1);
    for (i = 0; i < END_N; i++) f[i] = (int8_t)(tmp1[i] * tmp2[i]);
}
#endif

/* ------------------------------------------------------------------------- */
/* key pair                                                                   */
/* ------------------------------------------------------------------------- */
int END_NS(crypto_kem_keypair_derand)(uint8_t *pk, uint8_t *sk, const uint8_t seed[END_SEEDBYTES],
                                      void *tmp, size_t tmplen)
{
    end_xof_state state;
    int16_t f_inv[END_N] END_ALIGN, t1[END_N] END_ALIGN, t2[END_N] END_ALIGN;
    int8_t *f = (int8_t *)(sk + END_SK_OFF_F), *g = (int8_t *)(sk + END_SK_OFF_G);
    uintptr_t utmp1 = (uintptr_t)tmp, utmp2 = (utmp1 + 7) & ~(uintptr_t)7;
    uint32_t *tt32 = (uint32_t *)utmp2;
    int8_t *tG = (int8_t *)tt32;
    int8_t *tF = tG + END_N;
    uint32_t *tw = (uint32_t *)(tt32 + END_N);
    unsigned i;

    if (tmplen < (size_t)(utmp2 - utmp1) + 22 * END_N)
        return -1;
    end_xof_absorb_once(&state, seed, END_SEEDBYTES);
    end_xof_squeeze(sk + END_SK_OFF_K, END_MSGBYTES, &state);        /* implicit-rejection secret k */
    for (;;) {
        do {
            ternary_sample(f, &state);
            for (i = 0; i < END_N; i++) t1[i] = f[i];
            poly_ntt257(t1);
        } while (poly_baseinv257(f_inv, t1));
        do {
            ternary_sample(g, &state);
            for (i = 0; i < END_N; i++) t1[i] = g[i];
            poly_ntt257(t1);
        } while (poly_baseinv257(t2, t1));
        poly_basemul257(t2, f_inv, t1);                                /* h = g f^-1 (NTT domain) */
        pk_encode(pk, t2);
        KTRACE("keygen: f,g sampled (sum f=%d, sum g=%d)", (int)0, (int)0);
        /* solver convention: solve_NTRU(prof, logn, "f" := g, "g" := f) yields g*tF - f*tG = q */
        {
            int err = solve_NTRU(&SOLVE_END_PROFILE, END_LOGN, g, f, tt32);
            KTRACE("keygen: solve_NTRU -> %d", err);
            if (err != SOLVE_OK)
                continue;
        }
        if (!compute_w(END_LOGN, f, g, tF, tG, tw)) {
            KTRACE("keygen: compute_w failed");
            continue;
        }
        KTRACE("keygen: compute_w ok");
        memcpy(sk + END_SK_OFF_FF, tF, END_N);
        memcpy(sk + END_SK_OFF_GG, tG, END_N);
        for (i = 0; i < END_N; i++) {
            int32_t wi = ((int32_t)tw[i]) / 20478;                     /* remove R = 2^32 mod p (exact) */
            sk[END_SK_OFF_W + 2 * i] = (uint8_t)wi;
            sk[END_SK_OFF_W + 2 * i + 1] = (uint8_t)(wi >> 8);
        }
        break;
    }
    end_xof_release(&state);
    memcpy(sk + END_SK_OFF_PK, pk, END_PKBYTES);
    return 0;
}

int END_NS(crypto_kem_keypair)(uint8_t *pk, uint8_t *sk)
{
    uint8_t seed[END_SEEDBYTES];
#ifdef END_STATIC_KEYGEN_TMP
    static uint8_t tmp[END_KEYGEN_TMPBYTES];
#else
    uint8_t tmp[END_KEYGEN_TMPBYTES];
#endif
    if (END_NS(randombytes)(seed, sizeof seed) != 0)
        return -1;
    return END_NS(crypto_kem_keypair_derand)(pk, sk, seed, tmp, sizeof tmp);
}
