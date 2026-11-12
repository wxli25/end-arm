/*
 * END KEM -- parameter sets and sizes.
 *
 * Select the parameter set at compile time with -DEND_N=512 (NIST-I, default) or -DEND_N=1024 (NIST-V).
 * Sources: paper Table 2 (sizes, bounds) and the reference implementation commit dde0a77
 * (compression constants, gamma weights, hash choices, message layout).
 */
#ifndef END_PARAMS_H
#define END_PARAMS_H

#include <stdint.h>
#include <stddef.h>

#ifndef END_N
#define END_N 512
#endif

#define END_Q      257                 /* decryption modulus q              */
#define END_QP     12289               /* precision modulus Q               */
#define END_QQP    3158273             /* q*Q                               */
#define END_2QQP   6316546             /* 2*q*Q (modulus of the doubled GF residue) */
#define END_P3     3329                /* auxiliary NTT prime (s~' and r*f, exact for |x| < 1664) */

#if END_N == 512
#define END_LOGN        9
#define END_LEAF        4              /* degree of the NTT-257 / NTT-3329 leaves (x^4 - zeta) */
#define END_LOGLEAF     2
#define END_DE          64             /* ciphertext compression base d_e   */
#define END_DELTA       1              /* centring offset of the compression error (mu~ = 2 mu' - delta) */
#define END_MSGBYTES    32             /* |M| = n/2 bits                    */
#define END_PKBYTES     514
#define END_CTBYTES     384
#define END_GAMMA2      6              /* weight used in w = round(Q (G g~ + GAMMA2 F f~)/(g g~ + GAMMA2 f f~)) (reference: 6, paper: 4) */
#define END_BGF0        2              /* gf-stage candidate weight bound B_gf,0  (Table 2) */
#define END_BGF2_0      2              /* GF-stage candidate weight bound B_GF,0  (reference uses 2) */
#define END_HASH1_BYTES 32             /* H1 = H2 = SHA3-256 */
#define END_NAME        "END-512"
#elif END_N == 1024
#define END_LOGN        10
#define END_LEAF        8
#define END_LOGLEAF     3
#define END_DE          90
#define END_DELTA       0
#define END_MSGBYTES    64
#define END_PKBYTES     1027
#define END_CTBYTES     832
#define END_GAMMA2      3
#define END_BGF0        2
#define END_BGF2_0      2
#define END_HASH1_BYTES 64             /* H1 = H2 = SHA3-512 */
#define END_NAME        "END-1024"
#else
#error "END_N must be 512 or 1024"
#endif

#define END_B1          1              /* coefficient-magnitude bound of candidate wraps (B_1) */
#ifndef END_M0_GF
#define END_M0_GF       2              /* FCL output size m0 of the gf stage */
#endif
#ifndef END_M0_GF2
#define END_M0_GF2      2              /* FCL output size m0 of the GF stage */
#endif
#define END_FCL_MAX 8
#if END_M0_GF < 1 || END_M0_GF > END_FCL_MAX || END_M0_GF2 < 1 || END_M0_GF2 > END_FCL_MAX
#error "END_M0_GF and END_M0_GF2 must be in [1, 8]"
#endif
#define END_NCAND_GF    (1 + 2 * END_M0_GF + 4 * ((END_M0_GF * (END_M0_GF - 1)) / 2))   /* |E_cand| for B0 = 2, B1 = 1 */
#define END_NCAND_GF2   (1 + 2 * END_M0_GF2 + 4 * ((END_M0_GF2 * (END_M0_GF2 - 1)) / 2))

/* alignment of polynomial buffers (NEON 128-bit loads, 32-bit packed words on Cortex-M) */
#if defined(__GNUC__)
#define END_ALIGN __attribute__((aligned(16)))
#else
#define END_ALIGN
#endif

#define END_SEEDBYTES   32
#define END_SSBYTES     32             /* shared secret = SHA3-256(m || ct)  */
#define END_FG_LIMIT    31             /* |F|_inf, |G|_inf <= 31 guaranteed by the solver profile */

/* compact secret key: k (MSGBYTES) || pk || f || g || F || G (int8, n each) || w (int16 LE, n) */
#define END_SKBYTES     (END_MSGBYTES + END_PKBYTES + 4 * END_N + 2 * END_N)
#define END_SK_OFF_K    0
#define END_SK_OFF_PK   (END_MSGBYTES)
#define END_SK_OFF_F    (END_MSGBYTES + END_PKBYTES)
#define END_SK_OFF_G    (END_SK_OFF_F + END_N)
#define END_SK_OFF_FF   (END_SK_OFF_G + END_N)   /* capital F */
#define END_SK_OFF_GG   (END_SK_OFF_FF + END_N)  /* capital G */
#define END_SK_OFF_W    (END_SK_OFF_GG + END_N)

/* key generation scratch (bytes), as in the reference: 22*n for the solver + compute_w */
#ifndef END_KEYGEN_TMPBYTES
#define END_KEYGEN_TMPBYTES (22 * END_N + 16)
#endif

/* symbol namespacing so that both parameter sets can live in one binary */
#define END_CAT2(a, b) a##b
#define END_CAT(a, b) END_CAT2(a, b)
#if END_N == 512
#define END_NS(s) END_CAT(end512_, s)
#else
#define END_NS(s) END_CAT(end1024_, s)
#endif

#endif /* END_PARAMS_H */
