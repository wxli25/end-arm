/*
 * Benchmarks: KEM operations and the hot kernels of the selected backend.
 * Prints the median of N runs in cycles (platform counter) and in wall-clock microseconds (when available).
 *   bench [N] [--insn]     --insn : single-iteration mode for QEMU instruction counting (see scripts/qemu_icount.sh)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "end_internal.h"
#include "cycles.h"

static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

#define MAXN 1024
static uint64_t tc[MAXN], tw[MAXN];

#define BENCH(label, N, code) do {                                                   \
        int it_;                                                                     \
        for (it_ = 0; it_ < (N); it_++) {                                            \
            uint64_t c0 = cycles_now(), w0 = wall_ns();                              \
            code;                                                                    \
            tc[it_] = cycles_now() - c0; tw[it_] = wall_ns() - w0;                   \
        }                                                                            \
        qsort(tc, (N), sizeof tc[0], cmp_u64); qsort(tw, (N), sizeof tw[0], cmp_u64); \
        printf("  %-28s %12lu %s   %10lu us\n", label, (unsigned long)tc[(N) / 2], cycles_unit(), (unsigned long)(tw[(N) / 2] / 1000)); /* %lu only: newlib-nano has no %llu/%f */ \
    } while (0)

int main(int argc, char **argv)
{
    int N = argc > 1 ? atoi(argv[1]) : 100;
    int insn = argc > 2 && strcmp(argv[2], "--insn") == 0;
    static uint8_t pk[END_PKBYTES], sk[END_SKBYTES], ct[END_CTBYTES], ss[END_SSBYTES], ss2[END_SSBYTES];
    static uint8_t seed[END_SEEDBYTES], m[END_MSGBYTES];
    static uint8_t tmp[END_KEYGEN_TMPBYTES];
    static end_expkey ek;
    static int16_t a[END_N] END_ALIGN, b[END_N] END_ALIGN, c[END_N] END_ALIGN, sbar[END_N] END_ALIGN;
    static uint8_t buf[1024];
    int i, NK;
    if (N > MAXN) N = MAXN;
    NK = N < 10 ? N : 10;      /* keygen is slow */
    for (i = 0; i < END_SEEDBYTES; i++) seed[i] = (uint8_t)i;
    for (i = 0; i < END_MSGBYTES; i++) m[i] = (uint8_t)(3 * i + 1);
    printf("bench %s backend=%s N=%d\n", END_NAME, END_BACKEND_NAME, N);
    if (insn) {
        /* one call per operation, selected by argv[3]; used with the QEMU 'libinsn' plugin */
        const char *op = argc > 3 ? argv[3] : "all";
        END_NS(crypto_kem_keypair_derand)(pk, sk, seed, tmp, sizeof tmp);
        END_NS(expkey_from_sk)(&ek, sk);
        END_NS(crypto_kem_enc_derand)(ct, ss, pk, m);
        if (!strcmp(op, "keypair")) END_NS(crypto_kem_keypair_derand)(pk, sk, seed, tmp, sizeof tmp);
        else if (!strcmp(op, "enc")) END_NS(crypto_kem_enc_derand)(ct, ss, pk, m);
        else if (!strcmp(op, "dec")) END_NS(crypto_kem_dec)(ss2, ct, sk);
        else if (!strcmp(op, "dec_exp")) END_NS(crypto_kem_dec_exp)(ss2, ct, &ek);
        else if (!strcmp(op, "expkey")) END_NS(expkey_from_sk)(&ek, sk);
        else if (!strcmp(op, "ntt257")) poly_ntt257(a);
        else if (!strcmp(op, "intt257")) poly_intt257(a);
        else if (!strcmp(op, "basemul257")) poly_basemul257(c, a, b);
        else if (!strcmp(op, "ntt3329")) poly_ntt3329(a);
        else if (!strcmp(op, "intt3329")) poly_intt3329(a);
        else if (!strcmp(op, "basemul3329")) poly_basemul3329(c, a, b);
        else if (!strcmp(op, "ntt12289")) poly_ntt12289(a);
        else if (!strcmp(op, "intt12289")) poly_intt12289(a);
        else if (!strcmp(op, "rot16")) poly_rot16(c, a, 123);
        else if (!strcmp(op, "sha3_256")) sha3_256(buf, buf, 416);
        else if (!strcmp(op, "decode_trap")) { for (i = 0; i < END_N; i++) c[i] = (int16_t)(4 * (i % 64) + 3); decode_trap(sbar, &ek, c, NULL); }
        printf("insn mode: op=%s done\n", op);
        return 0;
    }
    BENCH("keypair", NK, END_NS(crypto_kem_keypair_derand)(pk, sk, seed, tmp, sizeof tmp));
    BENCH("enc", N, END_NS(crypto_kem_enc_derand)(ct, ss, pk, m));
    BENCH("dec (compact sk)", N, END_NS(crypto_kem_dec)(ss2, ct, sk));
    if (memcmp(ss, ss2, END_SSBYTES)) printf("  !! dec mismatch\n");
    BENCH("expkey_from_sk", N, END_NS(expkey_from_sk)(&ek, sk));
    BENCH("dec_exp (expanded sk)", N, END_NS(crypto_kem_dec_exp)(ss2, ct, &ek));
    if (memcmp(ss, ss2, END_SSBYTES)) printf("  !! dec_exp mismatch\n");
    for (i = 0; i < END_N; i++) { a[i] = (int16_t)(i * 7 % 257); b[i] = (int16_t)(i * 13 % 257); }
    BENCH("poly_ntt257", N, poly_ntt257(a));
    BENCH("poly_intt257", N, poly_intt257(a));
    BENCH("poly_basemul257", N, poly_basemul257(c, a, b));
    BENCH("poly_ntt3329", N, poly_ntt3329(a));
    BENCH("poly_intt3329", N, poly_intt3329(a));
    BENCH("poly_basemul3329", N, poly_basemul3329(c, a, b));
    BENCH("poly_ntt12289", N, poly_ntt12289(a));
    BENCH("poly_intt12289", N, poly_intt12289(a));
    BENCH("poly_pointwise12289", N, poly_pointwise12289(c, a, b));
    BENCH("poly_rot16", N, poly_rot16(c, a, (uint32_t)(i & (END_N - 1))));
    BENCH("sha3_256 (416 B)", N, sha3_256(buf, buf, 416));
    for (i = 0; i < END_N; i++) c[i] = (int16_t)(4 * (i % 64) + 3);
    BENCH("decode_trap", N, decode_trap(sbar, &ek, c, NULL));
#ifdef END_PLATFORM_STACK
    {
        extern uint32_t platform_stack_used(void);
        printf("  stack high-water mark (all operations): %u bytes", (unsigned)platform_stack_used()); putchar(10);
    }
#endif
    return 0;
}
