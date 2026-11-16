/*
 * Host speed test in the style of pq-crystals' test_speed (and the END paper's main.c): every KEM
 * operation is timed with the platform cycle counter over NTESTS runs, and both the median and the
 * mean are printed.  KeyGen uses fresh randomness on every run, so its mean includes the
 * seed-dependent retries of the NTRU solver.
 *   speed_host [NTESTS] [KEYGEN_TESTS]       defaults 10000 and 200
 * Used by scripts/host_baseline.sh to set this implementation next to the paper's reference code
 * and ML-KEM on the same machine.
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

static void report(const char *label, uint64_t *t, int n)
{
    uint64_t sum = 0;
    int i;
    for (i = 0; i < n; i++) sum += t[i];
    qsort(t, (size_t)n, sizeof t[0], cmp_u64);
    printf("%-22s median: %10llu  average: %10llu  (%d runs)\n", label,
           (unsigned long long)t[n / 2], (unsigned long long)(sum / (uint64_t)n), n);
}

int main(int argc, char **argv)
{
    int n = argc > 1 ? atoi(argv[1]) : 10000, nk = argc > 2 ? atoi(argv[2]) : 200, i;
    static uint8_t pk[END_PKBYTES], sk[END_SKBYTES], ct[END_CTBYTES], ss[END_SSBYTES], ss2[END_SSBYTES];
    static end_expkey ek;
    uint64_t *t = malloc(sizeof *t * (size_t)(n > nk ? n : nk)), c0;
    if (!t || n < 1 || nk < 1) return 1;
    printf("%s (%s backend), %s\n", END_NAME, END_BACKEND_NAME, cycles_unit());
    for (i = 0; i < nk; i++) {
        c0 = cycles_now();
        if (END_NS(crypto_kem_keypair)(pk, sk)) return 1;
        t[i] = cycles_now() - c0;
    }
    report("crypto_kem_keypair", t, nk);
    for (i = 0; i < n; i++) {
        c0 = cycles_now();
        END_NS(crypto_kem_enc)(ct, ss, pk);
        t[i] = cycles_now() - c0;
    }
    report("crypto_kem_enc", t, n);
    {   /* public-key decoding, part of every encapsulation (the paper's Encaps takes a decoded key) */
        static int16_t h[END_N];
        for (i = 0; i < n; i++) {
            c0 = cycles_now();
            pk_decode(h, pk);
            t[i] = cycles_now() - c0;
        }
        report("pk_decode", t, n);
    }
    for (i = 0; i < n; i++) {
        c0 = cycles_now();
        END_NS(crypto_kem_dec)(ss2, ct, sk);
        t[i] = cycles_now() - c0;
    }
    report("crypto_kem_dec", t, n);
    if (memcmp(ss, ss2, sizeof ss)) { printf("decapsulation mismatch\n"); return 1; }
    for (i = 0; i < n; i++) {
        c0 = cycles_now();
        END_NS(expkey_from_sk)(&ek, sk);
        t[i] = cycles_now() - c0;
    }
    report("expkey_from_sk", t, n);
    for (i = 0; i < n; i++) {
        c0 = cycles_now();
        END_NS(crypto_kem_dec_exp)(ss2, ct, &ek);
        t[i] = cycles_now() - c0;
    }
    report("crypto_kem_dec_exp", t, n);
    if (memcmp(ss, ss2, sizeof ss)) { printf("decapsulation mismatch\n"); return 1; }
    free(t);
    return 0;
}
