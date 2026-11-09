/*
 * KAT test: replay the reference implementation's known-answer file (kat/kat_end{512,1024}_ref.txt)
 *   keypair_derand(seed)  -> pk, f, g, F, G, w, k must match
 *   enc_derand(pk, m)     -> ct, ss must match
 *   dec(sk, ct)           -> ss must match (also through the expanded-key path)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "end_api.h"

static int hexval(int c)
{
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}
static size_t unhex(uint8_t *out, size_t max, const char *s)
{
    size_t n = 0;
    while (s[0] && s[1] && n < max) {
        int a = hexval(s[0]), b = hexval(s[1]);
        if (a < 0 || b < 0) break;
        out[n++] = (uint8_t)(a * 16 + b);
        s += 2;
    }
    return n;
}
static int read_field(FILE *fp, const char *name, uint8_t *out, size_t len)
{
    static char line[16384];
    size_t nl = strlen(name);
    while (fgets(line, sizeof line, fp)) {
        if (strncmp(line, name, nl) == 0 && line[nl] == ' ' && line[nl + 1] == '=') {
            const char *p = line + nl + 3;
            return unhex(out, len, p) == len ? 0 : -1;
        }
    }
    return -1;
}

int main(int argc, char **argv)
{
    const char *path = argc > 1 ? argv[1] : "kat/kat_end512_ref.txt";
    int maxn = argc > 2 ? atoi(argv[2]) : 1000000;
    FILE *fp = fopen(path, "r");
    static uint8_t pk[END_PKBYTES], sk[END_SKBYTES], ct[END_CTBYTES], ss[END_SSBYTES], ss2[END_SSBYTES];
    static uint8_t seed[32], m[END_MSGBYTES], kpk[END_PKBYTES], kf[END_N], kg[END_N], kF[END_N], kG[END_N], kw[2 * END_N],
                   kk[END_MSGBYTES], kct[END_CTBYTES], kss[END_SSBYTES];
    static uint8_t tmp[END_KEYGEN_TMPBYTES];
    static end_expkey ek;
    int count = 0, fails = 0;
    if (!fp) { printf("cannot open %s\n", path); return 2; }
    printf("test_kat %s: %s\n", END_NAME, path);
    while (count < maxn && read_field(fp, "keygen_seed", seed, 32) == 0) {
        int bad = 0;
        if (read_field(fp, "pk", kpk, END_PKBYTES) || read_field(fp, "f", kf, END_N) || read_field(fp, "g", kg, END_N) ||
            read_field(fp, "F", kF, END_N) || read_field(fp, "G", kG, END_N) || read_field(fp, "w", kw, 2 * END_N) ||
            read_field(fp, "k", kk, END_MSGBYTES) || read_field(fp, "m", m, END_MSGBYTES) ||
            read_field(fp, "ct", kct, END_CTBYTES) || read_field(fp, "ss", kss, END_SSBYTES)) {
            printf("malformed KAT entry %d\n", count); return 2;
        }
        END_NS(crypto_kem_keypair_derand)(pk, sk, seed, tmp, sizeof tmp);
        if (memcmp(pk, kpk, END_PKBYTES)) { printf("entry %d: pk mismatch\n", count); bad = 1; }
        if (memcmp(sk + END_SK_OFF_K, kk, END_MSGBYTES)) { printf("entry %d: k mismatch\n", count); bad = 1; }
        if (memcmp(sk + END_SK_OFF_F, kf, END_N)) { printf("entry %d: f mismatch\n", count); bad = 1; }
        if (memcmp(sk + END_SK_OFF_G, kg, END_N)) { printf("entry %d: g mismatch\n", count); bad = 1; }
        if (memcmp(sk + END_SK_OFF_FF, kF, END_N)) { printf("entry %d: F mismatch\n", count); bad = 1; }
        if (memcmp(sk + END_SK_OFF_GG, kG, END_N)) { printf("entry %d: G mismatch\n", count); bad = 1; }
        if (memcmp(sk + END_SK_OFF_W, kw, 2 * END_N)) { printf("entry %d: w mismatch\n", count); bad = 1; }
        END_NS(crypto_kem_enc_derand)(ct, ss, kpk, m);
        if (memcmp(ct, kct, END_CTBYTES)) { printf("entry %d: ct mismatch\n", count); bad = 1; }
        if (memcmp(ss, kss, END_SSBYTES)) { printf("entry %d: ss (enc) mismatch\n", count); bad = 1; }
        END_NS(crypto_kem_dec)(ss2, kct, sk);
        if (memcmp(ss2, kss, END_SSBYTES)) { printf("entry %d: ss (dec) mismatch\n", count); bad = 1; }
        END_NS(expkey_from_sk)(&ek, sk);
        END_NS(crypto_kem_dec_exp)(ss2, kct, &ek);
        if (memcmp(ss2, kss, END_SSBYTES)) { printf("entry %d: ss (dec_exp) mismatch\n", count); bad = 1; }
        fails += bad;
        count++;
    }
    fclose(fp);
    printf("test_kat %s: %d entries, %d failures -> %s\n", END_NAME, count, fails, fails ? "FAIL" : "ALL OK");
    return fails != 0 || count == 0;
}
