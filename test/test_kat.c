/*
 * KAT replay: test_kat [file] [count] [--import | --generate]
 *   --import uses serialized keys from an older vector set.
 *   --generate takes seeds and messages from the input and prints fresh vectors.
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
            return line[nl + 2] == ' ' && strcspn(p, "\r\n") == 2 * len &&
                   unhex(out, len, p) == len ? 0 : -1;
        }
    }
    return ferror(fp) ? -1 : 1;
}

static void print_field(const char *name, const uint8_t *buf, size_t len)
{
    size_t i;
    printf("%s = ", name);
    for (i = 0; i < len; i++)
        printf("%02x", (unsigned)buf[i]);
    putchar('\n');
}

int main(int argc, char **argv)
{
    const char *path = argc > 1 ? argv[1] :
        (END_N == 512 ? "kat/kat_end512_shake256.txt" : "kat/kat_end1024_ref.txt");
    int maxn = argc > 2 ? atoi(argv[2]) : (END_N == 512 ? 100 : 60);
    int import = argc > 3 && strcmp(argv[3], "--import") == 0;
    int generate = argc > 3 && strcmp(argv[3], "--generate") == 0;
    FILE *fp = fopen(path, "r");
    static uint8_t pk[END_PKBYTES], sk[END_SKBYTES], ct[END_CTBYTES], ss[END_SSBYTES], ss2[END_SSBYTES];
    static uint8_t seed[32], m[END_MSGBYTES], kpk[END_PKBYTES], kf[END_N], kg[END_N], kF[END_N], kG[END_N], kw[2 * END_N],
                   kk[END_MSGBYTES], kct[END_CTBYTES], kss[END_SSBYTES];
    static uint8_t tmp[END_KEYGEN_TMPBYTES];
    static end_expkey ek;
    int count = 0, fails = 0;
    if (!fp) { printf("cannot open %s\n", path); return 2; }
    if (maxn < 1 || (argc > 3 && !import && !generate)) {
        fclose(fp);
        return 2;
    }
    if (!generate)
        printf("test_kat %s: %s%s\n", END_NAME, path, import ? " (imported keys)" : "");
    while (count < maxn) {
        int bad = 0;
        if (read_field(fp, "keygen_seed", seed, sizeof seed) != 0)
            break;
        if (generate) {
            if (read_field(fp, "m", m, sizeof m) ||
                END_NS(crypto_kem_keypair_derand)(pk, sk, seed, tmp, sizeof tmp) ||
                END_NS(crypto_kem_enc_derand)(ct, ss, pk, m)) {
                fclose(fp);
                return 2;
            }
            printf("count = %d\n", count++);
            print_field("keygen_seed", seed, sizeof seed);
            print_field("pk", pk, sizeof pk);
            print_field("f", sk + END_SK_OFF_F, END_N);
            print_field("g", sk + END_SK_OFF_G, END_N);
            print_field("F", sk + END_SK_OFF_FF, END_N);
            print_field("G", sk + END_SK_OFF_GG, END_N);
            print_field("w", sk + END_SK_OFF_W, 2 * END_N);
            print_field("k", sk + END_SK_OFF_K, END_MSGBYTES);
            print_field("m", m, sizeof m);
            print_field("ct", ct, sizeof ct);
            print_field("ss", ss, sizeof ss);
            putchar('\n');
            continue;
        }
        if (read_field(fp, "pk", kpk, END_PKBYTES) || read_field(fp, "f", kf, END_N) || read_field(fp, "g", kg, END_N) ||
            read_field(fp, "F", kF, END_N) || read_field(fp, "G", kG, END_N) || read_field(fp, "w", kw, 2 * END_N) ||
            read_field(fp, "k", kk, END_MSGBYTES) || read_field(fp, "m", m, END_MSGBYTES) ||
            read_field(fp, "ct", kct, END_CTBYTES) || read_field(fp, "ss", kss, END_SSBYTES)) {
            printf("malformed KAT entry %d\n", count);
            fclose(fp);
            return 2;
        }
        if (import) {
            memcpy(sk + END_SK_OFF_PK, kpk, END_PKBYTES);
            memcpy(sk + END_SK_OFF_K, kk, END_MSGBYTES);
            memcpy(sk + END_SK_OFF_F, kf, END_N);
            memcpy(sk + END_SK_OFF_G, kg, END_N);
            memcpy(sk + END_SK_OFF_FF, kF, END_N);
            memcpy(sk + END_SK_OFF_GG, kG, END_N);
            memcpy(sk + END_SK_OFF_W, kw, 2 * END_N);
        } else {
            if (END_NS(crypto_kem_keypair_derand)(pk, sk, seed, tmp, sizeof tmp)) {
                fclose(fp);
                return 2;
            }
            if (memcmp(pk, kpk, END_PKBYTES)) { printf("entry %d: pk mismatch\n", count); bad = 1; }
            if (memcmp(sk + END_SK_OFF_K, kk, END_MSGBYTES)) { printf("entry %d: k mismatch\n", count); bad = 1; }
            if (memcmp(sk + END_SK_OFF_F, kf, END_N)) { printf("entry %d: f mismatch\n", count); bad = 1; }
            if (memcmp(sk + END_SK_OFF_G, kg, END_N)) { printf("entry %d: g mismatch\n", count); bad = 1; }
            if (memcmp(sk + END_SK_OFF_FF, kF, END_N)) { printf("entry %d: F mismatch\n", count); bad = 1; }
            if (memcmp(sk + END_SK_OFF_GG, kG, END_N)) { printf("entry %d: G mismatch\n", count); bad = 1; }
            if (memcmp(sk + END_SK_OFF_W, kw, 2 * END_N)) { printf("entry %d: w mismatch\n", count); bad = 1; }
        }
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
    if (count != maxn) {
        fprintf(stderr, "KAT expected %d entries, read %d\n", maxn, count);
        return 2;
    }
    if (generate)
        return fflush(stdout) != 0 || ferror(stdout) != 0;
    printf("test_kat %s: %d entries, %d failures -> %s\n", END_NAME, count, fails, fails ? "FAIL" : "ALL OK");
    return fails != 0;
}
