/*
 * Wrap-injection test for the C decoder: replays the instances exported by golden/test_wraps.py --dump
 * (kat/wraps_end{512,1024}.txt) with the key of KAT entry 0 and checks that decode_trap recovers s
 * exactly where the golden model does (and fails where it must fail).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "end_internal.h"

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
static char line[16384];
static int read_field(FILE *fp, const char *name, char **val)
{
    size_t nl = strlen(name);
    while (fgets(line, sizeof line, fp)) {
        if (strncmp(line, name, nl) == 0 && line[nl] == ' ' && line[nl + 1] == '=') {
            *val = line + nl + 3;
            return 0;
        }
    }
    return -1;
}

int main(int argc, char **argv)
{
    const char *katpath = argc > 1 ? argv[1] :
        (END_N == 512 ? "kat/kat_end512_ref.txt" : "kat/kat_end1024_ref.txt");
    const char *wrpath = argc > 2 ? argv[2] : (END_N == 512 ? "kat/wraps_end512.txt" : "kat/wraps_end1024.txt");
    FILE *fp;
    char *v;
    static uint8_t pkb[END_PKBYTES], kf[END_N], kg[END_N], kF[END_N], kG[END_N], kw[2 * END_N], kk[END_MSGBYTES];
    static int16_t h_ntt[END_N] END_ALIGN, w[END_N] END_ALIGN, cprime[END_N] END_ALIGN, sbar[END_N] END_ALIGN;
    static uint8_t cpb[2 * END_N], sb[END_N];
    static end_expkey ek;
    char name[256];
    unsigned i;
    int ntests = 0, fails = 0;

    fp = fopen(katpath, "r");
    if (!fp) { printf("cannot open %s\n", katpath); return 2; }
    if (read_field(fp, "pk", &v) || unhex(pkb, END_PKBYTES, v) != END_PKBYTES) return 2;
    if (read_field(fp, "f", &v) || unhex(kf, END_N, v) != END_N) return 2;
    if (read_field(fp, "g", &v) || unhex(kg, END_N, v) != END_N) return 2;
    if (read_field(fp, "F", &v) || unhex(kF, END_N, v) != END_N) return 2;
    if (read_field(fp, "G", &v) || unhex(kG, END_N, v) != END_N) return 2;
    if (read_field(fp, "w", &v) || unhex(kw, 2 * END_N, v) != 2 * END_N) return 2;
    if (read_field(fp, "k", &v) || unhex(kk, END_MSGBYTES, v) != END_MSGBYTES) return 2;
    fclose(fp);
    pk_decode(h_ntt, pkb);
    for (i = 0; i < END_N; i++) w[i] = (int16_t)((uint16_t)kw[2 * i] | ((uint16_t)kw[2 * i + 1] << 8));
    expkey_from_parts(&ek, kk, h_ntt, (const int8_t *)kf, (const int8_t *)kg, (const int8_t *)kF, (const int8_t *)kG, w);

    fp = fopen(wrpath, "r");
    if (!fp) { printf("cannot open %s\n", wrpath); return 2; }
    printf("test_wraps %s: %s\n", END_NAME, wrpath);
    while (read_field(fp, "name", &v) == 0) {
        int expect_ok, ok, same;
        uint32_t diag[8] = {0};
        strncpy(name, v, sizeof name - 1); name[sizeof name - 1] = 0;
        name[strcspn(name, "\r\n")] = 0;
        if (read_field(fp, "cprime", &v) || unhex(cpb, 2 * END_N, v) != 2 * END_N) { printf("bad cprime\n"); return 2; }
        if (read_field(fp, "s", &v) || unhex(sb, END_N, v) != END_N) { printf("bad s\n"); return 2; }
        if (read_field(fp, "expect_ok", &v)) { printf("bad expect_ok\n"); return 2; }
        expect_ok = atoi(v);
        for (i = 0; i < END_N; i++) cprime[i] = (int16_t)((uint16_t)cpb[2 * i] | ((uint16_t)cpb[2 * i + 1] << 8));
        ok = decode_trap(sbar, &ek, cprime, diag);
        same = 1;
        for (i = 0; i < END_N; i++) if (((sbar[i] + 1) >> 1) != sb[i]) same = 0;
        ntests++;
        if (expect_ok) {
            if (!(ok && same)) { fails++; printf("  [FAIL] %s (found=%d same=%d gf_cand=%u GF_cand=%u S_gf=(%u,%u) S_GF=(%u,%u))\n",
                                                name, ok, same, (unsigned)diag[0], (unsigned)diag[1], (unsigned)diag[2], (unsigned)diag[3], (unsigned)diag[4], (unsigned)diag[5]); }
            else printf("  [PASS] %s (gf_cand=%u GF_cand=%u)\n", name, (unsigned)diag[0], (unsigned)diag[1]);
        } else {
            if (ok && same) { fails++; printf("  [FAIL] %s: decoder succeeded where it must not\n", name); }
            else printf("  [PASS] %s (found=%d same=%d)\n", name, ok, same);
        }
    }
    fclose(fp);
    printf("test_wraps %s: %d instances, %d failures -> %s\n", END_NAME, ntests, fails, fails ? "FAIL" : "ALL OK");
    return fails != 0 || ntests == 0;
}
