/*
 * END KEM -- randombytes() for the host builds (Linux getrandom / Windows CryptGenRandom).
 * Bare-metal targets define END_RANDOMBYTES_EXTERNAL and provide the symbol themselves
 * (platform/cortex-m/..., e.g. from a TRNG or a DRBG seeded at boot).
 */
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include "end_api.h"

#ifndef END_RANDOMBYTES_EXTERNAL
#if defined(_WIN32)
#include <windows.h>
#include <wincrypt.h>
void END_NS(randombytes)(uint8_t *out, size_t outlen)
{
    HCRYPTPROV ctx;
    if (!CryptAcquireContext(&ctx, NULL, NULL, PROV_RSA_FULL, CRYPT_VERIFYCONTEXT)) abort();
    while (outlen > 0) {
        size_t len = outlen > 1048576 ? 1048576 : outlen;
        if (!CryptGenRandom(ctx, (DWORD)len, (BYTE *)out)) abort();
        out += len; outlen -= len;
    }
    CryptReleaseContext(ctx, 0);
}
#else
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
void END_NS(randombytes)(uint8_t *out, size_t outlen)
{
    static int fd = -1;
    if (fd == -1) {
        for (;;) {
            fd = open("/dev/urandom", O_RDONLY);
            if (fd != -1) break;
            if (errno != EINTR) abort();
        }
    }
    while (outlen > 0) {
        ssize_t r = read(fd, out, outlen);
        if (r < 0) { if (errno == EINTR) continue; abort(); }
        out += r; outlen -= (size_t)r;
    }
}
#endif
#endif /* !END_RANDOMBYTES_EXTERNAL */
