/*
 * END KEM -- randombytes() for the host builds (/dev/urandom or Windows CryptGenRandom).
 * Bare-metal targets define END_RANDOMBYTES_EXTERNAL and provide the symbol themselves
 * (platform/cortex-m/..., e.g. from a TRNG or a DRBG seeded at boot).
 */
#include <stdint.h>
#include <stddef.h>
#include "end_api.h"

#ifndef END_RANDOMBYTES_EXTERNAL
#if defined(_WIN32)
#include <windows.h>
#include <wincrypt.h>
int END_NS(randombytes)(uint8_t *out, size_t outlen)
{
    HCRYPTPROV ctx;
    if (!CryptAcquireContext(&ctx, NULL, NULL, PROV_RSA_FULL, CRYPT_VERIFYCONTEXT))
        return -1;
    while (outlen > 0) {
        size_t len = outlen > 1048576 ? 1048576 : outlen;
        if (!CryptGenRandom(ctx, (DWORD)len, (BYTE *)out)) {
            CryptReleaseContext(ctx, 0);
            return -1;
        }
        out += len; outlen -= len;
    }
    CryptReleaseContext(ctx, 0);
    return 0;
}
#else
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
int END_NS(randombytes)(uint8_t *out, size_t outlen)
{
    int fd;
    do {
        fd = open("/dev/urandom", O_RDONLY | O_CLOEXEC);
    } while (fd < 0 && errno == EINTR);
    if (fd < 0)
        return -1;
    while (outlen > 0) {
        size_t len = outlen > 1048576 ? 1048576 : outlen;
        ssize_t r = read(fd, out, len);
        if (r < 0 && errno == EINTR)
            continue;
        if (r <= 0) {
            close(fd);
            return -1;
        }
        out += r; outlen -= (size_t)r;
    }
    close(fd);
    return 0;
}
#endif
#endif /* !END_RANDOMBYTES_EXTERNAL */
