/* END KEM -- AArch64 NEON kernels for the seven-layer NTT and leaf products modulo 257. */
#include "poly.h"
#include "neon_common.h"

/* forward and inverse transform: three merged layers per pass (neon_ntt7.h).
 * forward: |in| <= 16384 and seven layers add less than 7 * 257, no input reduction; canonical output.
 * inverse: sums of layer 7 are folded to [-127, 383]; two groups of three layers multiply this by 64 at most. */
#define NTT7_FWD        poly_ntt257
#define NTT7_INV        poly_intt257
#define NTT7_P          257
#define NTT7_PINV       PINV16_257
#define NTT7_ZETAS      end_zetas257_mont
#define NTT7_ZETAS_INV  end_zetas257_inv_mont
#define NTT7_NINV       END_INV128_MONT_MOD257
#define NTT7_FQMUL(a, b) fqmul257(a, b)
#define NTT7_RED_IN(x)  (x)
#define NTT7_FWD_OUT(x) neon_red257(x)
#define NTT7_RED_SUM(x) neon_fold257(x)
#define NTT7_RED_MID    0
#define NTT7_INV_OUT(x) neon_red257(x)
#include "neon_ntt7.h"
#include "neon_basemul.h"

void poly_basemul257(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N])
{
    neon_basemul(r, a, b, end_zetas257_mont, 257, PINV16_257);
}
