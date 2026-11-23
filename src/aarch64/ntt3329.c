/* END KEM -- AArch64 NEON kernels for the seven-layer NTT and leaf products modulo 3329. */
#include "poly.h"
#include "neon_common.h"

#define P 3329

/* forward and inverse transform: three merged layers per pass (neon_ntt7.h).
 * forward: inputs reduced to |a| <= 1665 when they are read, seven lazy layers (< 1665 + 7 * 3329), canonical output.
 * inverse: sums of layer 7 reduced to |a| <= 1665; after three more layers the four sums (< 8 * 3329) are reduced
 *          again, so the last three layers stay below 8 * 3329 as well; exact centred output. */
static inline int16x8_t barrett3329v(int16x8_t x)
{
    return neon_barrett(x, vdupq_n_s16(BARRETT_V_3329), vdupq_n_s16(P));
}
static inline int16x8_t canon3329v(int16x8_t x)
{
    return neon_canon(barrett3329v(x), vdupq_n_s16(P));
}
static inline int16x8_t center3329v(int16x8_t x)
{
    return neon_center_exact(barrett3329v(x), vdupq_n_s16(P), vdupq_n_s16(1664));
}
#define NTT7_FWD        poly_ntt3329
#define NTT7_INV        poly_intt3329
#define NTT7_P          3329
#define NTT7_PINV       PINV16_3329
#define NTT7_ZETAS      end_zetas3329_mont
#define NTT7_ZETAS_INV  end_zetas3329_inv_mont
#define NTT7_NINV       END_INV128_MONT_MOD3329
#define NTT7_FQMUL(a, b) fqmul3329(a, b)
#define NTT7_RED_IN(x)  barrett3329v(x)
#define NTT7_FWD_OUT(x) canon3329v(x)
#define NTT7_RED_SUM(x) barrett3329v(x)
#define NTT7_RED_MID    1
#define NTT7_INV_OUT(x) center3329v(x)
#include "neon_ntt7.h"
#include "neon_basemul.h"

void poly_basemul3329(int16_t r[END_N], const int16_t a[END_N], const int16_t b[END_N])
{
    neon_basemul(r, a, b, end_zetas3329_mont, 3329, PINV16_3329);
}
