/*
 * END KEM -- AArch64 NEON template: the 7-layer NTT of one modulus, forward and inverse, with merged layers.
 *
 * Included by ntt257.c and ntt3329.c after defining
 *   NTT7_FWD, NTT7_INV        names of the two functions
 *   NTT7_P, NTT7_PINV         modulus and p^-1 mod 2^16
 *   NTT7_ZETAS, NTT7_ZETAS_INV  twiddle tables (Montgomery form, index k = number of blocks + block)
 *   NTT7_NINV                 128^-1 in Montgomery form
 *   NTT7_FQMUL(a, b)          scalar Montgomery product (to merge the scaling into a twiddle)
 *   NTT7_RED_IN(x)            reduction of the forward input (identity if 7 layers cannot overflow)
 *   NTT7_FWD_OUT(x)           final reduction of the forward transform
 *   NTT7_RED_SUM(x)           reduction of a Gentleman-Sande sum
 *   NTT7_RED_MID              1 if the large sums must be reduced between the two merged groups
 *   NTT7_INV_OUT(x)           final reduction of the inverse transform
 *
 * Three layers are merged on eight vectors held in registers: the 512 (1024) coefficients are read and written
 * three times per transform instead of once per layer.
 *   forward: layers 1-3 on vectors 8S apart (S = END_N/512), layers 4-6 inside blocks of 8S vectors,
 *            layer 7 together with the final reduction;   inverse: the mirror image, scaling merged into the end.
 * Range discipline (Montgomery products are in (-p, p), every Cooley-Tukey layer adds less than p):
 *   forward : |in| after NTT7_RED_IN plus 7p must stay below 2^15 (3329: 1665 + 7p = 24968; 257: 16384 + 1799).
 *   inverse : |in| <= 4095, so the sums and differences of layer 7 are below 8190. Layer 7 reduces its sums with
 *             NTT7_RED_SUM; its products are in (-p, p). Three merged Gentleman-Sande layers multiply a bound by
 *             at most 8:
 *               3329: values <= 3328 before layers 6-4, sums <= 26624 after them; NTT7_RED_MID reduces these
 *                     four sums again, so layers 3-1 also start from <= 3328;
 *               257 : values <= 383 before layers 6-4, <= 3064 after them, <= 24512 after layers 3-1.
 *             The final multiplication by 1/128 accepts any int16.
 */

#define NTT7_NV (END_N / 8)             /* vectors per polynomial            */
#define NTT7_S  (END_N / 512)           /* vectors per leaf pair: 1 or 2     */

/* twiddles of a merged group whose coarsest layer uses index k0: z[0] | z[1], z[2] | z[3..6] */
static inline void NTT7_CAT(tw_fwd_, NTT7_P)(int16x8_t z[7], int16x8_t zq[7], unsigned k0)
{
    const int16x8_t pinv = vdupq_n_s16(NTT7_PINV);
    unsigned i;
    z[0] = vdupq_n_s16(NTT7_ZETAS[k0]);
    z[1] = vdupq_n_s16(NTT7_ZETAS[2 * k0]);
    z[2] = vdupq_n_s16(NTT7_ZETAS[2 * k0 + 1]);
    for (i = 0; i < 4; i++) z[3 + i] = vdupq_n_s16(NTT7_ZETAS[4 * k0 + i]);
    for (i = 0; i < 7; i++) zq[i] = vmulq_s16(z[i], pinv);
}

/* inverse twiddles of the same group, finest layer first: z[0..3] | z[4], z[5] | z[6] */
static inline void NTT7_CAT(tw_inv_, NTT7_P)(int16x8_t z[7], int16x8_t zq[7], unsigned k0)
{
    const int16x8_t pinv = vdupq_n_s16(NTT7_PINV);
    unsigned i;
    for (i = 0; i < 4; i++) z[i] = vdupq_n_s16(NTT7_ZETAS_INV[4 * k0 + i]);
    z[4] = vdupq_n_s16(NTT7_ZETAS_INV[2 * k0]);
    z[5] = vdupq_n_s16(NTT7_ZETAS_INV[2 * k0 + 1]);
    z[6] = vdupq_n_s16(NTT7_ZETAS_INV[k0]);
    for (i = 0; i < 7; i++) zq[i] = vmulq_s16(z[i], pinv);
}

void NTT7_FWD(int16_t a[END_N])
{
    const int16x8_t pv = vdupq_n_s16(NTT7_P), pinv = vdupq_n_s16(NTT7_PINV);
    int16x8_t v[8], z[7], zq[7];
    unsigned i, j, blk;
    /* layers 1-3: eight vectors 8S vectors apart */
    NTT7_CAT(tw_fwd_, NTT7_P)(z, zq, 1);
    for (j = 0; j < 8 * NTT7_S; j++) {
        int16_t *p = a + 8 * j;
        for (i = 0; i < 8; i++) v[i] = NTT7_RED_IN(vld1q_s16(p + 64 * NTT7_S * i));
        neon_ct3(v, z, zq, pv);
        for (i = 0; i < 8; i++) vst1q_s16(p + 64 * NTT7_S * i, v[i]);
    }
    /* layers 4-6: inside the eight blocks of 8S vectors, eight vectors S vectors apart */
    for (blk = 0; blk < 8; blk++) {
        NTT7_CAT(tw_fwd_, NTT7_P)(z, zq, 8 + blk);
        for (j = 0; j < NTT7_S; j++) {
            int16_t *p = a + 64 * NTT7_S * blk + 8 * j;
            for (i = 0; i < 8; i++) v[i] = vld1q_s16(p + 8 * NTT7_S * i);
            neon_ct3(v, z, zq, pv);
            for (i = 0; i < 8; i++) vst1q_s16(p + 8 * NTT7_S * i, v[i]);
        }
    }
    /* layer 7 and the final reduction, two vectors at a time */
    for (i = 0; i < END_N; i += 16) {
        int16x8_t v0 = vld1q_s16(a + i), v1 = vld1q_s16(a + i + 8), t, zz, zzq;
#if END_LEAF == 4
        /* leaves of 4: the butterflies pair the two halves of a vector; one twiddle per vector */
        int16x8_t lo = neon_zip64_lo(v0, v1), hi = neon_zip64_hi(v0, v1);
        zz = vcombine_s16(vdup_n_s16(NTT7_ZETAS[64 + i / 8]), vdup_n_s16(NTT7_ZETAS[65 + i / 8]));
        zzq = vmulq_s16(zz, pinv);
        t = neon_montmul(hi, zz, zzq, pv);
        hi = vsubq_s16(lo, t);
        lo = vaddq_s16(lo, t);
        v0 = neon_zip64_lo(lo, hi);
        v1 = neon_zip64_hi(lo, hi);
#else
        /* leaves of 8: the butterflies pair two adjacent vectors */
        zz = vdupq_n_s16(NTT7_ZETAS[64 + i / 16]);
        zzq = vmulq_s16(zz, pinv);
        t = neon_montmul(v1, zz, zzq, pv);
        v1 = vsubq_s16(v0, t);
        v0 = vaddq_s16(v0, t);
#endif
        vst1q_s16(a + i, NTT7_FWD_OUT(v0));
        vst1q_s16(a + i + 8, NTT7_FWD_OUT(v1));
    }
}

void NTT7_INV(int16_t a[END_N])
{
    const int16x8_t pv = vdupq_n_s16(NTT7_P), pinv = vdupq_n_s16(NTT7_PINV);
    int16x8_t v[8], z[7], zq[7], zn, znq;
    unsigned i, j, blk;
    /* layer 7 with reduced sums: every value is then below max(p, bound of NTT7_RED_SUM) */
    for (i = 0; i < END_N; i += 16) {
        int16x8_t v0 = vld1q_s16(a + i), v1 = vld1q_s16(a + i + 8), d, zz, zzq;
#if END_LEAF == 4
        int16x8_t lo = neon_zip64_lo(v0, v1), hi = neon_zip64_hi(v0, v1);
        zz = vcombine_s16(vdup_n_s16(NTT7_ZETAS_INV[64 + i / 8]), vdup_n_s16(NTT7_ZETAS_INV[65 + i / 8]));
        zzq = vmulq_s16(zz, pinv);
        d = vsubq_s16(lo, hi);
        lo = NTT7_RED_SUM(vaddq_s16(lo, hi));
        hi = neon_montmul(d, zz, zzq, pv);
        v0 = neon_zip64_lo(lo, hi);
        v1 = neon_zip64_hi(lo, hi);
#else
        zz = vdupq_n_s16(NTT7_ZETAS_INV[64 + i / 16]);
        zzq = vmulq_s16(zz, pinv);
        d = vsubq_s16(v0, v1);
        v0 = NTT7_RED_SUM(vaddq_s16(v0, v1));
        v1 = neon_montmul(d, zz, zzq, pv);
#endif
        vst1q_s16(a + i, v0);
        vst1q_s16(a + i + 8, v1);
    }
    /* layers 6-4: three layers multiply a bound by at most 8 */
    for (blk = 0; blk < 8; blk++) {
        NTT7_CAT(tw_inv_, NTT7_P)(z, zq, 8 + blk);
        for (j = 0; j < NTT7_S; j++) {
            int16_t *p = a + 64 * NTT7_S * blk + 8 * j;
            for (i = 0; i < 8; i++) v[i] = vld1q_s16(p + 8 * NTT7_S * i);
            neon_gs3(v, z, zq, pv);
#if NTT7_RED_MID
            for (i = 0; i < 4; i++) v[i] = NTT7_RED_SUM(v[i]);       /* the sums; v[4..7] are products in (-p, p) */
#endif
            for (i = 0; i < 8; i++) vst1q_s16(p + 8 * NTT7_S * i, v[i]);
        }
    }
    /* layers 3-1 and the scaling by 1/128: merged into the last twiddle for the products, a multiplication for the sums */
    NTT7_CAT(tw_inv_, NTT7_P)(z, zq, 1);
    z[6] = vdupq_n_s16(NTT7_FQMUL(NTT7_ZETAS_INV[1], NTT7_NINV));
    zq[6] = vmulq_s16(z[6], pinv);
    zn = vdupq_n_s16(NTT7_NINV);
    znq = vmulq_s16(zn, pinv);
    for (j = 0; j < 8 * NTT7_S; j++) {
        int16_t *p = a + 8 * j;
        for (i = 0; i < 8; i++) v[i] = vld1q_s16(p + 64 * NTT7_S * i);
        neon_gs3(v, z, zq, pv);
        for (i = 0; i < 4; i++) v[i] = neon_montmul(v[i], zn, znq, pv);
        for (i = 0; i < 8; i++) vst1q_s16(p + 64 * NTT7_S * i, NTT7_INV_OUT(v[i]));
    }
}

#undef NTT7_NV
#undef NTT7_S
