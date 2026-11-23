/*
 * END KEM -- Cortex-M4 kernels: inner loops of the trapdoor decoder and the constant-time vector helpers
 * (Thumb-2 + DSP).  They replace the portable versions of src/common/end_deckern.c and of polyvec.c
 * (see end_backend.h) and compute exactly the same values.
 *
 *   poly_crt_c2        reconstruct the centred CRT coefficient using integer DSP only
 *   poly_cand_round    rounding vector and residue of a gf-stage candidate, 32-bit arithmetic only
 *   poly_cand_combine  s~'' of a candidate and its squared norm, two coefficients per instruction
 *   poly_cand_check    GF-stage candidate and its "all entries are +-1" test
 *   poly_lift257       parity lift of a residue mod 257
 *   poly_cmov16/32     masked copy with SEL
 *   poly_rot16         negacyclic rotation by a secret amount: barrel shifter, one SEL per word and stage
 *
 * Constant time: the candidate signs s0, s1 are public (the candidate list is fixed); every data value and the
 * rotation amount are secret. There are no secret-dependent branches or addresses; selections use masks or the
 * GE flags (UADD16 of an all-ones / all-zero mask sets or clears all four) with SEL. s0 is used as scratch.
 */
#include "asm_defs.h"

    .syntax unified
    .thumb

#define NBYTES   (2*END_N)                  /* bytes of an int16 polynomial */

#if defined(__ARM_FEATURE_DSP) && !defined(END_M_PORTABLE_CRT)
    .text
    .align 2

/* C2 = Q*k - t12. Since Q^-1 = -175 mod 514, with w the constant
 * coefficient and p the desired parity of u, we may take
 * k = c2 - 175*(t12 - 2*t257 + w) + 257*((w ^ p) & 1) mod 514.
 * Indeed 2*t257 - w + 257*((w ^ p) & 1) lifts u mod 257 to parity p.
 * Choose k in [-257,256] when t12 <= 0, [-256,257] otherwise. */
.macro crt_lane high
    .if \high
    asr     lr, r9, #16                   /* signed t12 */
    asr     r11, r10, #16                /* c2 plus parity correction */
    smlatb  r11, r8, r6, r11             /* k before centring */
    .else
    sxth    lr, r9
    sxth    r11, r10
    smlabb  r11, r8, r6, r11
    .endif
    /* Round (2*k +/- 1)/1028 to the nearest integer. The sign resolves
     * the endpoint at k = 257 mod 514. The numerator is odd and below
     * 2^22 in magnitude, so rounding the reciprocal cannot cross a tie. */
    sub     r12, lr, #1
    sub     r12, r11, r12, asr #31       /* k + (t12 <= 0) */
    lsl     r12, r12, #1
    sub     r12, r12, #1
    smmulr  r12, r12, r7
    add     r12, r12, r12, lsl #8
    sub     r11, r11, r12, lsl #1        /* subtract 514 times the quotient */
    add     r12, r11, r11, lsl #1
    add     r11, r11, r12, lsl #12      /* 12289*k */
    sub     r11, r11, lr
    str     r11, [r0], #4
.endm

    .global END_SYM(poly_crt_c2)
    .type END_SYM(poly_crt_c2), %function
END_SYM(poly_crt_c2):
    push    {r4-r11, lr}
    ldr     r4, [sp, #36]                /* wconst */
    ldr     r5, [sp, #40]                /* pwc1 */
    orr     r5, r5, r5, lsl #16
    mvn     r6, #174                     /* -175 = 12289^-1 mod 514 */
    movw    r7, #0xc040
    movt    r7, #0x003f                  /* round(2^32 / 1028) */
    add     r8, r0, #(4*END_N)
    push    {r8}
.Lcrt_pair:
    ldr     r8, [r1], #4                 /* t257 */
    ldr     r9, [r2], #4                 /* t12289 */
    ldr     r10, [r4], #4                /* wconst */
    ldr     r11, [r3], #4                /* c2 */
    sadd16  r8, r8, r8
    ssub16  r8, r9, r8
    sadd16  r8, r8, r10                 /* signed difference, fits a halfword */
    eor     r10, r10, r5
    and     r10, r10, #0x00010001
    add     r10, r10, r10, lsl #8
    sadd16  r10, r10, r11
    crt_lane 0
    crt_lane 1
    ldr     r8, [sp]
    cmp     r0, r8
    bne     .Lcrt_pair
    add     sp, sp, #4
    pop     {r4-r11, pc}
    .size END_SYM(poly_crt_c2), .-END_SYM(poly_crt_c2)
#endif

#ifdef END_M4_ASM
#define Q        257
#define QQP      3158273                    /* q Q                          */
#define K2QQP    6316546                    /* 2 q Q                        */
#define ROFF     (QQP+32*K2QQP)             /* offset that makes the dividend positive */

    .text

/* ------------------------------------------------------------------------------------------------------------
 * void poly_cand_round(int16_t r[N], int32_t c2t[N], const int32_t C2[N], const int16_t w0[N], const int16_t w1[N],
 *                      int32_t s0, int32_t s1);
 *   v = C2 - 2q (s0 w0 + s1 w1);  x = v + qQ + 32*2qQ in (0, 2^29);  q0 = (x * 679) >> 32 = floor(x / 2qQ) or one less
 *   rem = x - q0*2qQ - 2qQ;  r = q0 + (rem < 0 ? 0 : 1) - 32;  c2t = (rem < 0 ? rem + 2qQ : rem) - qQ
 *   (the same computation as round_div_2qqp_rem in end_arith.h)
 * ------------------------------------------------------------------------------------------------------------ */
.macro round_tail
    add     r10, r10, r7
    smmul   r11, r10, r8
    mls     r10, r11, r9, r10
    sub     r10, r10, r9
    add     r11, r11, r10, asr #31
    sub     r11, r11, #31
    strh    r11, [r0], #2
    and     r11, r9, r10, asr #31
    add     r10, r10, r11
    sub     r10, r10, r9, lsr #1
    str.w   r10, [r1], #4
.endm

    .global END_SYM(poly_cand_round)
    .type   END_SYM(poly_cand_round), %function
    .align  2
END_SYM(poly_cand_round):
    push    {r4-r11, lr}
    ldr     r4, [sp, #36]               /* w1 */
    ldr     r5, [sp, #40]               /* s0 */
    ldr     r6, [sp, #44]               /* s1 */
    movw    r7, #(2*Q)
    mul     r5, r5, r7
    mul     r6, r6, r7
    rsb     r5, r5, #0                  /* -2q s0 */
    rsb     r6, r6, #0                  /* -2q s1 */
    movw    r7, #(ROFF & 0xFFFF)
    movt    r7, #(ROFF >> 16)
    movw    r8, #679                    /* floor(2^32 / 2qQ) */
    movw    r9, #(K2QQP & 0xFFFF)
    movt    r9, #(K2QQP >> 16)
    add.w   r10, r0, #NBYTES
    vmov    s0, r10
1:
    .rept   2
    ldr.w   r12, [r3], #4               /* two coefficients of w0 */
    ldr.w   r14, [r4], #4               /* two coefficients of w1 */
    ldr.w   r10, [r2], #4
    smlabb  r10, r12, r5, r10
    smlabb  r10, r14, r6, r10
    round_tail
    ldr.w   r10, [r2], #4
    smlatb  r10, r12, r5, r10
    smlatb  r10, r14, r6, r10
    round_tail
    .endr
    vmov    r10, s0
    cmp     r0, r10
    bne.w   1b
    pop     {r4-r11, pc}
    .size   END_SYM(poly_cand_round), .-END_SYM(poly_cand_round)

/* ------------------------------------------------------------------------------------------------------------
 * uint64_t poly_cand_combine(int16_t t[N], const int16_t sp[N], const int16_t F0[N], const int16_t F1[N],
 *                            const int16_t rf[N], int32_t s0, int32_t s1);
 *   t = sp + 2 (s0 F0 + s1 F1) + 2 rf, returns the sum of the squares (SMLALD)
 * ------------------------------------------------------------------------------------------------------------ */
    .global END_SYM(poly_cand_combine)
    .type   END_SYM(poly_cand_combine), %function
    .align  2
END_SYM(poly_cand_combine):
    push    {r4-r11, lr}
    ldr     r4, [sp, #36]               /* rf */
    ldr     r5, [sp, #40]               /* s0 */
    ldr     r6, [sp, #44]               /* s1 */
    lsl     r5, r5, #1
    lsl     r6, r6, #1
    movs    r7, #0
    movs    r8, #0
    add.w   r9, r0, #NBYTES
    vmov    s0, r9
1:
    .rept   2
    ldr.w   r9, [r2], #4
    ldr.w   r10, [r3], #4
    ldr.w   r11, [r4], #4
    ldr.w   r12, [r1], #4
    smulbb  r14, r9, r5
    smlabb  r14, r10, r6, r14
    smultb  r9, r9, r5
    smlatb  r9, r10, r6, r9
    pkhbt   r14, r14, r9, lsl #16       /* 2 (s0 F0 + s1 F1) */
    sadd16  r11, r11, r11
    sadd16  r14, r14, r11
    sadd16  r14, r14, r12
    smlald  r7, r8, r14, r14
    str.w   r14, [r0], #4
    .endr
    vmov    r9, s0
    cmp     r0, r9
    bne.w   1b
    mov     r0, r7
    mov     r1, r8
    pop     {r4-r11, pc}
    .size   END_SYM(poly_cand_combine), .-END_SYM(poly_cand_combine)

/* ------------------------------------------------------------------------------------------------------------
 * uint32_t poly_cand_check(int16_t t[N], const int16_t b[N], const int16_t f0[N], const int16_t f1[N],
 *                          int32_t s0, int32_t s1);
 *   t = b - 2 (s0 f0 + s1 f1); a half x is +-1 exactly when (x + 1) & ~2 == 0; returns the OR of these values
 * ------------------------------------------------------------------------------------------------------------ */
    .global END_SYM(poly_cand_check)
    .type   END_SYM(poly_cand_check), %function
    .align  2
END_SYM(poly_cand_check):
    push    {r4-r11, lr}
    ldr     r4, [sp, #36]               /* s0 */
    ldr     r5, [sp, #40]               /* s1 */
    lsl     r4, r4, #1
    lsl     r5, r5, #1
    mov.w   r6, #0x00010001
    movs    r7, #0
    add.w   r12, r0, #NBYTES
1:
    .rept   2
    ldr.w   r8, [r2], #4
    ldr.w   r9, [r3], #4
    ldr.w   r10, [r1], #4
    smulbb  r11, r8, r4
    smlabb  r11, r9, r5, r11
    smultb  r8, r8, r4
    smlatb  r8, r9, r5, r8
    pkhbt   r11, r11, r8, lsl #16       /* 2 (s0 f0 + s1 f1) */
    ssub16  r11, r10, r11
    str.w   r11, [r0], #4
    sadd16  r11, r11, r6
    bic     r11, r11, #0x00020002
    orr     r7, r7, r11
    .endr
    cmp     r0, r12
    bne.w   1b
    mov     r0, r7
    pop     {r4-r11, pc}
    .size   END_SYM(poly_cand_check), .-END_SYM(poly_cand_check)

/* ------------------------------------------------------------------------------------------------------------
 * void poly_lift257(int16_t out[N], const int16_t t[N], const int16_t cst[N], uint16_t parity);
 *   x = (2 t - cst) mod 257 in [0, 257)   (add 257, fold the two bytes, add 257 where negative)
 *   x += 257 where the parity of x differs from `parity`;  x -= 514 where x >= 257
 * ------------------------------------------------------------------------------------------------------------ */
    .global END_SYM(poly_lift257)
    .type   END_SYM(poly_lift257), %function
    .align  2
END_SYM(poly_lift257):
    push    {r4-r8}
    orr     r3, r3, r3, lsl #16
    mov.w   r4, #0x01010101             /* 257 twice */
    lsl     r5, r4, #1                  /* 514 twice */
    add.w   r12, r0, #NBYTES
1:
    ldr.w   r6, [r1], #4
    ldr.w   r7, [r2], #4
    sadd16  r6, r6, r6
    ssub16  r6, r6, r7                  /* 2t - cst in [-256, 512] */
    uadd16  r6, r6, r4                  /* [1, 769] */
    uxtb16  r7, r6
    uxtb16  r6, r6, ror #8
    ssub16  r6, r7, r6                  /* [-3, 255] */
    uadd16  r7, r6, r4
    sel     r6, r7, r6                  /* [0, 257) */
    eor     r7, r6, r3
    and     r7, r7, #0x00010001
    add     r7, r7, r7, lsl #8          /* 257 where the parity is wrong */
    add     r6, r6, r7                  /* [0, 514), no carry between the halves */
    ssub16  r7, r6, r5
    ssub16  r8, r6, r4                  /* GE: x >= 257 */
    sel     r6, r7, r6
    str.w   r6, [r0], #4
    cmp     r0, r12
    bne.w   1b
    pop     {r4-r8}
    bx      lr
    .size   END_SYM(poly_lift257), .-END_SYM(poly_lift257)

/* ------------------------------------------------------------------------------------------------------------
 * void poly_cmov16(int16_t r[N], const int16_t a[N], uint16_t mask);     r = mask ? a : r   (mask 0 or 0xFFFF)
 * void poly_cmov32(int32_t r[N], const int32_t a[N], uint32_t mask);     (mask 0 or 0xFFFFFFFF)
 * ------------------------------------------------------------------------------------------------------------ */
    .global END_SYM(poly_cmov16)
    .type   END_SYM(poly_cmov16), %function
    .global END_SYM(poly_cmov32)
    .type   END_SYM(poly_cmov32), %function
    .align  2
END_SYM(poly_cmov16):
    orr     r2, r2, r2, lsl #16
    add.w   r3, r0, #NBYTES
    b       2f
END_SYM(poly_cmov32):
    add.w   r3, r0, #(2*NBYTES)
2:
    push    {r4-r11}
    uadd16  r12, r2, r2                 /* GE = 1111 iff the mask is all-ones */
1:
    ldr.w   r4, [r1], #4
    ldr.w   r5, [r1], #4
    ldr.w   r6, [r1], #4
    ldr.w   r7, [r1], #4
    ldr.w   r8, [r0]
    ldr.w   r9, [r0, #4]
    ldr.w   r10, [r0, #8]
    ldr.w   r11, [r0, #12]
    sel     r8, r4, r8
    sel     r9, r5, r9
    sel     r10, r6, r10
    sel     r11, r7, r11
    str.w   r8, [r0], #4
    str.w   r9, [r0], #4
    str.w   r10, [r0], #4
    str.w   r11, [r0], #4
    cmp     r0, r3
    bne.w   1b
    pop     {r4-r11}
    bx      lr
    .size   END_SYM(poly_cmov16), .-END_SYM(poly_cmov16)
    .size   END_SYM(poly_cmov32), .-END_SYM(poly_cmov32)

/* ------------------------------------------------------------------------------------------------------------
 * void poly_rot16(int16_t r[N], const int16_t a[N], uint32_t k);     r = x^k a (negacyclic), k secret, r != a
 *
 * Barrel shifter with log2(N) stages. Stage 0 shifts by one coefficient (bit 0 of k), stage b >= 1 by 2^(b-1)
 * words (bit b). Every stage reads one buffer and writes another: the buffers alternate between r and a stack
 * frame of one polynomial, starting so that the last stage writes r.
 *   stage b >= 1, shift of sw words:   dst[i] = bit ? -src[i + NW - sw] : src[i]      for i <  sw  (masks)
 *                                      dst[i] = bit ?  src[i - sw]      : src[i]      for i >= sw  (SEL)
 * ------------------------------------------------------------------------------------------------------------ */
    .global END_SYM(poly_rot16)
    .type   END_SYM(poly_rot16), %function
    .align  2
END_SYM(poly_rot16):
    push    {r4-r11, lr}
    sub.w   sp, sp, #NBYTES
#if (END_N == 512)
    mov     r3, r0                      /* 9 stages: stage 0 writes r */
    mov     r4, sp
#else
    mov     r3, sp                      /* 10 stages: stage 0 writes the frame */
    mov     r4, r0
#endif
    /* ---- stage 0: one coefficient ---- */
    and     r5, r2, #1
    rsb     r5, r5, #0
    uadd16  r5, r5, r5                  /* GE from bit 0 */
    ldr.w   r6, [r1, #(NBYTES-4)]
    asr     r6, r6, #16
    rsb     r6, r6, #0
    uxth    r6, r6                      /* -a[N-1]: the coefficient that wraps around */
    mov     r7, r3
    add.w   r8, r1, #NBYTES
1:
    .rept   4
    ldr.w   r9, [r1], #4
    orr     r10, r6, r9, lsl #16        /* [previous high half | low half] */
    sel     r10, r10, r9
    str.w   r10, [r7], #4
    lsr     r6, r9, #16
    .endr
    cmp     r1, r8
    bne.w   1b
    /* ---- stages 1 .. log2(N)-1: r3 = source, r4 = destination, r5 = shift in bytes, r2 = remaining bits of k ---- */
    lsr     r2, r2, #1
    movs    r5, #4
2:
    and     r6, r2, #1
    rsb     r6, r6, #0                  /* mask of this stage */
    add.w   r7, r3, #NBYTES
    sub     r7, r7, r5                  /* source of the wrapped words */
    mov     r8, r3
    mov     r9, r4
    add     r10, r4, r5
    movs    r11, #0
3:
    ldr.w   r12, [r7], #4
    ldr.w   r14, [r8], #4
    ssub16  r12, r11, r12               /* negate both halves */
    and     r12, r12, r6
    bic     r14, r14, r6
    orr     r12, r12, r14
    str.w   r12, [r9], #4
    cmp     r9, r10
    bne.w   3b
    uadd16  r12, r6, r6                 /* GE from the mask */
    mov     r7, r3
    rsb     r10, r5, #NBYTES            /* bytes left */
    and     r10, r10, #28               /* (words left mod 8) * 4 */
    add     r10, r9, r10
    cmp     r9, r10
    beq.w   5f
4:
    ldr.w   r12, [r7], #4
    ldr.w   r14, [r8], #4
    sel     r14, r12, r14
    str.w   r14, [r9], #4
    cmp     r9, r10
    bne.w   4b
5:
    add.w   r10, r4, #NBYTES
    cmp     r9, r10
    beq.w   7f
    vmov    s0, r2
    vmov    s1, r5
6:
    .rept   2
    ldr.w   r0, [r7], #4
    ldr.w   r1, [r7], #4
    ldr.w   r2, [r7], #4
    ldr.w   r5, [r7], #4
    ldr.w   r6, [r8], #4
    ldr.w   r11, [r8], #4
    ldr.w   r12, [r8], #4
    ldr.w   r14, [r8], #4
    sel     r6, r0, r6
    sel     r11, r1, r11
    sel     r12, r2, r12
    sel     r14, r5, r14
    str.w   r6, [r9], #4
    str.w   r11, [r9], #4
    str.w   r12, [r9], #4
    str.w   r14, [r9], #4
    .endr
    cmp     r9, r10
    bne.w   6b
    vmov    r2, s0
    vmov    r5, s1
7:
    mov     r12, r3                     /* swap the buffers */
    mov     r3, r4
    mov     r4, r12
    lsr     r2, r2, #1
    lsl     r5, r5, #1
    cmp     r5, #NBYTES
    bne.w   2b
    add.w   sp, sp, #NBYTES
    pop     {r4-r11, pc}
    .size   END_SYM(poly_rot16), .-END_SYM(poly_rot16)

#endif /* END_M4_ASM */
